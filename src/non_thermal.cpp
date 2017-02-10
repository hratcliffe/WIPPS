

//  Created by Heather Ratcliffe on 11/02/2016.
/** Small class to hold a non-thermal electron distribution we can operate on. Mainly contains routine to parse from deck.status using nonthermal.conf defined params or a generic assumed deck*/

#include <stdio.h>
#include <cmath>
#include <boost/math/special_functions.hpp>
//Provides Bessel functions, erf, and many more
#include <functional>
#include <boost/math/special_functions/sign.hpp>
#include "support.h"
#include "non_thermal.h"
#include "data_array.h"

extern deck_constants my_const;
extern const mpi_info_struc mpi_info;


non_thermal::non_thermal(std::string file_prefix){
/** \brief Construct non-thermal distrib
*
*Sets default params
*/
  ref_dens = 1.0;
  ref_B = 1.0;
  fraction = 0.0;
  v_par = 0.0;
  v_perp = 0.0;
  norm = 1.0;
  ncomps = 1;
  total_dens = 0;
  norely = 0;
  //Defaults which are deliberately meaningless....
  
  dp = 10.0;
  //Sensible dp for electrons (m/s)
  
  bool err = configure_from_file(file_prefix);
  if(err) my_print("WARNING: Error getting electron data from .status file", mpi_info.rank);

}

bool non_thermal::configure_from_file(std::string file_prefix){
/** \brief Setup non-thermal distribution from file
*
*Reads a [prefix]deck.status file and parses out the neccessary parameters. Note that this depends on the specific deck format used to get the parameter names \todo Multithermal
*/

  this->total_dens = 1.0;

  std::ifstream infile;
  infile.open(file_prefix+"deck.status");
  std::string line;
  std::vector<std::string> lines;
  bool found=false;
  if(!infile.is_open()){
    my_print("No deck.status file found, aborting", mpi_info.rank);
    safe_exit();
  
  }
  while(infile){
    getline(infile, line);
    if(!(strcmp(line.substr(0, CONSTANTS.size()).c_str(), CONSTANTS.c_str()))){
      found = true;
      break;
    }
  }
  if(!found){
    my_print("No constants dump found in deck.status, aborting", mpi_info.rank);
    safe_exit();
  }

  while(infile){
    getline(infile, line);
    if(!(strcmp(line.substr(0, CONSTANTS_END.size()).c_str(), CONSTANTS_END.c_str()))) break;
    lines.push_back(line);
  }
  
  std::string name, val;
  bool parse_err;
  float val_f;
  //parse all params into map
  for(size_t i=0; i< lines.size(); i++){

    parse_err = parse_name_val(lines[i], name, val);
    if(parse_err) continue;
    if(name == "lookup"){
    //Any stringy params here
      extra_parameters[name] = val;
    }else{
      val_f = atof(val.c_str());
      
      parameters[name] = val_f;
    }
  }
  infile.close();
  /*If electron setup file found, identify components. File should look like this:
  n_comps = 2
  warm:
  function = max
  dens = dens_rat
  vpar = vpar
  hot:
  function = max
  dens=dens_rath
  etc where the RHS's are the string names in the deck and the function is either maxwell, kappa etc. ALL LOWER CASE. this then defines a plasma from deck constants*/
  /** \todo And in plasma, fix lower casing*/
  
  my_print("Reading "+file_prefix+"nonthermal.conf");
  infile.open(file_prefix+"nonthermal.conf");
  if(!infile){
    my_print("Using generic bimaxwellian electrons");

    auto tmp_fn = std::bind(bimax, std::placeholders::_1, std::placeholders::_2, parameters[VPAR], parameters[VPERP], parameters[DENS_RAT]);
    f_p_private.push_back(tmp_fn);
    total_dens += parameters[DENS_RAT];

  }else{
    
    std::string line, name, val;
    int block_num = -1;
    bool parse_err;

    //very naive parsing. We spin through until we find a ":"
    //if we don't find n_comps such blocks, we report and continue
    std::string function="", dens="", vpar="", vperp="", kappa="", lookup="";
    std::function<calc_type(calc_type, calc_type)> tmp_fn;
    bool function_set = false;

    while(getline(infile, line)){
      function_set= false;
      if(line.find(':') != std::string::npos){
        //Create a plasma component function bound to these parameters and add to vector
        calc_type norm;
        if(function !=""){
          //Some general setups
          
          if(parameters.count(dens)!=0) total_dens += parameters[dens];
        }
        if(function == "max"){
          //Check we have all params
          my_print("Adding component "+dens);
          for(auto nam : {dens, vpar, vperp}) if(!parameters.count(nam)) my_print("Param "+nam+" not found!!!");
          this->v_par = parameters[vpar];
          this->v_perp = parameters[vperp];
          calc_type a_par = std::sqrt(2.0)*parameters[vpar];
          if(!this->norely) a_par/= std::sqrt(1.0 - std::pow(parameters[vpar]/v0, 2));
          
          calc_type a_perp_sq = 2.0*std::pow(parameters[vperp], 2);
          if(!this->norely) a_perp_sq /= (1.0 - 2.0*std::pow(parameters[vperp]/v0, 2));
          norm = 1.0/(pi*std::sqrt(pi)*a_par*a_perp_sq);
          tmp_fn = std::bind(bimax, std::placeholders::_1, std::placeholders::_2, a_par, std::sqrt(a_perp_sq), norm*parameters[dens]);
          f_p_private.push_back(tmp_fn);
          function_set = true;
        }
        else if(function =="kappa"){
          for(auto nam : {dens, vpar, vperp, kappa}) if(!parameters.count(nam)) my_print("Param "+nam+" not found!!!");

          tmp_fn = std::bind(bikappa, std::placeholders::_1, std::placeholders::_2, parameters[kappa], parameters[vpar], parameters[vperp], parameters[dens]);
          f_p_private.push_back(tmp_fn);
          function_set = true;

        }else if(function=="lookup"){
          tmp_fn = configure_lookup(file_prefix, lookup, this);
          f_p_private.push_back(tmp_fn);
          function_set = true;
        }
        if(function_set) my_print("Configured function "+function);
        block_num ++;
        function=""; dens=""; vpar=""; vperp="";kappa="", lookup="";
        if(block_num >= ncomps) break;
        //Reset lookups
        continue;
        // is next block, skip this header line
      }
      if(block_num == -1 && line.find('=')){
      //This might be an n_comps definition or nonrely flag!
        parse_err = parse_name_val(line, name, val);
        if(!parse_err){
          if(name =="ncomps") this->ncomps = atoi(val.c_str());
          else if(name =="nonrely"){
            if(val.size() >0 && val[0] == '1'){
              this->norely=1;
              my_print("Using nonrelativistic calculation!");
            }
          }else{
            my_print("!!!!!Unknown name "+name);
          }
          
        }
      }
      if(block_num >= 0){
        //this line is a valid input one, probably!
        parse_err = parse_name_val(line, name, val);
        if(!parse_err){
          //Grab the strings for each parameter
          if(name == "function") function = val;
          else if(name =="dens") dens = val;
          else if(name =="vpar") vpar = val;
          else if(name =="vperp") vperp = val;
          else if(name =="kappa") kappa = val;
          else if(name =="lookup") lookup = val;
          else my_print("!!!!Unknown name "+name);
        }
      }
      
    }
    if(block_num > ncomps){
    my_print("Too many blocks in plasma file, ignoring", mpi_info.rank);
    }
    else if(block_num < ncomps-1){
      my_print(mk_str(block_num), mpi_info.rank);
  
      my_print("Insufficient blocks in config file", mpi_info.rank);
    }
  }
  return 0;

}

void non_thermal::write(std::ofstream &outfile){
/**\brief Record non-thermal params
*
*Write params to file for e.g. log, or source for a redo
*/
  outfile<<DENS<<" "<<ref_dens<<"\n";
  outfile<<"B_ref "<<ref_B<<"\n";
  outfile<<DENS_RAT<<" "<<fraction<<"\n";
  outfile<<VPAR<<" "<<v_par<<"\n";
  outfile<<VPERP<<" "<<v_perp<<"\n";

}

void non_thermal::dump(std::fstream &outfile){
/**\brief Record non-thermal params
*
*Dump without tags etc
*/
  outfile.write((char*)& ref_dens, sizeof(calc_type));
  outfile.write((char*)& fraction, sizeof(calc_type));
  outfile.write((char*)& v_par, sizeof(calc_type));
  outfile.write((char*)& v_perp, sizeof(calc_type));
  
}

calc_type non_thermal::d_f_p(calc_type p_par, calc_type p_perp, bool parallel){
/** Simple Numerical derivative*/
  
  if(parallel) return (f_p(p_par+dp, p_perp) -f_p(p_par, p_perp))/dp;
  else return (f_p(p_par, p_perp+dp) -f_p(p_par, p_perp))/dp;
  
}

calc_type non_thermal::f_p(calc_type p_par, calc_type p_perp){
/**\brief Return f(p)
*
*Evaluates function at p_perp, p_par and returns result. Calls bound function so we can use multiple functional forms for f
*/

  calc_type ret = 0;
  for(size_t i=0; i<f_p_private.size(); i++) ret +=f_p_private[i](p_par, p_perp);
  return ret;

}


//The following are the generic functions which we can bind. p is always the first parameter

calc_type max(calc_type p, calc_type p_th, calc_type A){
/** Exp function*/

  return A*std::exp(-p*p/p_th/p_th);
}

calc_type bimax(calc_type p, calc_type p2, calc_type p_th, calc_type p_th2, calc_type A){
  return A*std::exp(-p*p/p_th/p_th - p2*p2/p_th2/p_th2);
}

calc_type double_max(calc_type p, calc_type v_th, calc_type v_th2, calc_type A, calc_type A2){
/** Double exp function*/

  return A*std::exp(-p*p/v_th/v_th/me/me) + A2*std::exp(-p*p/v_th2/v_th2/me/me);
}

calc_type bikappa(calc_type p, calc_type p2, calc_type kappa, calc_type v_kpar, calc_type v_kperp, calc_type A){

  return 0;
}

calc_type pancake(calc_type p, calc_type p2, size_t px_sz, calc_type * px_vals){

  return 0.0;
}

calc_type lookup(calc_type p_par, calc_type p_perp, my_type * data, size_t par_sz, size_t perp_sz,calc_type dp_par_ax, calc_type p_par_ax_min, calc_type dp_perp_ax, calc_type p_perp_ax_min){

/** Lookup table returning values assuming LINEAR axes for p_par and perp and doing a weighed linear interpolation. Data is assumed to be 1-d with Fortran ordering, mapping to 2-d. Thus f(i, j) = data[j*par_sz+i] */
/** NOTE p_par is following Xiao and is NOT actually momentum*/
/** NB NB we assume 0 density outside lookup as we have no indication how to continue*/

  calc_type p_par_ind_decimal, p_perp_ind_decimal;
  long long p_par_ind, p_perp_ind;
  calc_type f_tmp_plus, f_tmp_minus, par_diff, perp_diff, interp_val;
  //Calculate exact axis position
  if(p_par > p_par_ax_min){
    p_par_ind_decimal = (p_par - p_par_ax_min)/dp_par_ax;
    p_par_ind = std::floor(p_par_ind_decimal);
  }else{
    p_par_ind = -1;
  }
  if(p_perp > p_perp_ax_min){
    p_perp_ind_decimal = (p_perp - p_perp_ax_min)/dp_perp_ax;
    p_perp_ind = std::floor( p_perp_ind_decimal);
  }else{
    p_perp_ind = -1;
  }

  if(p_par_ind >= par_sz || p_perp_ind >= perp_sz || p_par_ind <0 || p_perp_ind < 0){
    //If out of range, assume 0. Would also be reasonable to assume last in-range value.
    return 0.0;
  }
  par_diff = p_par_ind_decimal - p_par_ind;
  perp_diff = p_perp_ind_decimal - p_perp_ind;
  //4 point linear interpolation to "exact" location
  //Calculate f at lower and upper perp values interpolated on par lines
  f_tmp_minus = par_diff *data[par_sz*p_perp_ind + p_par_ind + 1] + (1.0 - par_diff)* data[par_sz*p_perp_ind + p_par_ind];
  f_tmp_plus = par_diff *data[par_sz*(p_perp_ind + 1) + p_par_ind + 1] + (1.0 - par_diff)* data[par_sz*(p_perp_ind + 1) + p_par_ind];
  
  interp_val = (perp_diff * f_tmp_plus + (1.0 - perp_diff)*f_tmp_minus);

  return interp_val;

}

calc_type seperable_lookup(calc_type p_par, calc_type p_perp, my_type * data, size_t par_sz, size_t perp_sz,calc_type dp_par_ax, calc_type p_par_ax_min, calc_type dp_perp_ax, calc_type p_perp_ax_min){

/** Lookup table for seperable functions (h(p_par)*g(p_perp)) returning values assuming LINEAR axes for p_par and perp and doing a weighed linear interpolation. Data is assumed to be sucessive 1-d arrays, f then g. Thus f(i, j) = data[i]*data[par_sz+j]*/

  calc_type p_par_ind_decimal, p_perp_ind_decimal;
  size_t p_par_ind, p_perp_ind;
  calc_type par_diff, perp_diff, interp_val;
  
  //Calculate exact axis position
  p_par_ind_decimal = (p_par - p_par_ax_min)/dp_par_ax;
  p_perp_ind_decimal = (p_perp - p_perp_ax_min)/dp_perp_ax;
  
  p_par_ind = std::floor(p_par_ind_decimal);
  p_perp_ind = std::floor( p_perp_ind_decimal);
  
  par_diff = p_par_ind_decimal - p_par_ind;
  perp_diff = p_perp_ind_decimal - p_perp_ind;
  
  //2 point linear interpolation to "exact" location on par and perp
  
  interp_val = (par_diff *data[p_par_ind + 1] + (1.0 - par_diff)*data[p_par_ind]) * (perp_diff *data[par_sz + p_perp_ind + 1] + (1.0 - perp_diff)*data[par_sz+p_perp_ind]);

  return interp_val;

}


std::function<calc_type(calc_type p_par, calc_type p_perp)> configure_lookup(std::string file_prefix, std::string file, non_thermal * my_nonth){
/** Read lookup data from file(s) (supply comma-separated list) and bind lookup function
*
*We're going to be simple and require a pre-space averaged distribution. Either one file, in which case we assume f(x, y), or two, when we assume f(x)*g(y)

*/

  my_print("Configuring lookup from "+file_prefix+file);
  std::function<calc_type(calc_type p_par, calc_type p_perp)> bound_lookup;

  calc_type p_par_ax_min, p_perp_ax_min;
  size_t par_sz, perp_sz;

  //Open and read file sizes

  bool sep = (file.find_first_of(',')!=std::string::npos);
  size_t tot_els;
  std::ifstream infile;
  
  if(sep){
    my_print("Detected seperable lookup");
    std::string file_par, file_perp;
    file_par = file.substr(0, file.find_first_of(','));
    file_perp = file.substr(file.find_first_of(',')+1, file.size());
    trim_string(file_par);
    trim_string(file_perp);

    data_array tmp_array = data_array(file_prefix+file_par);
    data_array tmp_array_perp = data_array(file_prefix+file_perp);

    {
      std::vector<size_t> par_dims, perp_dims;

      for(size_t i=0; i <tmp_array.get_dims(); i++) par_dims.push_back(tmp_array.get_dims(i));
      for(size_t i=0; i <tmp_array_perp.get_dims(); i++) perp_dims.push_back(tmp_array_perp.get_dims(i));
      //Throw error if dims are not 1 or all larger dims = 1
      size_t total =1;
      for(size_t i=1; i< par_dims.size(); i++) total *=par_dims[i];
      if(par_dims.size() != 1 || total > 1){
        my_print("Error, seperable array read is not 1-D");
        exit(1);
      }
      total=1;
      for(size_t i=1; i< perp_dims.size(); i++) total *=perp_dims[i];
      if(perp_dims.size() != 1 || total > 1){
        my_print("Error, seperable array read is not 1-D");
        exit(1);
      }
      par_sz = par_dims[0];
      perp_sz = perp_dims[0];
    }
    my_type * par_axes = tmp_array.disown_axes();
    my_type * par_data = tmp_array.disown_data();
    my_type * perp_axes = tmp_array_perp.disown_axes();
    my_type * perp_data = tmp_array_perp.disown_data();

    //Quickly check axes are linear binned, else throw an error. Retain min and inc if good
    p_par_ax_min = par_axes[0];
    p_perp_ax_min = perp_axes[0];
    my_type dp_par_ax = par_axes[1]-par_axes[0], dp_perp_ax = perp_axes[1]-perp_axes[0];
    bool err=0;
    for(size_t j=0; j<par_sz-1; j++){
        if( (par_axes[j+1]-par_axes[j] - dp_par_ax) > GEN_PRECISION*v0) err=1;
    }
    for(size_t j=0; j<perp_sz-1; j++){
        if( (perp_axes[j+1]-perp_axes[j] - dp_perp_ax) > GEN_PRECISION*v0) err=1;
    }
    if(err){
      my_print("Error, axes are not linear");
      exit(1);
    }
    free(par_axes);
    free(perp_axes);
    //Now we have data in two places. Realloc the parallel pointer to hold both and copy into place
    tot_els = par_sz + perp_sz;
    my_type * tmp = (my_type *) realloc((void*) par_data, tot_els*sizeof(my_type));
    if(!tmp){
      my_print("Error allocating memory");
      exit(1);
    }
    std::copy(perp_data, perp_data+perp_sz, tmp+ par_sz);
    //Copy perp data over
    my_nonth->lookup_data = tmp;
    //Keep pointer
    free(perp_data);
    //Free perp data. Par_data has been repurposed
    
    bound_lookup = std::bind(seperable_lookup, std::placeholders::_1, std::placeholders::_2, my_nonth->lookup_data, par_sz, perp_sz, dp_par_ax, p_par_ax_min, dp_perp_ax, p_perp_ax_min);

    //Bind the lookup function

  }else{
    trim_string(file);

    data_array tmp_array = data_array(file_prefix+file);
//    my_print("Array n_dims "+ mk_str(tmp_array.get_dims()));
    {
      std::vector<size_t> dims;

      for(size_t i=0; i <tmp_array.get_dims(); i++) dims.push_back(tmp_array.get_dims(i));
      //Throw error if dims are not 2 or all larger dims = 1
      size_t total =1;
      for(size_t i=2; i< dims.size(); i++) total *=dims[i];
      if(dims.size() != 2 || total > 1){
        my_print("Error, array read is not 2-D");
        exit(1);
      }
      par_sz = dims[0];
      perp_sz = dims[1];
    }
    my_type * axes = tmp_array.disown_axes();
    my_type * data = tmp_array.disown_data();

    //Quickly check axes are linear binned, else throw an error. Retain min and inc if good
    p_par_ax_min = axes[0];
    p_perp_ax_min = axes[par_sz+0];
    my_type dp_par_ax = axes[1]-axes[0], dp_perp_ax = axes[par_sz+1]-axes[par_sz+0];
    bool err=0;
    for(size_t j=0; j<par_sz-1; j++){
        if( (axes[j+1]-axes[j] - dp_par_ax) > GEN_PRECISION*v0) err=1;
    }
    for(size_t j=0; j<perp_sz-1; j++){
        if( (axes[par_sz+j+1]-axes[par_sz+j] - dp_perp_ax) > GEN_PRECISION*v0) err=1;
    }
    if(err){
      my_print("Error, axes are not linear");
      exit(1);
    }
    free(axes);
    my_nonth->lookup_data = data;
    //Keep pointer
    bound_lookup = std::bind(lookup, std::placeholders::_1, std::placeholders::_2, my_nonth->lookup_data, par_sz, perp_sz, dp_par_ax, p_par_ax_min, dp_perp_ax, p_perp_ax_min);

    //Bind the lookup function

  }
  my_nonth->dims[0] = par_sz;
  my_nonth->dims[1] = perp_sz;
  return bound_lookup;
}

