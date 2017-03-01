

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

/********Basic setup and allocation functions ****/
bool non_thermal::configure_from_file(std::string file_prefix){
/** \brief Setup non-thermal distribution from file
*
*Reads a [prefix]deck.status file and parses all parameters from constants block. Then, reads [prefix]nonthermal.conf and sets up a nonthermal electron distribution accordingly, using the named parameters
*/

  //Parse all constants from deck.status file into params map
  //Open file and locate constants block
  std::ifstream infile;
  infile.open(file_prefix+"deck.status");
  std::string line;
  bool found = false;
  if(!infile.is_open()){
    my_error_print("No deck.status file found", mpi_info.rank);
    return 1;
  
  }
  while(infile){
    getline(infile, line);
    if(!(strcmp(line.substr(0, CONSTANTS.size()).c_str(), CONSTANTS.c_str()))){
      found = true;
      break;
    }
  }
  if(!found){
    my_error_print("No constants dump found in deck.status", mpi_info.rank);
    return 1;
  }

  //Read all constants lines and parse params out
  std::string name, val;
  bool parse_err;
  float val_f;
  while(infile){
    getline(infile, line);
    if(!(strcmp(line.substr(0, CONSTANTS_END.size()).c_str(), CONSTANTS_END.c_str()))) break;
    parse_err = parse_name_val(line, name, val);
    if(parse_err) continue;
    val_f = atof(val.c_str());
    parameters[name] = val_f;
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
  ...
  end:
  etc where the RHS's are the string names in the deck and the function is "max", "kappa", "lookup" etc. ALL LOWER CASE. this then defines a plasma from deck constants*/
  /** \todo And in plasma, fix lower casing*/
  
  my_print("Reading "+file_prefix+"nonthermal.conf");
  infile.open(file_prefix+"nonthermal.conf");
  if(!infile){
    //Default to using a Bimaxwellian.
    my_print("Using generic bimaxwellian electrons");

    auto tmp_fn = std::bind(bimax, std::placeholders::_1, std::placeholders::_2, parameters[VPAR], parameters[VPERP], parameters[DENS_RAT]);
    f_p_private.push_back(tmp_fn);
    total_dens += parameters[DENS_RAT];

  }else{
    //Read conf file and create f_p component for each block according to specs
    size_t ncomps = 0;
    std::string line, name, val;
    size_t block_num = 0; //0th block is header stuff. First species is block 1 etc
    bool parse_err;

    //very naive parsing. We spin through until we find a ":"
    //if we don't find n_comps such blocks, we report and continue
    std::string function="", dens="", vpar="", vperp="", kappa="", lookup="", block_name="";
    std::function<calc_type(calc_type, calc_type)> tmp_fn;
    bool function_set = false;

    while(getline(infile, line)){
      function_set = false;
      if(line.find(':') != std::string::npos){
        //Each time we find a new block start line we try and create the component specified by the previous block
        calc_type norm;
        if(function !=""){
          //Some general setups for ANY function
          if(parameters.count(dens)!=0) total_dens += parameters[dens];
        }
  //------Binding function free parameters to create f_p component ------
        function = str_to_lower(function);
        if(function == "max"){
          //Check we have all params
          my_print("Adding component '"+block_name+"'");
          for(auto nam : {dens, vpar, vperp}) if(!parameters.count(nam)) my_error_print("Param "+nam+" not found!!!");

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
          //Check we have all params
          my_print("Adding component '"+block_name+"'");
          for(auto nam : {dens, vpar, vperp, kappa}) if(!parameters.count(nam)) my_error_print("Param "+nam+" not found!!!");

          tmp_fn = std::bind(bikappa, std::placeholders::_1, std::placeholders::_2, parameters[kappa], parameters[vpar], parameters[vperp], parameters[dens]);
          f_p_private.push_back(tmp_fn);
          function_set = true;

        }else if(function=="lookup"){
          //Check we have all params
          my_print("Adding component '"+block_name+"'");
          tmp_fn = configure_lookup(file_prefix, lookup, this);
          //Lookup config can fail and return empty function: we don't want to keep this
          if(tmp_fn){
            f_p_private.push_back(tmp_fn);
            function_set = true;
          }
        }
        //Report and reset for next block
        if(function_set) my_print("Configured function "+function);
        block_num ++;
        function=""; dens=""; vpar=""; vperp="";kappa="", lookup="";
        //Stop after ncomps blocks
        if(block_num > ncomps) break;
        //Set name of current block
        block_name = line.substr(0, line.find(':'));
        //Skip this header line and jump to read next line
        continue;
      }
      if(block_num == 0 && line.find('=')){
      //This might be an n_comps definition or nonrely flag!
        parse_err = parse_name_val(line, name, val);
        if(!parse_err){
          name = str_to_lower(name);
          if(name =="ncomps") ncomps = atoi(val.c_str());
          else if(name =="nonrely"){
            if(val.size() >0 && val[0] == '1'){
              this->norely=1;
              my_print("Using nonrelativistic calculation!");
            }
          }else{
            my_error_print("!!!!!Unknown name "+name);
          }
          
        }
      }
      if(block_num > 0){
        //This line is a valid input one, probably!
        parse_err = parse_name_val(line, name, val);
        if(!parse_err){
          //Our strings are matched without case-sensitivity, but deck keys retain casing
          name = str_to_lower(name);
          //Grab the strings for each parameter
          if(name == "function") function = val;
          else if(name =="dens") dens = val;
          else if(name =="vpar") vpar = val;
          else if(name =="vperp") vperp = val;
          else if(name =="kappa") kappa = val;
          else if(name =="lookup") lookup = val;
          else my_error_print("!!!!Unknown name "+name);
        }
      }
      
    }
    if(block_num > ncomps+1){
    my_print("Too many blocks in plasma file, ignoring", mpi_info.rank);
    }
    else if(block_num < ncomps){
      my_print(mk_str(block_num), mpi_info.rank);
  
      my_print("Insufficient blocks in config file", mpi_info.rank);
    }
  }
  return 0;

}

non_thermal::non_thermal(std::string file_prefix){
/** \brief Construct non-thermal distrib
*
*Sets default params
*/
  //Base density is 1
  total_dens = 1.0;
  norely = 0;
  
  dp = 10.0; //Sensible dp for electrons (m/s)
  
  bool err = configure_from_file(file_prefix);
  if(err) my_error_print("WARNING: Error getting electron data from .status file", mpi_info.rank);

}

non_thermal::~non_thermal(){
/** \brief Clean up. 
*
*
Calls clean_lookup() which can be used to do anything needed to cleanup after a lookup function
*/

  clean_lookup(*this);

}

/********Primary interface functions ****/
calc_type non_thermal::d_f_p(calc_type p_par, calc_type p_perp, bool parallel){
/** Simple Numerical derivative*/
  
  if(parallel) return (f_p(p_par+dp, p_perp) -f_p(p_par, p_perp))/dp;
  else return (f_p(p_par, p_perp+dp) -f_p(p_par, p_perp))/dp;
  
}

calc_type non_thermal::f_p(calc_type p_par, calc_type p_perp){
/**\brief Return value of f(p_par, p_perp)
*
*Evaluates current specification of f at p_perp, p_par and returns result. If f has multiple components these are summed.
*/

  calc_type ret = 0;
  for(size_t i=0; i<f_p_private.size(); i++) ret +=f_p_private[i](p_par, p_perp);
  return ret;

}

/********Functions providing component specifications ****/
//These must all take at least 2 parameters. p_par is always the first parameter, p_perp the second

calc_type bimax(calc_type p, calc_type p2, calc_type p_th, calc_type p_th2, calc_type A){
/** Bimaxwellian (different parallel and perp temperatures)*/
  return A*std::exp(-p*p/p_th/p_th - p2*p2/p_th2/p_th2);
}

calc_type bikappa(calc_type p, calc_type p2, calc_type kappa, calc_type v_kpar, calc_type v_kperp, calc_type A){
/* Bikappa function*/
  return A*std::pow( 1.0 + (std::pow(p/v_kpar, 2) + std::pow(p2/v_kperp, 2))/kappa , -kappa-1);
}

calc_type lookup(calc_type p_par, calc_type p_perp, my_type * data, size_t par_sz, size_t perp_sz,calc_type dp_par_ax, calc_type p_par_ax_min, calc_type dp_perp_ax, calc_type p_perp_ax_min){

/** Lookup table returning values assuming LINEAR axes for p_par and perp and doing a weighed linear interpolation. Data is assumed to be 1-d with Fortran ordering, mapping to 2-d. Thus f(i, j) = data[j*par_sz+i] The axes run from p_par_ax_min and p_perp_ax_min up in steps of dp_par_ax and dp_perp_ax respectively. */
/** NOTE p_par is following Xiao and is NOT actually momentum*/
/** NB NB we assume 0 density outside lookup as we have no indication how to continue*/

  calc_type p_par_ind_decimal, p_perp_ind_decimal;
  long long p_par_ind, p_perp_ind;
  calc_type f_tmp_plus, f_tmp_minus, par_diff, perp_diff, interp_val;
  //Calculate axis index for par and perp
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
  //Check we're in range of the table
  //Casts here should only ever widen the type
  if(p_par_ind >= (long long) par_sz || p_perp_ind >= (long long)perp_sz || p_par_ind <0 || p_perp_ind < 0){
    //If out of range, assume 0. Would also be reasonable to assume last in-range value.
    return 0.0;
  }
  //Calculate the fraction of cell for interpolation
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

/** Lookup table for seperable functions (h(p_par)*g(p_perp)) returning values assuming LINEAR axes for p_par and perp and doing a weighed linear interpolation. Data is assumed to be sucessive 1-d arrays, f then g. Thus f(i, j) = data[i]*data[par_sz+j] The axes run from p_par_ax_min and p_perp_ax_min up in steps of dp_par_ax and dp_perp_ax respectively. */

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
  
  //2 point linear interpolation to "exact" location on par and perp, and product the two for final value
  if(p_par_ind < par_sz && p_perp_ind < perp_sz){
    interp_val = (par_diff *data[p_par_ind + 1] + (1.0 - par_diff)*data[p_par_ind]) * (perp_diff *data[par_sz + p_perp_ind + 1] + (1.0 - perp_diff)*data[par_sz+p_perp_ind]);
  }else{
    interp_val = 0.0;
  }
  return interp_val;

}

/********Helpers for component functions ****/
std::function<calc_type(calc_type p_par, calc_type p_perp)> configure_lookup(std::string file_prefix, std::string file, non_thermal * my_nonth){
/** Read lookup data from file(s) (supply comma-separated list) and bind lookup function
*
*We're going to be simple and require a pre-space averaged distribution. Either one file, in which case we assume f(x, y), or two, when we assume f(x)*g(y). We check for linear axes and store the dp values.

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
        my_error_print("Error, seperable array read is not 1-D");
        return bound_lookup;
      }
      total=1;
      for(size_t i=1; i< perp_dims.size(); i++) total *= perp_dims[i];
      if(perp_dims.size() != 1 || total > 1){
        my_error_print("Error, seperable array read is not 1-D");
        return bound_lookup;
      }
      par_sz = par_dims[0];
      perp_sz = perp_dims[0];
    }
    //Steal the data and axes memory from arrays
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
      my_error_print("Error, axes are not linear");
      return bound_lookup;
    }
    free(par_axes);
    free(perp_axes);
    //Now we have data in two places. Realloc the parallel pointer to hold both and copy into place
    tot_els = par_sz + perp_sz;
    my_type * tmp = (my_type *) realloc((void*) par_data, tot_els*sizeof(my_type));
    if(!tmp){
      my_error_print("Error allocating memory");
      return bound_lookup;
    }
    //Copy perp data over
    std::copy(perp_data, perp_data+perp_sz, tmp+ par_sz);
    //Keep pointer
    my_nonth->lookup_data = tmp;
    //Free perp data. Par_data has been repurposed
    free(perp_data);
    
    //Bind the lookup function
    bound_lookup = std::bind(seperable_lookup, std::placeholders::_1, std::placeholders::_2, my_nonth->lookup_data, par_sz, perp_sz, dp_par_ax, p_par_ax_min, dp_perp_ax, p_perp_ax_min);

  }else{
    //Non-separable lookup, i.e. 2-D array
    trim_string(file);

    data_array tmp_array = data_array(file_prefix+file);
    {
      std::vector<size_t> dims;

      for(size_t i=0; i <tmp_array.get_dims(); i++) dims.push_back(tmp_array.get_dims(i));
      //Throw error if dims are not 2 or all larger dims = 1
      size_t total =1;
      for(size_t i=2; i< dims.size(); i++) total *=dims[i];
      if(dims.size() != 2 || total > 1){
        my_error_print("Error, array read is not 2-D");
        return bound_lookup;
      }
      par_sz = dims[0];
      perp_sz = dims[1];
    }
    //Steal data and array memory from array
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
      my_error_print("Error, axes are not linear");
      return bound_lookup;
    }
    free(axes);
    //Keep pointer
    my_nonth->lookup_data = data;

    //Bind the lookup function
    bound_lookup = std::bind(lookup, std::placeholders::_1, std::placeholders::_2, my_nonth->lookup_data, par_sz, perp_sz, dp_par_ax, p_par_ax_min, dp_perp_ax, p_perp_ax_min);
  }
  return bound_lookup;
}

void clean_lookup(non_thermal & my_nonth){
/** \brief Cleanup lookup configuration
*
*Any cleanup from configure_lookup should go here. Currently we free the lookup data pointer in the non_thermal class itself. If configure lookup does anything which allocates memory etc, AND the only references to that memory are in the bound function, then we need to keep some other reference and free it here, or we'll leak.
*/
  if(my_nonth.lookup_data) free(my_nonth.lookup_data);
  
}

