

//  Created by Heather Ratcliffe on 11/02/2016.
/** Very small class to hold a non-thermal electron distribution we can operate on. Mainly contains routine to parse from deck.status assuming particular input deck format.*/

#include <stdio.h>
#include <cmath>
#include <boost/math/special_functions.hpp>
//Provides Bessel functions, erf, and many more
#include <functional>
#include <boost/math/special_functions/sign.hpp>
#include "support.h"
#include "non_thermal.h"

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
  //Defaults which are deliberately meaningless....
  
  bool err = configure_from_file(file_prefix);
  if(err) my_print("WARNING: Error getting electron data from .status file", mpi_info.rank);

}

bool non_thermal::configure_from_file(std::string file_prefix){
/** \brief Setup non-thermal distribution from file
*
*Reads a [prefix]deck.status file and parses out the neccessary parameters. Note that this depends on the specific deck format used to get the parameter names \todo Multithermal
*/


  std::ifstream infile;
  infile.open(file_prefix+"deck.status");
  std::string header_row, line;
  infile>> header_row;
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
  etc where the RHS's are the string names in the deck and the function is either maxwell, kappa etc. ALL LOWER CASE. this then defines a plasma from deck constants*/
  /** \todo And in plasma, fix lower casing*/
  
  std::cout<<file_prefix+"nonthermal.conf"<<'\n';
  infile.open(file_prefix+"nonthermal.conf");
  if(!infile){
    my_print("Using generic bimaxwellian electrons");

    auto tmp_fn = std::bind(bimax, std::placeholders::_1, std::placeholders::_2, parameters[VPAR], parameters[VPERP], parameters[DENS_RAT]);
    f_p_private.push_back(tmp_fn);

  }else{
    
    std::string line, name, val;
    int block_num = -1;
    bool parse_err;
    size_t pos;

    //very naive parsing. We spin through until we find a ":"
    //if we don't find n_comps such blocks, we report and continue
    std::string function="", dens="", vpar="", vperp="", kappa="";
    std::function<calc_type(calc_type, calc_type)> tmp_fn;

    while(getline(infile, line)){
      if(line.find(':') != std::string::npos){
        //Create a plasma component function bound to these parameters and add to vector
        
        if(function == "max"){
          //Check we have all params
          for(auto nam : {dens, vpar, vperp}) if(!parameters.count(nam)) my_print("Param "+nam+" not found!!!");
          
          tmp_fn = std::bind(bimax, std::placeholders::_1, std::placeholders::_2, parameters[vpar], parameters[vperp], parameters[dens]);
          f_p_private.push_back(tmp_fn);
        }
        else if(function =="kappa"){
          for(auto nam : {dens, vpar, vperp, kappa}) if(!parameters.count(nam)) my_print("Param "+nam+" not found!!!");

          tmp_fn = std::bind(bikappa, std::placeholders::_1, std::placeholders::_2, parameters[kappa], parameters[vpar], parameters[vperp], parameters[dens]);
          f_p_private.push_back(tmp_fn);
        
        }
        block_num ++;
        function=""; dens=""; vpar=""; vperp="";kappa="";
        if(block_num >= ncomps) break;
        //Reset lookups
        continue;
        // is next block, skip this header line
      }
      if(block_num == -1 && line.find('=')){
      //This might be an n_comps definition!
        parse_err = parse_name_val(line, name, val);
        if(!parse_err && name =="ncomps") this->ncomps = atoi(val.c_str());
      }
      if(block_num >= 0){
        //this line is a valid input one, probably!
        parse_err = parse_name_val(line, name, val);
        if(!parse_err){
          //Grab the strings for each parameter
          if(name == "function") function = val;
          if(name =="dens") dens = val;
          if(name =="vpar") vpar = val;
          if(name =="vperp") vperp = val;
          if(name =="kappa") kappa = val;
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

calc_type non_thermal::f_p(calc_type p_par, calc_type p_perp){
/**\brief Return f(p)
*
*Evaluates function at p_perp, p_par and returns result. Calls bound function so we can use multiple functional forms for f, assuming always seperable into f_par, f_perp
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


