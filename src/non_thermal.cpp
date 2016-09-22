

//  Created by Heather Ratcliffe on 11/02/2016.
/** Very small class to hold a non-thermal electron distribution we can operate on. Mainly contains routine to parse from deck.status assuming particular input deck format.*/

#include <stdio.h>
#include <cmath>
#include <boost/math/special_functions.hpp>
//Provides Bessel functions, erf, and many more
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
  //Defaults which are deliberately meaningless....
  
  bool err = configure_from_file(file_prefix);
  if(err) my_print("WARNING: Error getting electron data from .status file", mpi_info.rank);

}

non_thermal::~non_thermal(){



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
  for(size_t i=0; i< lines.size(); i++){

    parse_err = parse_name_val(lines[i], name, val);
    if(parse_err) continue;
    val_f = atof(val.c_str());
  
    if(name == OMEGA_CE) this->ref_B = val_f * me/std::abs(q0);
    else if(name == DENS) this->ref_dens = val_f;
    else if(name == VPAR) this->v_par = val_f;
    else if(name == VPERP) this->v_perp = val_f;
    else if(name == DENS_RAT) this->fraction = val_f;
    
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
