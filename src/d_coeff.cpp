//
//  d_coeff.cpp
//  
//
/** \file d_coeff.cpp This file will calculate the diffusion coefficients from a spectrum. We'll try first just going for it. Does it make sense to have this be any sort of object? Possibly, in that we want to keep info on how, where, etc, etc, and the axes. So possibly it is again descended from a data_array. For now, implement some basic calcs
*
* @author Heather Ratcliffe @date 23/09/2015.
*/

#include <math.h>
#include <boost/math/special_functions.hpp>
#include "support.h"
#include "my_array.h"
#include "spectrum.h"
#include "plasma.h"
#include "d_coeff.h"

extern mpi_info_struc mpi_info;

/**First we want to generate a particle velocity axis, and a pitch angle axis. Possibly the latter matches the one from spectrum?
*We need the refractive index. 
*Probably numerical integration/differentiation routines


*/

diffusion_coeff::diffusion_coeff(int nx, int n_ang):data_array(nx, n_ang){


}


void diffusion_coeff::set_ids(float time1, float time2, int space1, int space2, int wave_id, char block_id[10], int function_type){

}


bool diffusion_coeff::write_to_file(std::fstream &file){


  return 0;

}

void diffusion_coeff::calculate(spectrum * spect, plasma * my_mu){
//takes spectrum and plasma, and calls whatever auxilliarlies it needs to calc D. This is local as fn of n, x

//Which one? eneergy or angle. They're related trivially though. Pick one.
//When to sum over resonances?


/*Subsections:
To solve Eq 2 in Albert:
Get G1


Get phi, mu, dmu/domega which are used to:
  get G2 - needs mu, dmu/domega
  get denominator | 1 - (d omega / d kparallel)_x | from 6 - needs mu, dmu/domega
  get Eq 7 factor - needs mu, dmu/domega, dmu/dx

*/




}


