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

diffusion_coeff::diffusion_coeff(int nx, int n_angs):data_array(nx, n_angs){

}


void diffusion_coeff::set_ids(float time1, float time2, int space1, int space2, int wave_id, char block_id[10], int function_type){

}


bool diffusion_coeff::write_to_file(std::fstream &file){


  return 0;

}

void diffusion_coeff::calculate(spectrum * spect, plasma * plas){
//takes spectrum and plasma, and calls whatever auxilliarlies it needs to calc D. This is local as fn of n, x

//Which one? eneergy or angle. They're related trivially though. Pick one.
//When to sum over resonances?


/*Subsections:
To solve Eq 2 in Albert:
Get G1

Get phi (note also needs mu internally)
Get mu, dmu/domega which are used to:
  get G2 - needs mu, dmu/domega
  get denominator | 1 - (d omega / d kparallel)_x | from 6 - needs mu, dmu/domega
  get Eq 7 factor - needs mu, dmu/domega, dmu/dx

*/

//loop over th, w, psi, alpha, n

  //To DO!
  //get n ranges... and add in
  //resonant freq relation to x


  calc_type theta, omega, psi, alpha, omega_n, inc, x, x_inc, D_tmp;
  calc_type Eq6, mu_dom_mu, Eq7, dmudx;
  int n;
  alpha = 0.0;
  inc =  2.0*pi/dims[1]/4.0;
  //do one quadrant...
  
  x_inc = 0.1;
  // FAKENUMBERS
  //Allocate and construct a dx from theta.
  calc_type * dx, * D_theta;
  dx = (calc_type *) calloc(dims[1], sizeof(calc_type));
  D_theta = (calc_type *) calloc(dims[1], sizeof(calc_type));

  for(int i=0; i<dims[1]; i++) dx[i] = i*x_inc;
  
  //innermost loop should be n. Next theta, as we need mu at each theta. Omega and x are interchangeable...
  //Alpha remains, as does particle v.

  for(int i =0; i<1; i++){
    //This will be particle velocity
    for(int k =0; k<10; k++){
    //outer loop for longer for timing
    //This will be alpha loop
      alpha += inc;

      for(int j=0;j<this->dims[1]; j++){
      //theta loop
        x = x + dx[j];
        theta = atan(x);
        //calc_type th, calc_type w, calc_type psi)
        mu my_mu = plas->get_root(theta, omega, psi);
        mu_dom_mu = my_mu.mu + omega * my_mu.dmudom;
        dmudx = my_mu.dmudth; //transform theta to tan theta.
        // FAKENUMBERS
        
        for(n=0; n<1000; n++){

          //calc_type th, calc_type w, calc_type psi, calc_type alpha, int n
          calc_type phi = plas->get_phi(theta, omega, psi, alpha, n, omega_n);

          Eq6 = omega/(omega - omega_n)* my_mu.mu/mu_dom_mu;
          Eq7 = -1.0* (my_mu.mu*omega_n/(omega*(omega-omega_n)) - my_mu.dmudom)/(my_mu.mu *sin(psi)*cos(psi) - dmudx);
        
          D_tmp = 1.0; // FAKENUMBERS This will be the n summed D
        //now total up in n...
        
        }
        //Store into theta array
        D_theta[j] = D_tmp;
      }
      //now integrate in x = tan theta
      D_tmp = integrator(D_theta, dims[1], dx);
    }
  }
}


