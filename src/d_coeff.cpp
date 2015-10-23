//
//  d_coeff.cpp
//  
//
/** \file d_coeff.cpp This file will calculate the diffusion coefficients from a spectrum. We'll try first just going for it. Does it make sense to have this be any sort of object? Possibly, in that we want to keep info on how, where, etc, etc, and the axes. So possibly it is again descended from a data_array. For now, implement some basic calcs
*
* @author Heather Ratcliffe @date 23/09/2015.
*/

#include <math.h>
#include <algorithm>
#include <boost/math/special_functions.hpp>
#include "support.h"
#include "my_array.h"
#include "spectrum.h"
#include "plasma.h"
#include "d_coeff.h"

extern mpi_info_struc mpi_info;

extern deck_constants my_const;

/**First we want to generate a particle velocity axis, and a pitch angle axis. Possibly the latter matches the one from spectrum?
*We need the refractive index. 
*Probably numerical integration/differentiation routines


*/

diffusion_coeff::diffusion_coeff(int nx, int n_angs):data_array(nx, n_angs){

  n_thetas = 10;
  n_n = 0;
  // FAKENUMBERS
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
//HERE we want to end up with D_alpha_alpha

  //To DO!
  //get n ranges... and add in
  //resonant freq relation to x
  //We probably want to set the n and theta ranges in another routine... or using UI data

  //Setup v_par, alpha ranges in axes

  this->copy_ids(spect);
  //copy block id, ranges etc from spect.

  calc_type theta, omega,lat, alpha, omega_n, inc, x, x_inc, D_tmp;
  calc_type Eq6, mu_dom_mu, Eq7, dmudx, v_par, alpha_inc;

  alpha = 0.0;
  alpha_inc =  2.0*pi/dims[1]/4.0;
  //do one quadrant...
  
  x_inc = 0.02/n_thetas; //To cover range from 0 to 2...
  lat = 0.0;
  // FAKENUMBERS
  //Allocate and construct a dx from theta.
  calc_type * dx, * D_theta, * v_axis;
  calc_type om_res = 10;
  // FAKENUMBERS for om_res..
  calc_type D_consts = 0.5* pi*my_const.omega_ce *pow(v0, 3)/om_res;

  dx = (calc_type *) calloc(dims[1], sizeof(calc_type));
  D_theta = (calc_type *) calloc(dims[1], sizeof(calc_type));
  //Want theta res to match alpha??
  v_axis = (calc_type *) calloc(dims[0], sizeof(calc_type));
  
  mu_dmudom my_mu;

  for(int i=0; i<n_thetas; ++i) dx[i] = i*x_inc;
  
  //innermost loop should be n. Next theta, as we need mu at each theta. Omega and x are interchangeable...
  //Alpha remains, as does particle v.

//-------------------Main loops here----------------------------

  for(int i =0; i< ((1< dims[0]) ? 1:dims[0]); ++i){
    //particle parallel velocity
    v_par = v_axis[i];
    v_par = 0.15 * v0;
    // FAKENUMBERS
    for(int k =0; k< ((1 <dims[1]) ? 1: dims[1]); k++){
      //particle pitch angle
      alpha += alpha_inc;

      for(int j=0;j<n_thetas; ++j){
      //theta loop for wave angle or x=tan theta
        x = x + dx[j];
        theta = atan(x) + pi;
//        omega = plas->get_omega(x); // FIX!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        //calc_type th, calc_type w, calc_type psi)
//        mu my_mu = plas->get_root(lat, omega, theta);
        D_tmp = 0.0;
        //for(int n=-n_n; n<n_n; ++n)
        {int n = 1;
          // n is resonant number
          //calc_type th, calc_type w, calc_type psi, calc_type alpha, int n
          omega = plas->get_omega(x, v_par, (calc_type) n);
          std::cout<<omega<<std::endl;
//          calc_type x, calc_type v_par, calc_type wc, calc_type n
          my_mu = plas->get_phi_mu_om(lat, omega, theta, alpha, n, omega_n);
          std::cout<< "mu "<<my_mu.mu<<" "<< my_mu.dmudom<<std::endl;
         // std::cout<< n<<" "<<my_mu.err<<std::endl;
//          if(my_mu.err) std::cout<<"Halp! Dispersion went wrong!..."<<std::endl;
          mu_dom_mu = my_mu.mu + omega * my_mu.dmudom;
//          dmudx = my_mu.dmudtheta * pow(cos(theta), 2); ; //transform theta to tan theta.
        // FAKENUMBERS

          
//          calc_type phi = plas->get_phi(lat, omega, theta, alpha, n, omega_n);

          Eq6 = omega/(omega - omega_n)* my_mu.mu/mu_dom_mu;
          Eq7 = -1.0* (my_mu.mu*omega_n/(omega*(omega-omega_n)) - my_mu.dmudom)/(my_mu.mu *sin(theta)*cos(theta) - dmudx);
          //Need this iff we use second expression in Eq 5
          
          D_tmp += my_mu.phi; // FAKENUMBERS This will be the n summed D so add onto it each n iteration
        
        }
        //Store into theta array. Might not need tmp in the end
        //Except that it saves us one multiplication per iteration without acting on element in place after.
        D_theta[j] = D_tmp*pow(cos(theta), 2);
      }
      //now integrate in x = tan theta
      D_tmp = integrator(D_theta, dims[1], dx);
      //Leaves us with one D for each alpha and v.
      //Store into data at i, k
      set_element(i, k, D_tmp*D_consts);
    
    }
  }
//------------End main loops-----------------------------------

}



