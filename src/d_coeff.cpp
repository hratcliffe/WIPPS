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
#include "controller.h"
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

  my_controller = nullptr;

  n_thetas = 10;
  n_n = 4;
  // FAKENUMBERS
}

diffusion_coeff::~diffusion_coeff(){
/** Delete interconnections */


}

void diffusion_coeff::set_ids(float time1, float time2, int space1, int space2, int wave_id, char block_id[10]){

  this->time[0] = time1;
  this->time[1] = time2;
  this->space[0] = space1;
  this->space[1] = space2;
  strcpy(this->block_id, block_id);
  this->wave_id = wave_id;
}


bool diffusion_coeff::write_to_file(std::fstream &file){

  file<<"";
  return 0;

}

void diffusion_coeff::calculate(){
// calls whatever auxilliarlies it needs to calc D. This is local as fn of n, x

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

  plasma * plas;
  spectrum * spect;
  if(my_controller){
    plas = my_controller->my_plas;
    spect = my_controller->my_spect;
  }
  else{
    std::cout<<"No controller"<<std::endl;
    return;
  }
  
  this->copy_ids(spect);
  //copy block id, ranges etc from spect.

  calc_type theta, omega,lat, alpha, omega_n=0.0, x, x_inc, D_tmp;
  calc_type Eq6, mu_dom_mu, Eq7, dmudx, v_par, alpha_inc;

  std::vector<calc_type> omega_calc;

  alpha = 0.0;
  alpha_inc =  2.0*pi/dims[1]/4.0;
  //do one quadrant...
  x_inc = 4.0/n_thetas; //To cover range from 0 to 2...
  x= - n_thetas/ 2.0 * x_inc;

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

  for(int i=0; i<n_thetas; ++i) dx[i] = x_inc;
  //linear grid for now
  
  
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
        theta = atan(x);
        //theta = pi;
        //x = tan(theta);
        std::cout<<j<<" "<<theta/pi<<"+++++++++++++++++++++++++"<<std::endl;

        D_tmp = 0.0;
        for(int n=-n_n; n<n_n; ++n){
//        { int n=-1;//std::cout<<n<<"-----------------"<<std::endl;
          // n is resonant number
          omega_calc = plas->get_omega(x, v_par, (calc_type) n);
          //if(std::abs(omega/my_const.omega_ce) > 1.0) continue;
          if(omega_calc.size()==0) continue; //redundant?
          for(size_t ii =0; ii< omega_calc.size(); ++ii){
            omega = omega_calc[ii];
            std::cout<<"Freq is "<<omega/my_const.omega_ce<<std::endl;
            my_mu = plas->get_phi_mu_om(omega, theta, alpha, n, omega_n);

            if(!my_mu.err) std::cout<<"YAY!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! "<<n<<" "<<j<<std::endl;
            mu_dom_mu = my_mu.mu + omega * my_mu.dmudom;
            dmudx = 0.0;
            // FAKENUMBERS
            Eq6 = omega/(omega - omega_n)* my_mu.mu/mu_dom_mu;
            Eq7 = -1.0* (my_mu.mu*omega_n/(omega*(omega-omega_n)) - my_mu.dmudom)/(my_mu.mu *sin(theta)*cos(theta) - dmudx);
            //Need this iff we use second expression in Eq 5
            
            D_tmp += my_mu.phi; // FAKENUMBERS This will be the n summed D so add onto it each n iteration
          }
        }
        //Store into theta array. Might not need tmp in the end
        //Except that it saves us one multiplication per iteration without acting on element in place after.
        D_theta[j] = D_tmp*pow(cos(theta), 2);
        x+= dx[j];

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



