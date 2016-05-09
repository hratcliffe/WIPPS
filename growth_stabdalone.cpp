
#include <cmath>
#include <math.h>
#include <fstream>
#include <iostream>
//#include <limits>
#include <stdio.h>
#define calc_type double

const calc_type pi = 3.14159265359;/**< Pi */
const calc_type v0 = 2.997924e8; /**< Speed of light in m/s^2 */

template<typename T> T integrator(T * start, int len, T * increment);
calc_type * make_momentum_axis(int n_momenta, calc_type v_max_over_c);

int main(){
//Do a sample growth rate calculation to match Xiao et al paper
//This duplocates my main get_growth_rate but without using plasma etc etc class

  const int n_momenta = 10000;
  const int n_trials = 2000;

  calc_type * p_axis, *growth_rate;
  calc_type min_v = 0.0, max_v = 0.95;
  p_axis = make_momentum_axis(n_momenta, max_v);
  growth_rate = (calc_type*) malloc(n_trials* sizeof(calc_type));

/*// These reporoduce Xiao 2.a
  calc_type my_elec_fraction = 1.0000000000e-02;
  calc_type my_elec_v_par =  93776865.15999/sqrt(2.0);
  calc_type A = 2.0;
  //Anisotropy
*/

  calc_type my_elec_fraction = 1.0000000000e-02;
  calc_type my_elec_v_par =  44968868.700;
  calc_type A = 1.9;
//These are my parameters


  calc_type my_elec_v_perp = my_elec_v_par*sqrt(2.0*(A + 1.0));   
//  std::cout<<my_elec_v_par<<" "<<my_elec_v_perp/my_elec_v_par<<'\n';
  std::ofstream outfile;
  outfile.open("Growth.dat");
  outfile<<n_trials<<"\n";
  outfile<<"Omega \t Growth rate:"<<"\n";
  outfile<<"BEGIN"<<"\n";

  calc_type k;

  calc_type ck_om, om_ce = 10000.0, om_pe = 3.0*om_ce, om_diff, omega_in;

  calc_type f_tmp, a_par, a_perp, v_tmp, norm_f, S_tot=0.0, S_full_tot=0.0, dp, A_crit, A_rel, eta_rel;
  calc_type  *S, *S_full, * dp_ax;
  calc_type gamma, p_res, Delta_res;

  S = (calc_type *)malloc(n_momenta*sizeof(calc_type));
  S_full = (calc_type *)malloc(n_momenta*sizeof(calc_type));
  dp_ax = (calc_type *)malloc(n_momenta*sizeof(calc_type));
  
  calc_type d_om = std::abs(om_ce) / (float) (n_trials-1);
  omega_in = 0.0;
  int om_orders = 3;
  calc_type d_i = (calc_type) (n_trials -1)/(calc_type) om_orders;

  //Orders of magnitude to cover
  d_om = std::abs(om_ce)/ std::pow(10, om_orders);
  //"Logarithmic" axis of n_trials points up to 10^om_orders times om_ce

  for(int i=0; i< n_trials; i++){
    omega_in = std::pow(10, (calc_type) i /d_i)*d_om;

    //Calculate k using Eq 9 of Xiao
    k = std::sqrt( (std::pow(omega_in, 2) - std::pow(om_pe, 2)*omega_in/om_diff)/std::pow(v0, 2));
    ck_om = v0*k/omega_in;
    om_diff = omega_in - om_ce;
    
    //Get RMS momenta from velocities...
    // a_x = RMS p_x (Note factor of 2 in perp, not in par...)
    v_tmp = my_elec_v_par;
    a_par = std::sqrt(2.0)*v_tmp / std::sqrt(1.0 - (v_tmp/v0)*(v_tmp/v0));

    v_tmp = my_elec_v_perp;
    a_perp = v_tmp / std::sqrt(1.0 - (v_tmp/v0)*(v_tmp/v0));

    norm_f = 1.0/(a_perp*a_perp*a_par * pi * std::sqrt(pi));
//    std::cout<<a_par<<" "<<a_perp<<'\n';   
    for(int j=0; j< n_momenta; ++j){
    
      gamma = - 1.0 + ck_om * std::sqrt( (ck_om*ck_om -1.0 )*(1.0 + p_axis[j]*p_axis[j]/v0/v0)*(omega_in*omega_in/om_ce/om_ce) + 1.0 );
      gamma /= ((ck_om*ck_om - 1.0)*omega_in/om_ce);
      //14 in Xiao resonant gamma factor

      p_res = (gamma * omega_in - om_ce)/k;
      //Resonant momentum
      
      Delta_res = 1.0 - (omega_in*p_res / (v0*v0*k*gamma));
      //Xiao 15, no meaning given. Always +ve
      if(Delta_res < 1e-15) std::cout<<"ERROR!!"<<std::endl;
      
      //For f Maxwellian as Xiao 28: d f/ dp_x = 2 p_x a_x
      //Now p_par = p_res and p_perp is p_axis[j]
      f_tmp = norm_f * std::exp(- (p_res*p_res/(a_par*a_par)) - (p_axis[j]*p_axis[j]/(a_perp*a_perp)));

      S[j] = -2.0 * std::pow(p_axis[j], 3) * f_tmp / Delta_res;

      S_full[j] = 2.0 * std::pow(p_axis[j], 3) * f_tmp / Delta_res *(omega_in - om_ce/gamma) * (1.0 - a_perp*a_perp/a_par/a_par);
      //Both of these have removed factor of a_perp**2

    }

    dp_ax[0] = 0.0;
    for(int j=1; j<n_momenta; j++) dp_ax[j] = p_axis[j] - p_axis[j-1];

    S_tot = integrator(S, n_momenta, dp_ax);
    S_full_tot = integrator(S_full, n_momenta, dp_ax);

    calc_type ret = 0.0;
    
    if(std::abs(S_tot) > 0){
      A_crit = - omega_in /om_diff;
      A_rel = S_full_tot/S_tot/om_diff;
      
      eta_rel = pi * my_elec_fraction* om_diff/k * S_tot / a_perp/a_perp;
      
      ret = pi*om_pe*om_pe/(2.0*omega_in + om_pe*om_pe*om_ce/(std::pow(om_diff, 2))) * eta_rel * (A_rel - A_crit);
    
      //ret = A_rel;

    }
    growth_rate[i] = ret;
 //   outfile<<omega_in<<" "<<growth_rate[i]<<"\n";
    outfile<<omega_in<<" "<<ret<<"\n";
  }
  

  free(S);
  free(S_full);
  free(dp_ax);

}

calc_type * make_momentum_axis(int n_mom, calc_type v_max){

  //Max velocity to consider norm'd to c
  calc_type dp = v_max * v0/ std::sqrt(1.0- v_max*v_max) / (calc_type) (n_mom - 1);
  //Momentum step size, including gamma

  calc_type * p_axis;
  p_axis = (calc_type*) malloc(n_mom*sizeof(calc_type));
  if(!p_axis) return nullptr;

  p_axis[0] = 0.0;
  for(int i=1; i< n_mom; ++i) p_axis[i] = p_axis[i-1] + dp;
  //Set momenta from 0 to v_max. Do this once.
  if( std::abs(p_axis[n_mom-1] -  (v_max * v0/ std::sqrt(1.0- v_max*v_max)))> 1e-4) std::cout<<"Momentum axis errorrrrrrr";
  
  return p_axis;

}

template<typename T> T integrator(T * start, int len, T * increment){
/** \brief Basic numerical integrator
*
*Uses trapezium rule. WARNING this is working with contiguous memory. Not very C++ but faster.
*/

  T value=0.0;
  
  for(int i=0; i<len-1; i++){
  
    value += 0.5*(start[i] + start[i+1]) * increment[i];
    
  }
//  value += start[len-1]*increment[len-1];
  //top bnd we assume flat

 return value;

}


