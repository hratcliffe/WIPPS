/** \file main_growth.cpp \brief Helper program to get expected and/or actual growth rates of waves
*
* Can take sdf files or derived spectrum files, extract the approximate wave growth rates (peak or integrated) and output these, plus a theoretical rate, OR just the theoretical rate. Parameters are obtained as in Main, so requires a plasma.conf file and deck.status file.
* Depends on the SDF file libraries, the FFTW library. A set of test arguments is supplied. Call using ./growth `<growth_test_pars` to use these.
  \author Heather Ratcliffe \date 11/02/2016.
*/


#include <math.h>
#include <cmath>
#include <fstream>
#include <boost/math/special_functions.hpp>
//Provides Bessel functions, erf, and many more
#include <iostream>
#include <stdio.h>
#include "sdf.h"
//SDF file libraries
#include <mpi.h>
#include <complex.h>
#include <fftw3.h>
//FFTW3 Fourier transform libs

#include "support.h"
#include "main_support.h"
#include "reader.h"
#include "controller.h"
#include "plasma.h"
#include "my_array.h"
#include "d_coeff.h"
#include "spectrum.h"
#include "non_thermal.h"


deck_constants my_const;/**< Physical constants*/
extern const mpi_info_struc mpi_info;/**< Link to mpi_info as const*/

struct g_args{

  int src;/**<Whether to process no files (analytic only), sdf or spectrum 0, 1,2, respectively*/
  std::string spect_file;/**<File listing spectra if using this option*/
};


calc_type * make_momentum_axis(int n_momenta, calc_type v_max_over_c);

calc_type get_growth_rate(plasma * my_plas, non_thermal * my_elec, int n_momenta, calc_type * p_axis, calc_type omega_in);
g_args g_command_line(int argc, char * argv[]);

void write_growth_header(std::string in_file, plasma * my_plas, non_thermal * my_elec, int n_momenta, calc_type min_v, calc_type max_v, std::ofstream &outfile);

void write_growth(calc_type omega, calc_type growth, std::ofstream &outfile);


int main(int argc, char *argv[]){
/**
*In theory, these classes and functions should be named well enough that the function here is largely clear. Remains to be seen, eh?
*
*/

  int err;
  
  int ierr = local_MPI_setup(argc, argv);
  if(ierr){
    std::cout<< "Error initialising MPI. ABORTING!";
    return 1;
  }

  my_print(std::string("Code Version: ")+ VERSION, mpi_info.rank);
  my_print("Code is running on "+mk_str(mpi_info.n_procs)+" processing elements.", mpi_info.rank);

  MPI_Barrier(MPI_COMM_WORLD);

  setup_args cmd_line_args = process_command_line(argc, argv);
  g_args extra_args = g_command_line(argc, argv);

  if(mpi_info.rank == 0) get_deck_constants(cmd_line_args.file_prefix);
  share_consts();
  /** Get constants from deck and share to other procs*/




  if(extra_args.src == 1){
    //Read SDF files and derive...
    my_print("Processing sdf files", mpi_info.rank);
  
  
  
  }else if(extra_args.src == 2){
    //Read spectrum files and derive...
    my_print("Processing spectrum files", mpi_info.rank);
  
  
  
  
  }
  
  //Now we do the analytics.

  my_print("Calculating growth rates", mpi_info.rank);

  plasma * my_plas = new plasma(my_const.omega_ce * me/std::abs(q0), cmd_line_args.file_prefix);

  non_thermal * my_elec = new non_thermal(cmd_line_args.file_prefix);

  const int n_momenta = 1000;
  const int n_trials = 100;

  calc_type * p_axis, *growth_rate;
  calc_type min_v = 0.0, max_v = 0.9;
  p_axis = make_momentum_axis(n_momenta, max_v);
  growth_rate = (calc_type*) malloc(n_trials* sizeof(calc_type));

  if(!p_axis || ! growth_rate){
    my_print("Cannot allocate memory for arrays", mpi_info.rank);
    
    if(p_axis) free(p_axis);
    if(growth_rate) free(growth_rate);
    //Free whichever was allocated
    
    safe_exit();
  }

  calc_type omega;
  calc_type d_om = std::abs(my_plas->get_omega_ref("ce")) / (float) n_trials;
  omega = d_om;
  
  std::ofstream outfile;
  outfile.open("Growth.dat");

  if(!outfile){
     my_print("Error opening output file", mpi_info.rank);
     free(p_axis);
     free(growth_rate);
     safe_exit();
  
  }

  if(mpi_info.rank ==0) write_growth_header(cmd_line_args.file_prefix, my_plas, my_elec, n_momenta, min_v, max_v, outfile);
  for(int i=0; i<n_trials; ++i){
    
    growth_rate[i] = get_growth_rate(my_plas, my_elec, n_momenta, p_axis, omega);
//    std::cout<<my_plas<<" "<<my_elec<<std::endl;
//    std::cout<<n_momenta<<" "<<omega<<std::endl;
//     std::cout<<p_axis<<std::endl;
    if(mpi_info.rank ==0) write_growth(omega, growth_rate[i], outfile);
    omega += d_om;
  }

  outfile.close();

  safe_exit();

}

calc_type get_growth_rate(plasma * my_plas, non_thermal * my_elec, int n_momenta, calc_type * p_axis, calc_type omega_in){

  calc_type k = my_plas->get_dispersion(omega_in, WAVE_WHISTLER, 1);

  calc_type ck_om = v0*k/omega_in, om_ce = std::abs(my_plas->get_omega_ref("ce")), om_pe = std::abs(my_plas->get_omega_ref("pe")), om_diff = omega_in - om_ce;
  
  calc_type f_tmp, a_par, a_perp, v_tmp, norm_f, S_tot=0.0, S_full_tot=0.0, dp, A_crit, A_rel, eta_rel;
  calc_type  *S, *S_full, * dp_ax;
  calc_type gamma, p_res, Delta_res;

  S = (calc_type *)malloc(n_momenta*sizeof(calc_type));
  S_full = (calc_type *)malloc(n_momenta*sizeof(calc_type));
  dp_ax = (calc_type *)malloc(n_momenta*sizeof(calc_type));

  
  //Get RMS momenta from velocities...
  // a_x = RMS p_x (Note factor of 2 in perp, not in par...)
  v_tmp = my_elec->v_par;
  //std::cout<<std::sqrt(1.0 - (v_tmp/v0)*(v_tmp/v0))<<std::endl;
  a_par = 2.0*v_tmp / std::sqrt(1.0 - (v_tmp/v0)*(v_tmp/v0));
  //std::cout<<a_par<<std::endl;

  v_tmp = my_elec->v_perp;
  //std::cout<<v_tmp<<std::endl;
  //std::cout<<std::sqrt(1.0 - (v_tmp/v0)*(v_tmp/v0))<<std::endl;

  a_perp = v_tmp / std::sqrt(1.0 - (v_tmp/v0)*(v_tmp/v0));
  //These aren't right for high gamma, but nor is a damn Maxwellian...
//  std::cout<<a_perp<<std::endl;

  norm_f = 1.0/(a_perp*a_perp*a_par * pi * std::sqrt(pi));
  
  //std::cout<<norm_f<<std::endl;
  
  for(int j=0; j< n_momenta; ++j){
  
    gamma = - 1.0 + ck_om * std::sqrt( (ck_om*ck_om -1.0 )*(1.0 + p_axis[j]*p_axis[j]/v0/v0)*(omega_in*omega_in/om_ce/om_ce) + 1 );
    gamma /= ((ck_om*ck_om -1)*omega_in/om_ce);
    //14 in Xiao resonant gamma factor
//    std::cout<<gamma[j]<<std::endl;
    p_res = (gamma * om_diff)/k;
//    std::cout<<p_res[j]<<std::endl;

    //Resonant momentum
    Delta_res = 1.0 - (omega_in*p_res / (v0*v0*k*gamma));
    //Xiao 15, no meaning given
    
    //For f Maxwellian as Xiao 28: d f/ dp_x = 2 p_x a_x
    //Now p_par = p_res and p_perp is p_axis[j]
    f_tmp = norm_f * std::exp(- (p_res*p_res/(a_par*a_par)) - (p_axis[j]*p_axis[j]/(a_perp*a_perp)));
    S[j] = std::pow(p_axis[j], 3) * f_tmp / Delta_res;
    //Arrays in case p_axis is not uniform...
    S_full[j] = (om_ce/gamma - omega_in)*S[j];
    
//    S_tot += S[j];
  //  S_full_tot += S_full[j];
  
  }

  dp_ax[0] = 0.0;
  for(int j=1; j<n_momenta; j++) dp_ax[j] = p_axis[j] - p_axis[j-1];

  //dp = std::abs(p_axis[0] - p_axis[1]);
  //S_tot *=dp;
  //S_full_tot *=dp;
  //Assume for now linear grid.

  S_tot = integrator(S, n_momenta, dp_ax);
  S_full_tot = integrator(S_full, n_momenta, dp_ax);

  calc_type ret = 0.0;
  if(std::abs(S_tot) > GEN_PRECISION){
    A_crit = omega_in /om_diff;
    A_rel = ((a_perp*a_perp/(a_par*a_par)) - 1.0)*S_full_tot/om_diff/S_tot;
    eta_rel = -2.0 * pi * my_elec->fraction* om_diff/k * S_tot / a_perp/a_perp;
    
    
    ret = pi*om_pe*om_pe/(2.0*omega_in + om_pe*om_pe*om_ce/(std::pow(om_diff, 2))) * eta_rel * (A_rel - A_crit);
  
  }

   std::cout<<ret/om_ce<<std::endl;

  free(S);
  free(S_full);
  free(dp_ax);

  return ret;
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
  if( std::abs(p_axis[n_mom-1] -  (v_max * v0/ std::sqrt(1.0- v_max*v_max)))> 1e-4) my_print("Momentum axis errorrrrrrr", mpi_info.rank);
  
  return p_axis;

}

g_args g_command_line(int argc, char * argv[]){
/** Check whether to process no files (analytic only), sdf or spectrum 0, 1,2, respectively. Looping through again is silly, but we're stealing from main in chunks here... But if we have a -f arg we use sdf files, if a -s we use the spectra listed in that file, if neither, we output analytic only and if both the last one is used*/

  g_args extra_cmd_line;

  extra_cmd_line.src = 0;
  extra_cmd_line.spect_file = "";
  
  for(int i=1; i< argc; i++){
    if(strcmp(argv[i], "-s")==0 && i < argc-1){
      extra_cmd_line.spect_file = atoi(argv[i+1]);
      extra_cmd_line.src = 2;
    }
    if(strcmp(argv[i], "-f")==0 && i < argc-1) extra_cmd_line.src = 1;
  }
  return extra_cmd_line;


}

void write_growth_header(std::string in_file, plasma * my_plas, non_thermal * my_elec, int n_momenta, calc_type min_v, calc_type max_v, std::ofstream &outfile){

  /** Write general parameters for the growth calcs...*/

  outfile<<"Source: "<<in_file<<"\n";
  
  outfile<<"Number of momentum points "<<n_momenta<<"\n";
  outfile<<"Min and max velocity "<<min_v<<" "<<max_v<<"\n";
  outfile<<"Electron parameters: "<<"\n";
  my_elec->write(outfile);
  outfile<<"Omega \t Growth rate:"<<"\n";
  outfile<<"BEGIN"<<"\n";
}

void write_growth(calc_type omega, calc_type growth, std::ofstream &outfile){

  /** Write growth rates*/
  outfile<<omega<<" "<<growth<<"\n";

}

