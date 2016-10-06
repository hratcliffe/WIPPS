/** \file calculate_growth.cpp \brief Helper program to get expected and/or actual growth rates of waves
*
* Can take derived spectrum files, extract the approximate wave growth rates (peak or integrated) (assumes spectra are in time-order and all have the same axes) and output these, plus a theoretical rate, OR just the theoretical rate. In former case we use the last files axes, in latter we use a wide coverage log axis with n_trials elements. Parameters are obtained as in Main, so requires a plasma.conf file and deck.status file. Output uses same format as our arrays etc A set of test arguments is supplied. Supply a file containing paths to spectra to use, in order. Paths should be relative to given -f parameter Call using ./calculate_growth `<growth_test_pars` to use these. Or try ./calculate_growth -h
  \author Heather Ratcliffe \date 11/02/2016.
*/


#include <math.h>
#include <cmath>
#include <fstream>
#include <boost/math/special_functions.hpp>
//Provides Bessel functions, erf, and many more
#include <iostream>
#include <limits>
#include <stdio.h>
#include <mpi.h>
#include <complex.h>

#include "support.h"
#include "main_support.h"
#include "controller.h"
#include "plasma.h"
#include "my_array.h"
#include "d_coeff.h"
#include "spectrum.h"
#include "non_thermal.h"


const int n_trials = 2000;/**< Number of data points for analytic (if no real data)*/

deck_constants my_const;/**< Physical constants*/
extern const mpi_info_struc mpi_info;/**< Link to mpi_info as const*/

struct g_args{

  bool real;/**<Whether to derive real growth from spectra list too*/
  std::string spect_file;/**<File listing spectra if using this option*/
  std::string outfile;/**< Output filename*/
};


calc_type * make_momentum_axis(int n_momenta, calc_type v_max_over_c);

calc_type get_growth_rate(plasma * my_plas, non_thermal * my_elec, int n_momenta, calc_type * p_axis, calc_type omega_in);
g_args g_command_line(int argc, char * argv[]);

bool write_growth_closer(std::string in_file, plasma * my_plas, non_thermal * my_elec, size_t n_momenta, calc_type min_v, calc_type max_v, std::fstream &file);

std::vector<std::string> read_filelist(std::string infile);

int main(int argc, char *argv[]){

  int err=0;
  
  int ierr = local_MPI_setup(argc, argv);
  if(ierr){
    std::cout<< "Error initialising MPI. ABORTING!";
    return 1;
  }

  my_print(std::string("Code Version: ")+ VERSION, mpi_info.rank);
  my_print("Code is running on "+mk_str(mpi_info.n_procs)+" processing elements.", mpi_info.rank);

  g_args extra_args = g_command_line(argc, argv);
  setup_args cmd_line_args = process_command_line(argc, argv);

  if(mpi_info.rank == 0) get_deck_constants(cmd_line_args.file_prefix);
  share_consts();
  /** Get constants from deck and share to other procs*/

  // Setup output file
  std::string filename = cmd_line_args.file_prefix + extra_args.outfile;
  std::fstream file;
  file.open(filename.c_str(),std::ios::out|std::ios::in|std::ios::binary|std::ios::trunc);
  if(!file.is_open()) return 1;

  //Setup other things we need
  plasma * my_plas = new plasma(cmd_line_args.file_prefix);

  data_array numeric_data, analytic_data;

  if(extra_args.real){
    //Read spectrum files
    my_print("Processing spectrum files", mpi_info.rank);
    std::vector<std::string> filelist;
    filelist = read_filelist(cmd_line_args.file_prefix+extra_args.spect_file);
    if(filelist.size() >= 2){

      spectrum * my_spect = new spectrum(cmd_line_args.file_prefix+ filelist[0]);
      if(my_spect->get_B_dims() <1){
      //Wrap in quotes so we highlight stray trailing whitespace etc
        my_print("Invalid or missing spectrum file '"+cmd_line_args.file_prefix+filelist[0]+"'");
      }
      data_array B_prev;
      B_prev = my_spect->copy_out_B();
      my_type delta_t = 0.0;
      for(size_t i=1; i<filelist.size(); i++){
        delete my_spect;
        my_spect = new spectrum(cmd_line_args.file_prefix+filelist[i]);
        if(my_spect->get_B_dims() <1){
          my_print("Invalid or missing spectrum file '"+cmd_line_args.file_prefix+filelist[i]+"'");
          continue;
        }
        numeric_data = my_spect->copy_out_B();
        //Now do difference between spectra
        numeric_data.subtract(B_prev);
        //Average prev and current time vals to get mid-point and then take difference
        delta_t = (numeric_data.time[1] + numeric_data.time[0] -B_prev.time[1] - B_prev.time[0])/2;
        
        numeric_data.divide(delta_t);
        
        if(B_prev.get_dims() > 0 && numeric_data.get_dims() > 0) err |= numeric_data.write_to_file(file, false);

        B_prev = numeric_data;
      }
      delete my_spect;
    }
  }
  
  //Now we do the analytics.

  my_print("Calculating analytic growth rates", mpi_info.rank);

  non_thermal * my_elec = new non_thermal(cmd_line_args.file_prefix);

  const size_t n_momenta = 10000;
  calc_type * p_axis;
  calc_type min_v = 0.0, max_v = 0.95;
  p_axis = make_momentum_axis(n_momenta, max_v);
  calc_type growth_rate = 0.0;

  if(!p_axis){
    my_print("Cannot allocate memory for arrays", mpi_info.rank);
    if(p_axis) free(p_axis);
    safe_exit();
  }

  calc_type omega;
  //Setup an array with desired omega axis
  if(!extra_args.real){
    //Log axis over 3 orders of magnitude
    calc_type d_om = std::abs(my_plas->get_omega_ref("ce")) / (float) (n_trials-1);
    omega = 0.0;
    int om_orders = 3;
    calc_type d_i = (calc_type) (n_trials -1)/(calc_type) om_orders;
    //Orders of magnitude to cover
    d_om = std::abs(my_plas->get_omega_ref("ce"))/ std::pow(10, om_orders);
    analytic_data = data_array(n_trials);
    for(int i=0; i<n_trials; ++i){
      omega = std::pow(10, (calc_type) i /d_i)*d_om;
      analytic_data.set_axis_element(0, i, omega);
    }
  }else{
    //Copy axis from the last numeric result
    size_t sz = numeric_data.get_dims(0);
    analytic_data = data_array(sz);
    for(size_t i=0; i< sz; i++) analytic_data.set_axis_element(0, i, numeric_data.get_axis_element(0, i));
  
  }
  strcpy(analytic_data.block_id, "g_an");
  //for(int j=1; j<n_momenta; j++) std::cout<<my_elec->f_p(p_axis[j], 0)<<", ";

  for(int i=0; i<analytic_data.get_dims(0); ++i){
    omega = analytic_data.get_axis_element(0, i);
    growth_rate = get_growth_rate(my_plas, my_elec, n_momenta, p_axis, omega);
    analytic_data.set_element(i, growth_rate);
  }

  err |= analytic_data.write_to_file(file, false);
  
  //Special finish for the file  
  err |=write_growth_closer(cmd_line_args.file_prefix, my_plas, my_elec, n_momenta, min_v, max_v, file);

  file.close();
  if(!err) my_print("Written in "+filename);
  else my_print("Error writing to "+filename);
  safe_exit();

}

calc_type get_growth_rate(plasma * my_plas, non_thermal * my_elec, int n_momenta, calc_type * p_axis, calc_type omega_in){

  calc_type k = my_plas->get_dispersion(omega_in, WAVE_WHISTLER, 1), k_paper;

  calc_type ck_om = v0*k/omega_in, om_ce = std::abs(my_plas->get_omega_ref("ce")), om_pe = std::abs(my_plas->get_omega_ref("pe")), om_diff = omega_in - om_ce;
  
  
  //Calculate k using Eq 9 of Xiao k_paper = std::sqrt( (std::pow(omega_in, 2) - std::pow(om_pe, 2)*omega_in/om_diff)/std::pow(v0, 2));
  
  calc_type S_tot=0.0, S_full_tot=0.0, dp, A_crit, A_rel, eta_rel;
  calc_type  *S, *S_full, * dp_ax;
  calc_type gamma, p_res, Delta_res;

  S = (calc_type *)malloc(n_momenta*sizeof(calc_type));
  S_full = (calc_type *)malloc(n_momenta*sizeof(calc_type));
  dp_ax = (calc_type *)malloc(n_momenta*sizeof(calc_type));

  
  //Get RMS momenta from velocities...
  // a_x = RMS p_x (Note factor of 2 in perp, not in par...)

  for(int j=0; j< n_momenta; ++j){
  
    gamma = - 1.0 + ck_om * std::sqrt( (ck_om*ck_om -1.0 )*(1.0 + p_axis[j]*p_axis[j]/v0/v0)*(omega_in*omega_in/om_ce/om_ce) + 1.0 );
    gamma /= ((ck_om*ck_om - 1.0)*omega_in/om_ce);
    //14 in Xiao resonant gamma factor

    p_res = (gamma * omega_in - om_ce)/k;
    //Resonant momentum
    Delta_res = 1.0 - (omega_in*p_res / (v0*v0*k*gamma));
    //Xiao 15, no meaning given. Always +ve
    if(Delta_res < GEN_PRECISION) std::cout<<"ERROR!!"<<std::endl;
    
    //Now p_par = p_res and p_perp is p_axis[j]
    
    S[j] = std::pow(p_axis[j], 2) * my_elec->d_f_p(p_res, p_axis[j], 0) / Delta_res;

    S_full[j] = - std::pow(p_axis[j], 2)/ Delta_res*(omega_in - om_ce/gamma) * (p_res*my_elec->d_f_p(p_res, p_axis[j], 0) - p_axis[j]*my_elec->d_f_p(p_res, p_axis[j], 1))/p_res;

  }

  dp_ax[0] = 0.0;
  for(int j=1; j<n_momenta; j++) dp_ax[j] = p_axis[j] - p_axis[j-1];

  S_tot = integrator(S, n_momenta, dp_ax);
  S_full_tot = integrator(S_full, n_momenta, dp_ax);

  calc_type ret = 0.0;
  
  if(std::abs(S_tot) > std::numeric_limits<calc_type>::min()){
    A_crit = - omega_in /om_diff;
    A_rel = S_full_tot/S_tot/om_diff;
    
    eta_rel = pi * om_diff/k * S_tot;
    
    ret = pi*om_pe*om_pe/(2.0*omega_in + om_pe*om_pe*om_ce/(std::pow(om_diff, 2))) * eta_rel * (A_rel - A_crit);
  
  }

  free(S);
  free(S_full);
  free(dp_ax);

  return ret;
}


calc_type * make_momentum_axis(int n_mom, calc_type v_max){
/** \todo This is v????*/
  //Max velocity to consider norm'd to c
  calc_type dp = v_max * v0/ std::sqrt(1.0- v_max*v_max) / (calc_type) (n_mom - 1);
  //Momentum step size, including gamma

  calc_type * p_axis;
  p_axis = (calc_type*) malloc(n_mom*sizeof(calc_type));
  if(!p_axis) return nullptr;

  p_axis[0] = 0.0;
  for(int i=1; i< n_mom; ++i) p_axis[i] = p_axis[i-1] + dp;
  //Set momenta from 0 to v_max. Do this once.
  if( std::abs(p_axis[n_mom-1] -  (v_max * v0/ std::sqrt(1.0- v_max*v_max)))> 1e-4) my_print("Momentum axis error", mpi_info.rank);
  
  return p_axis;

}

std::vector<std::string> read_filelist(std::string infile){
/** Read file
*
*Reads a list of strings from given infile
*/
  std::cout<<infile<<'\n';
  std::ifstream file;
  file.open(infile);
  std::vector<std::string> names;
  std::string name;

  while(file){
    //Copy all lines that are not just whitespace into vector
    std::getline(file, name);
    if(name.find_first_not_of("\t\n ") !=std::string::npos) names.push_back(name);
  }

  file.close();
  return names;
}

g_args g_command_line(int argc, char * argv[]){
/** Check whether to handle real spectra. If -s we use the spectra (plural) listed in that file, if not we output analytic only*/

  g_args extra_cmd_line;

  extra_cmd_line.real = false;
  extra_cmd_line.spect_file = "";
  extra_cmd_line.outfile = "growth.dat";
  
  for(int i=1; i< argc; i++){
    if(strcmp(argv[i], "-s")==0 && i < argc-1){
      extra_cmd_line.spect_file = argv[i+1];
      extra_cmd_line.real = true;
      strcpy(argv[i], HANDLED_ARG);
      strcpy(argv[i+1], HANDLED_ARG);
      i++;
    }
    else if(strcmp(argv[i], "-out")==0 && i < argc-1){
      extra_cmd_line.outfile = argv[i+1];
      strcpy(argv[i], HANDLED_ARG);
      strcpy(argv[i+1], HANDLED_ARG);
      i++;
    }
    else if(strcmp(argv[i], "-h")==0){
      strcpy(argv[i], HANDLED_ARG);
      print_help('w');
      exit(0);
    }
  }
  return extra_cmd_line;

}

bool write_growth_closer(std::string in_file, plasma * my_plas, non_thermal * my_elec, size_t n_momenta, calc_type min_v, calc_type max_v, std::fstream &file){
  /** Write the general parameters as a file footer with offsets*/
  
  
  size_t ftr_start = file.tellg();
  bool write_err =false;
  size_t next_location = ftr_start + (in_file.size()+1)*sizeof(char)+sizeof(size_t)*3 +sizeof(calc_type)*6;
  //Nonthermal dumps 4 calc types

  file.write((char*) & next_location, sizeof(size_t));
  //Position of next section
  file.write((char*)& n_momenta, sizeof(size_t));
  file.write((char*) &min_v, sizeof(calc_type));
  file.write((char*) &max_v, sizeof(calc_type));
  
  my_elec->dump(file);
  
  size_t len =in_file.size()+1;
  file.write((char*)&len, sizeof(size_t));
  file.write(in_file.c_str(), len*sizeof(char));

  if((size_t)file.tellg() != next_location) write_err=1;
  if(write_err) my_print("Error writing offset positions", mpi_info.rank);
  file.write((char*) & ftr_start, sizeof(size_t));

  return write_err;
}
