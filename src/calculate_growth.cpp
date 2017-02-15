/** \file calculate_growth.cpp \brief Helper program to get expected and/or actual growth rates of waves
*
* Can take derived spectrum files, extract the approximate wave growth rates (peak or integrated) (assumes spectra are in time-order and all have the same axes) and output these, plus a theoretical rate, OR just the theoretical rate. In former case we use the last files axes, in latter we use a wide coverage log axis with n_trials elements. Parameters are obtained as in Main, so requires a plasma.conf file and deck.status file. Output uses same format as our arrays etc A set of test arguments is supplied. Supply a file containing paths to spectra to use, in order. Paths should be relative to given -f parameter Call using ./calculate_growth `<growth_test_pars` to use these. Or try ./calculate_growth -h 
  \author Heather Ratcliffe \date 11/02/2016.
*/


#include <math.h>
#include <cmath>
#include <fstream>
//Provides Bessel functions, erf, and many more
#include <iostream>
#include <limits>
#include <stdio.h>
#include <mpi.h>
#include <complex.h>

#include "support.h"
#include "main_support.h"
#include "plasma.h"
#include "my_array.h"
#include "spectrum.h"
#include "non_thermal.h"


/** \defgroup utils Utility programs
*@{ */

/** \defgroup growth_util Growth rate calculation utility
*@{ */

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

bool write_growth_closer(std::string in_file, non_thermal * my_elec, size_t n_momenta, calc_type min_v, calc_type max_v, std::fstream &file);

std::vector<std::string> read_filelist(std::string infile);

void dump_distrib(non_thermal * my_elec, std::string filename);

calc_type estimate_spectrum_noise(spectrum & spec_in);

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
    if(filelist.size() >= 1){

      spectrum * my_spect = new spectrum(cmd_line_args.file_prefix+ filelist[0]);
      if(my_spect->get_B_dims() < 1){
      //Wrap in quotes so we highlight stray trailing whitespace etc
        my_error_print("Invalid or missing spectrum file '"+cmd_line_args.file_prefix+filelist[0]+"'");
      }
      //We read spetra which may be either E or B, so call it F
      data_array F_prev, F_curr;
      //Read and ln F for 0th step
      std::function<calc_type(calc_type)> log_function = [](calc_type el) -> calc_type { return log(el); } ;
      F_prev = my_spect->copy_out_B();
      F_prev.apply(log_function);
      my_type delta_t = 0.0;
      
      for(size_t i=0; i<filelist.size(); i++){
        if(i > 0) my_print("Differencing "+filelist[i-1]+" and "+filelist[i]);
        else my_print("Differencing noise estimate and "+filelist[i]);
        delete my_spect;
        my_spect = new spectrum(cmd_line_args.file_prefix+filelist[i]);
        if(my_spect->get_B_dims() < 1){
          my_error_print("Invalid or missing spectrum file '"+cmd_line_args.file_prefix+filelist[i]+"'");
          continue;
        }else if(my_spect->get_B_dims(0) < 3){
          //We use first few els in special case below, and < 3 els will never be useful growth rate
          my_error_print("Invalid field in '"+cmd_line_args.file_prefix+filelist[i]+"'");
          continue;
        }
        //Read and ln F for this step
        F_curr = my_spect->copy_out_B();
        F_curr.apply(log_function);
        //Now calc gamma from \gamma = \ln (F_2/F_1) / (t_2 - t_1)
        numeric_data = F_curr;
        if(i > 0){
          //Average prev and current time vals to get mid-point and then take difference
          delta_t = (numeric_data.time[1] + numeric_data.time[0] -F_prev.time[1] - F_prev.time[0])/2.0;
          numeric_data.subtract(F_prev);

        }else{
        //We need some sort of estimate of the spectrum "noise" to eliminate this in the first file. Then we assume F(t_0) = ref_val. Perhaps should average raw not logged, but this will do.
//          calc_type ref_val = (F_prev.get_element(6) + F_prev.get_element(7))/2.0;
          calc_type ref_val = estimate_spectrum_noise(*my_spect);
          ref_val = log_function(ref_val);
          my_print("Noise estimate: "+mk_str(ref_val));
          std::function<calc_type(calc_type)> minus_const_function = [ref_val](calc_type el) -> calc_type { return el - ref_val; } ;
          numeric_data.apply(minus_const_function);
          delta_t =(numeric_data.time[1] + numeric_data.time[0])/2.0;
        }
        //Since we're working with F^2 we also divide by 2
        numeric_data.divide(delta_t*2.0);
        
        if(numeric_data.get_dims() > 0) err |= numeric_data.write_to_file(file, false);

        F_prev = F_curr;
      }
      delete my_spect;
    }else{
      my_error_print("Empty or missing spectrum file '"+cmd_line_args.file_prefix+extra_args.spect_file+"'", mpi_info.rank);
    }
  }
  
  //Now we do the analytics.

  my_print("Calculating analytic growth rates", mpi_info.rank);

  non_thermal * my_elec = new non_thermal(cmd_line_args.file_prefix);
  if(my_elec->get_norely()) my_print("Selected non relativistic calculation", mpi_info.rank);

  const size_t n_momenta = 10000;
  calc_type * p_axis;
  calc_type min_v = 0.0, max_v = 0.99;
  p_axis = make_momentum_axis(n_momenta, max_v);
  calc_type growth_rate = 0.0;

  if(!p_axis){
    my_error_print("Cannot allocate memory for arrays", mpi_info.rank);
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

  for(size_t i=0; i<analytic_data.get_dims(0); ++i){
    omega = analytic_data.get_axis_element(0, i);
    growth_rate = get_growth_rate(my_plas, my_elec, n_momenta, p_axis, omega);
    analytic_data.set_element(i, growth_rate);
  }

  err |= analytic_data.write_to_file(file, false);
  
  //Special finish for the file  
  err |=write_growth_closer(cmd_line_args.file_prefix, my_elec, n_momenta, min_v, max_v, file);

  file.close();
  if(!err) my_print("Written in "+filename);
  else my_error_print("Error writing to "+filename);
  dump_distrib(my_elec, cmd_line_args.file_prefix+"NonThermalDump.dat");
  safe_exit();

}

calc_type get_growth_rate(plasma * my_plas, non_thermal * my_elec, int n_momenta, calc_type * p_axis, calc_type omega_in){

  omega_in = std::abs(omega_in);
  calc_type k = my_plas->get_dispersion(omega_in, WAVE_WHISTLER, 1);

  calc_type ck_om = v0*k/omega_in, om_ce = std::abs(my_plas->get_omega_ref("ce")), om_pe = std::abs(my_plas->get_omega_ref("pe")), om_diff = omega_in - om_ce;
  //This corrects omega_pe to match the fast electrons we're using, which may not match the deck exactly
  om_pe = om_pe*std::sqrt(my_elec->get_total_dens()/my_const.dens_factor);
  
  //Calculate k using Eq 9 of Xiao  k_paper = std::sqrt( (std::pow(omega_in, 2) - std::pow(om_pe, 2)*omega_in/om_diff)/std::pow(v0, 2));

  calc_type S_tot=0.0, S_full_tot=0.0, A_crit, A_rel, eta_rel;
  calc_type  *S, *S_full, * dp_ax;
  calc_type gamma, p_res, Delta_res, v_res;
  bool norely = my_elec->get_norely();

  S = (calc_type *)malloc(n_momenta*sizeof(calc_type));
  S_full = (calc_type *)malloc(n_momenta*sizeof(calc_type));
  dp_ax = (calc_type *)malloc(n_momenta*sizeof(calc_type));

  v_res = om_diff/k;
  
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
    if(!norely){
      // Relativistic version---------------------------------
      S[j] = std::pow(p_axis[j], 2) * my_elec->d_f_p(p_res, p_axis[j], 0) / Delta_res;

      S_full[j] = std::pow(p_axis[j], 2)/ Delta_res/gamma * (p_axis[j]*my_elec->d_f_p(p_res, p_axis[j], 1) - p_res*my_elec->d_f_p(p_res, p_axis[j], 0));
    }else{
      //---Non-rely version--------------------------

      S[j] = p_axis[j]*my_elec->f_p(v_res, p_axis[j]);
    
      S_full[j] = - std::pow(p_axis[j], 2)* (p_axis[j]*my_elec->d_f_p(v_res, p_axis[j], 1) - v_res*my_elec->d_f_p(v_res, p_axis[j], 0));
    }
  }

  dp_ax[0] = 0.0;
  for(int j=1; j<n_momenta; j++) dp_ax[j] = p_axis[j] - p_axis[j-1];

  S_tot = integrator(S, n_momenta, dp_ax);
  S_full_tot = integrator(S_full, n_momenta, dp_ax);

  calc_type ret = 0.0;

  if(std::abs(S_tot) > std::numeric_limits<calc_type>::min()){
    A_crit = - omega_in /om_diff;

    if(!norely){
      // Relativistic version---------------------------------
      A_rel = k*S_full_tot/S_tot/om_diff;
      eta_rel = pi * om_diff/k * S_tot;
    }else{
      //---Non-rely version--------------------------
      A_rel = S_full_tot/S_tot/2.0/v_res;
      eta_rel = -2.0*pi * v_res * S_tot;
    }

    ret = pi*om_pe*om_pe/(2.0*omega_in + om_pe*om_pe*om_ce/(std::pow(om_diff, 2))) * eta_rel * (A_rel - A_crit);
  }

  free(S);
  free(S_full);
  free(dp_ax);

  return ret;
}


calc_type * make_momentum_axis(int n_mom, calc_type v_max){

  //Max velocity to consider norm'd to c
  calc_type dp = v_max * v0/ (calc_type) (n_mom - 1)/std::sqrt(1.0- v_max*v_max);
  //Momentum step size, including gamma

  calc_type * p_axis;
  p_axis = (calc_type*) malloc(n_mom*sizeof(calc_type));
  if(!p_axis) return nullptr;

  p_axis[0] = 0.0;
  for(int i=1; i< n_mom; ++i) p_axis[i] = p_axis[i-1] + dp;
  //Set momenta from 0 to v_max. Do this once.
  if( std::abs(p_axis[n_mom-1] -  (v_max * v0/ std::sqrt(1.0- v_max*v_max)))> 1e-3) my_error_print("Momentum axis error of "+mk_str(std::abs(p_axis[n_mom-1] -  (v_max * v0/ std::sqrt(1.0- v_max*v_max)))), mpi_info.rank);
  
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
/** \todo -spec is not an error??? Nor is -in*/
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

bool write_growth_closer(std::string in_file, non_thermal * my_elec, size_t n_momenta, calc_type min_v, calc_type max_v, std::fstream &file){
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
  if(write_err) my_error_print("Error writing offset positions", mpi_info.rank);
  file.write((char*) & ftr_start, sizeof(size_t));

  return write_err;
}

void dump_distrib(non_thermal * my_elec, std::string filename){
  //Write out a table of the current non_thermal distribution
  
  size_t x_len = 500, y_len=500;
  my_type x_ax_max = 0.99*v0, y_ax_max = 0.99*v0;
  std::string block_id = "nonth";

  //Create array and axes
  data_array created_array = data_array(500, 500);
  strncpy(created_array.block_id, block_id.c_str(), ID_SIZE);
  for(size_t i=0; i<x_len; i++) created_array.set_axis_element(0, i, -x_ax_max + (float)i*(x_ax_max*2.0/(float) x_len));
  for(size_t i=0; i<y_len; i++) created_array.set_axis_element(1, i, -y_ax_max + (float)i*(y_ax_max*2.0/(float) y_len));
  
  for(size_t i=0; i< x_len; i++){
    for(size_t j=0; j< y_len; j++){
   // std::cout<<i<<' '<<j<<' '<<created_array.get_axis_element(0, i)<<' '<< created_array.get_axis_element(1, j)<<'\n';
      created_array.set_element(i, j, my_elec->f_p(created_array.get_axis_element(0, i), created_array.get_axis_element(1, j)));
    }
  }
  std::fstream outfile;
  outfile.open(filename,std::ios::binary|std::ios::out);
  if(outfile) created_array.write_to_file(outfile);
}

calc_type estimate_spectrum_noise(spectrum & spec_in){
/** Estimates the "noise" in a given spectrum. 
*
* There are many ways we might do this, using varying amounts of additional information about the "spectrum". Averaging the end values works poorly. Averaging the largest two values that are greater than something? Fit a straight line to the 0.8 to 1 region?
*/
/** \todo Find better noise estimate...*/
  calc_type noise_val = (spec_in.get_B_element(6) + spec_in.get_B_element(7))/2.0;

  return noise_val;
}

/** @} */
/** @} */
