//
//  d_coeff.cpp
//  
//

#include <math.h>
#include <algorithm>
#include <boost/math/special_functions.hpp>
#include "support.h"
#include "d_coeff.h"

diffusion_coeff::diffusion_coeff(int n_momenta, int n_angs):data_array(n_momenta, n_angs){
/** \brief Create coefficient
*
*Creates rectangular data_array and sets additional parameters to defaults @param nx Number of points in momentum to use @param n_thetas Number of particle pitch angle points
*/
  my_controller = nullptr;

  n_thetas = 100;
  n_n = 15;
  wave_id = WAVE_WHISTLER;
  latitude = 0;
  tag="";
}

void diffusion_coeff::set_ids(float time1, float time2, int space1, int space2, int wave_id, char block_id[ID_SIZE]){
/**\brief Set parameters
*
*Sets the time and space ranges, wave type etc attached to the spectrum. Times should be in terms of file output time. Space in terms of grid points.
*/

  this->time[0] = time1;
  this->time[1] = time2;
  this->space[0] = space1;
  this->space[1] = space2;
  strcpy(this->block_id, block_id);
  this->wave_id = wave_id;
}


bool diffusion_coeff::write_to_file(std::fstream &file){
/** \brief Write diffusion coeff to file
*
* Writes file dump of a diffusion coefficient. NB the file passed in must be opened for input and output. \todo Read routine
*/

  if(!file.is_open()) return 1;
  
  data_array::write_to_file(file);
  file.seekg(-1*sizeof(size_t), std::ios::cur);
  size_t ftr_start = 0, ftr_end=0;

  file.read((char*) &ftr_start, sizeof(size_t));
  //Ftr start is position of footer start.
  file.seekg(-1*sizeof(size_t), std::ios::cur);

  //Now we add the other tags
  //First the wave id, then the first 10 chars of the tag string
  file.write((char*) &wave_id, sizeof(int));

  char buffer[10]="";//Initialise to empty

  size_t n_char = std::min((size_t)10, tag.size());
  strncpy(buffer, tag.c_str(), n_char);
  file.write(buffer, 10*sizeof(char));

  ftr_end = (size_t) file.tellg();
  file.write((char*)&ftr_start, sizeof(size_t));
  //Finish file with position of ftr start
  
  file.seekg(ftr_start);
  file.write((char*) &ftr_end, sizeof(size_t));
  //Go back and correct the next_block entry

/*    std::cout<<"end is "<<ftr_end<<'\n';
  file.seekg(ftr_start);

  int wave;
  file.read((char*) &ftr_start, sizeof(size_t));
  std::cout<<"start is "<<ftr_start<<'\n';

  file.read((char*) &buffer, 10*sizeof(char));
std::cout<<buffer<<'\n';
  file.read((char*) &wave, sizeof(int));
  std::cout<<wave<<'\n';
  file.read((char*) &buffer, 10*sizeof(char));
std::cout<<buffer<<'\n';*/
  
  return 0;

}

void diffusion_coeff::make_velocity_axis(){
/**\brief Set velocity axis
*
*Makes linear axis between V_MIN and V_MAX
*/
  calc_type res = (V_MAX - V_MIN)/(this->get_dims(0)-1);
  long offset = std::abs(V_MIN)/res;
  make_linear_axis(0, res, offset);
}

void diffusion_coeff::make_pitch_axis(){
/**\brief Set pitch angle axis (tan theta? alpha???)
*
*Makes linear axis between ANG_MIN and ANG_MAX
*/

  calc_type res = (ANG_MAX - ANG_MIN)/this->get_dims(1); //To cover range from Min to Max
  long offset = ANG_MIN/res;
  make_linear_axis(1, res, offset);
}

d_report diffusion_coeff::calculate(D_type_spec type_of_D, bool quiet){
/** \brief Calculate D from wave spectrum and plasma
*
*Uses the data available via my_controller to calculate D, the raw diffusion coefficient as function of particle velocity and pitch angle. For more details of the calculation see Derivations#Calculation_of_D Note that here we use a default number of wave-normal angle points, NOT the number present in the parent spectrum. See get_G2 for details of how interpolation is done.
*/

//----Initialise d_report --------
  d_report report;
  report.n_solutions = 0;
  report.n_fails = 0;
  report.error = true;
  report.n_av = 0;
  report.n_max = 0;
  report.n_min = 0;
  
//----- Get the plasma and spectrum bits
  plasma plas;
  spectrum * spect;
  if(my_controller){
    plas = my_controller->get_plasma();
    spect = my_controller->get_current_spectrum();
  }
  else{
    my_error_print("No controller", mpi_info.rank);
    return report;
  }
  this->copy_ids(spect);//Copy block id, ranges etc from spect.
  calc_type k_thresh = spect->check_upper();
  my_print("Using omega upper threshold of "+ mk_str(plas.get_dispersion(k_thresh, WAVE_WHISTLER)/plas.get_omega_ref("ce")), 1);

//---- Set type of D to calculate
  bool is_mixed=false, is_pp=false;
  if(type_of_D == D_type_spec::alpha_p || type_of_D == D_type_spec::p_alpha) is_mixed = true;
  else if(type_of_D == D_type_spec::p_p) is_pp = true;

//----- Major definitions------------------
  calc_type theta, omega_n = 0.0, D_part_n_sum, D_final_without_consts;
  calc_type alpha, mod_v, c2th, s2alpha, tan_alpha, cos_alpha;
  calc_type Eq6, numerator, gamma_particle, D_conversion_factor;
  int n_min, n_max;
  std::vector<calc_type> omega_calc;
  mu_dmudom my_mu;

  calc_type om_ce_ref = plas.get_omega_ref("ce"), omega_solution;
  calc_type D_consts = 0.5* pi*om_ce_ref*om_ce_ref*v0 * std::pow(1.0/plas.get_B0(), 2) * spect->get_norm_B();//Constant part of D ( pi/2 Om_ce*c^3) * m^2 / B_0^2. But we then divide by (m^2 c^2) so that we get D/p^2 by just dividing by (gamma^2 - 1)

//-------- Wave angle temporaries -----------------
  //D has to be integrated over x, so we need a temporary and the axis
  calc_type *D_theta = (calc_type *) calloc(n_thetas, sizeof(calc_type));
  calc_type *x = (calc_type *) calloc(n_thetas, sizeof(calc_type));
  calc_type *dx = (calc_type *) calloc(n_thetas, sizeof(calc_type));
  for(int i = 0; i < n_thetas; ++i) x[i] = i* TAN_MAX/n_thetas;
  for(int i = 0; i < n_thetas - 1; ++i) dx[i] = x[i+1] - x[i];
  

//-------- Some reporting info --------
  size_t last_report = 0;
  size_t report_interval = dims[0]/10; //10 prints per round
  if(report_interval > 20) report_interval = 20;
  if(report_interval < 1) report_interval = 1;

  int counter = 0, non_counter = 0;
  size_t n_av = 0, n_omega=0;

//-------------------Main loops here----------------------------
//We have deep nested loops. Move ANYTHING that can be as far up tree as possible

  for(size_t v_ind = 0; v_ind < dims[0]; v_ind++){
    //particle parallel velocity or momentum
    mod_v = get_axis_element(0, v_ind);
    gamma_particle = gamma_rel(mod_v);
    if(std::abs(mod_v) < 1.0) continue;//Skip if velocity is tiny
    
    if(!quiet){
      my_print("Velocity "+mk_str(mod_v/v0, true)+" c", mpi_info.rank);
      if((v_ind - last_report) >= report_interval){
        my_print("i "+mk_str(v_ind)+" (n_max, n_min "+mk_str(report.n_max)+' '+mk_str(report.n_min)+")", mpi_info.rank);
        last_report = v_ind;
      }
    }
    for(size_t part_pitch_ind = 0; part_pitch_ind < dims[1]; part_pitch_ind++){
      //particle pitch angle
      alpha = get_axis_element_ang(part_pitch_ind);
      tan_alpha = tan(alpha);
      cos_alpha = cos(alpha);
      s2alpha = std::pow(std::sin(alpha), 2);

      //Get limits on n for this velocity and angle
      n_min = get_min_n(mod_v, cos_alpha, k_thresh, om_ce_ref);
      n_max = get_max_n(mod_v, cos_alpha, k_thresh, om_ce_ref);
      n_av += n_max;//Track the average n_max over all iterations
      report.n_max = std::max(report.n_max, (size_t)n_max);//Track the extreme n_max and n_min
      report.n_min = std::max(report.n_min, (size_t) -n_min);

      for(int wave_ang_ind = 0; wave_ang_ind < n_thetas; wave_ang_ind++){
      //theta loop for wave angle or x=tan theta
        theta = atan(x[wave_ang_ind]);
        c2th = std::pow(cos(theta), 2);
        D_part_n_sum = 0.0;

        for(int n = n_min; n < n_max; ++n){
//          {int n = 0;
          // n is resonant number
          omega_calc = plas.get_resonant_omega(x[wave_ang_ind], mod_v*cos_alpha, gamma_particle, (calc_type) n);
          omega_n = (calc_type) n * om_ce_ref;
          n_omega = omega_calc.size();
          if(n_omega > 0){
          //With loop if is redundant but we might want to log the failure
            for(size_t om_solution_num = 0; om_solution_num < n_omega; ++om_solution_num){
              //Loop over solutions
              //Ignore sign info for now. This might lead to double or quad counting, that can be resolved later
              /** \todo Check for double counting*/
              omega_solution = std::abs(omega_calc[om_solution_num]);
              
              my_mu = plas.get_high_dens_phi_mu_om(omega_solution, theta, alpha, n, gamma_particle);

              my_mu.err ? non_counter++ : counter++;
              //If err is true, then once angle is included we found no solution. Keep track for final report

              Eq6 = omega_solution/(omega_solution - omega_n)* my_mu.mu/(my_mu.mu + omega_solution * my_mu.dmudom);
              numerator = std::pow( -s2alpha + omega_n/omega_solution, 2);
              D_conversion_factor = sin(alpha)*cos(alpha)/(-s2alpha + omega_n/omega_solution);
              //Get the part of D summed over n. Note additional factors below
              //NB NB get_G_1 precancels the \Delta\omega and B_wave factors
              D_part_n_sum += numerator * my_mu.phi / std::abs(1.0 - Eq6) * get_G1(spect, omega_solution) * get_G2(spect, omega_solution, x[wave_ang_ind]);
              //Convert alpha_alpha to requested D type
              D_part_n_sum = (is_mixed ? D_part_n_sum*D_conversion_factor : D_part_n_sum);
              D_part_n_sum = (is_pp ? D_part_n_sum*D_conversion_factor*D_conversion_factor : D_part_n_sum);
            }
          }
        }
        //Store into temporary theta array
        D_theta[wave_ang_ind] = D_part_n_sum * c2th * x[wave_ang_ind];
      }
      //now integrate in x = tan theta, restore the velocity factor and save
      D_final_without_consts = integrator(D_theta, n_thetas, dx);
      D_final_without_consts /= mod_v*std::pow(cos_alpha, 3);
      set_element(v_ind, part_pitch_ind, D_final_without_consts*D_consts);
    }
  }
//------------End main loops-----------------------------------

  free(dx);
  free(x);
  free(D_theta);
  
//----Assemble the final report -----------------------------
  report.n_solutions = counter;
  report.n_fails = non_counter;
  report.error = false;
  report.n_av = n_av/dims[0]/dims[1];
  my_print(mk_str(counter) + " solutions vs "+ mk_str(non_counter), mpi_info.rank);
  
  return report;
  
}

int diffusion_coeff::get_min_n(calc_type mod_v, calc_type cos_alpha,  my_type k_thresh, calc_type om_ce){
/** \brief Limits on n 
*
* Uses the maximum k and the velocity to give min/max n to consider (note signs). Note always has abs value ge 1.
*/

  calc_type gamma = gamma_rel(mod_v);
  int n = std::max(-(int)(gamma * k_thresh * std::abs(mod_v*cos_alpha / om_ce)), -n_n);
  return std::min(-1, n-1);
}

int diffusion_coeff::get_max_n(calc_type mod_v, calc_type cos_alpha, my_type k_thresh, calc_type om_ce){
/** \brief Limits on n 
*
* Uses the maximum k and the velocity to give min/max n to consider (note signs) and cap to 1.
*/

  calc_type gamma = gamma_rel(mod_v);

  int n = std::min((int)(gamma * k_thresh * std::abs(mod_v*cos_alpha / om_ce)), n_n);
  return std::max(1, n+1);
}

void diffusion_coeff::copy_ids( spectrum * src){
/** Copies ID fields from src array to this*/

  strcpy(this->block_id, src->block_id);
  std::copy(src->time, src->time + 2, this->time);
  for(int i=0; i < 2; ++i) this->space[i] = src->space[i];
}

my_type diffusion_coeff::get_element_by_values(my_type p, my_type alpha){
/** Get D element by values of p and alpha. \todo Q? try interpolating? */
  if(this->get_dims() != 2){
#ifdef DEBUG_DIMS
    my_error_print("Wrong dimensions, attempting 2 with "+mk_str(this->get_dims()));
#endif
    return 0.0;
  }
  size_t len_a=0, len_p = 0;
  long alpha_ind, p_ind;

  my_type * d_axis = this->get_axis(0, len_p);
  p_ind = where(d_axis, len_p, p);
  d_axis = this->get_axis(1, len_a);
  alpha_ind = where(d_axis, len_a, angle_to_stored_angle(alpha));

  if(p_ind >=0 && alpha_ind >= 0 && p_ind < len_p && alpha_ind < len_a){
    return this->get_element(p_ind, alpha_ind);
  }else{
    return 0.0;
  }

}

my_type diffusion_coeff::get_axis_element_ang(size_t ind){
/** Return ANGLE value, rather than raw axis entry*/
  return stored_angle_to_angle(this->get_axis_element(1, ind));

}

