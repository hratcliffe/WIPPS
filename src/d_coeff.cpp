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
*Creates rectangular data_array and sets additional parameters to defaults 
@param n_momenta Number of points in momentum to use 
@param n_angs Number of particle pitch angle points
\caveat There is a hard limit on the number of resonances calculated, given by diffusion_coeff::n_n and set here
*/
  my_controller = nullptr;

  n_thetas = 100;
  n_n = 15;
  single_n = false;
  n_used = 0;
  wave_id = WAVE_WHISTLER;
  latitude = 0; //Doesn't do anything "yet"
  tag="";
}

void diffusion_coeff::set_ids(float time1, float time2, int space1, int space2, int wave_id, char block_id[ID_SIZE]){
/**\brief Set parameters
*
*Sets the time and space ranges, wave type etc attached to the spectrum used to calculate this D. Times in seconds. Space in terms of grid points.
@param time1 Initial time of data used
@param time2 End time of data used
@param space1 Start index of space range of data used
@param space2 End index of space range of data used
@param wave_id Wave type (see support.h)
@param block_id Name of block used to calculate this D
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
* Writes file dump of a diffusion coefficient. NB the file passed in must be opened for input and output.
@param file Filestream to write to
@return 0 for success, 1 for failure (usually file access)
/todo Write n_info or perhaps whole report?
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
  
  return 0;
}

bool diffusion_coeff::read_from_file(std::fstream &file){
/** \brief Read diffusion coeff from file
*
* Reads file dump of a diffusion coefficient. D should have been created to the correct size already
@param file Filestream to read from
@return 0 for success, 1 for failure (file access problem)
*/


  if(!file.is_open()) return 1;
  
  bool err = data_array::read_from_file(file);
  //If error we very likely can't continue, so don't try
  if(!err){
    //Grab wave id
    file.read((char*) &wave_id, sizeof(int));
    //Grab tag
    char buffer[10]="";//Initialise to empty
    file.read(buffer, sizeof(char)*10);
    this->tag = buffer;
  }
  return err;

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
/**\brief Set particle pitch angle axis
*
*Makes linear axis between ANG_MIN and ANG_MAX. In equations this is usally called alpha
*/

  calc_type res = (ANG_MAX - ANG_MIN)/this->get_dims(1); //To cover range from Min to Max
  long offset = ANG_MIN/res;
  make_linear_axis(1, res, offset);
}

void diffusion_coeff::first_running_report(size_t total_its, running_report &rept, bool quiet){
/** \brief Initialise progress info*/

  size_t report_interval = total_its/10;
  if(report_interval > 20) report_interval = 20;
  if(report_interval < 1) report_interval = 1;
  rept.last_report = 0;
  rept.report_interval = report_interval;
  rept.quiet = quiet;
}

void diffusion_coeff::do_running_report(size_t v_ind, calc_type mod_v, size_t n_min, size_t n_max, running_report &rept){
/** \brief Print progress info*/

  if(!rept.quiet){
    my_print("Velocity "+mk_str(mod_v/v0, true)+" c", mpi_info.rank);
    if((v_ind - rept.last_report) >= rept.report_interval){
      if(single_n){
        my_print("Calculating velocity number "+mk_str(v_ind)+" (n "+mk_str(n_used)+")", mpi_info.rank);
      }else{
        my_print("Calculating velocity number "+mk_str(v_ind)+" (n_max, n_min "+mk_str(n_max)+' '+mk_str(n_min)+")", mpi_info.rank);
      }
      rept.last_report = v_ind;
    }
  }
}

d_report diffusion_coeff::calculate(D_type_spec type_of_D, bool quiet){
/** \brief Calculate D from wave spectrum and plasma
*
*Uses the data available via my_controller to calculate D, the raw diffusion coefficient as function of particle velocity and pitch angle. For more details of the calculation see Derivations#Calculation_of_D Note that here we use a default number of wave-normal angle points, NOT the number present in the parent spectrum. See get_G2 for details of how interpolation is done.
@param type_of_D Which D to calculate (pitch angle, momentum, mixed)
@param quiet True to suppress display of progress
@return d_report structure containing info on calcs
\caveat We re-range the D calculation so what is actually calculated is D/p^2
\caveat We currently use plasma::get_high_dens_phi_mu_om which uses the high-density approximation to the plasma dispersion, as does the resonant frequency solver plasma::get_resonant_omega. This missed some solutions for smaller om_pe/om_ce. Above 3 or so seems to be alright
*/

//----Initialise d_report --------
  d_report report = {true, 0, 0, 0, false};
  
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
  bool is_mixed = false, is_pp = false;
  if(type_of_D == D_type_spec::alpha_p || type_of_D == D_type_spec::p_alpha) is_mixed = true;
  else if(type_of_D == D_type_spec::p_p) is_pp = true;

//----- Major definitions------------------
  //Velocity and angle temporaries
  calc_type tan_theta, theta, c2th, mod_v, gamma_particle, alpha, s2alpha, cos_alpha;
  //Temporaries for steps of D calculation
  calc_type D_part_n_sum, D_final_without_consts, D_conversion_factor;
  //More temporaries for parts of D
  calc_type Eq6, numerator;
  //Resonant freq and array of all frequency solutions
  calc_type omega_n = 0.0;
  std::vector<calc_type> omega_calc;
  //Resonant numbers
  int n_min, n_max;
  //Plasma dispersion
  mu_dmudom my_mu;

  calc_type om_ce_ref = plas.get_omega_ref("ce"), omega_solution;
  calc_type D_consts = 0.5* pi*om_ce_ref*om_ce_ref*v0 * std::pow(1.0/plas.get_B0(), 2) * spect->get_norm_B();//Constant part of D ( pi/2 Om_ce*c^3) * m^2 / B_0^2. But we then divide by (m^2 c^2) so that we get D/p^2 by just dividing by (gamma^2 - 1) which is done in the last step

//-------- Wave angle temporaries -----------------
  //D has to be integrated over x, so we need a temporary and the axis
  calc_type *D_theta = (calc_type *) calloc(n_thetas, sizeof(calc_type));
  calc_type *x = (calc_type *) calloc(n_thetas, sizeof(calc_type));
  calc_type *dx = (calc_type *) calloc(n_thetas, sizeof(calc_type));
  for(int i = 0; i < n_thetas; ++i) x[i] = i* TAN_MAX/n_thetas;
  for(int i = 0; i < n_thetas - 1; ++i) dx[i] = x[i+1] - x[i];
  
//-------- Initialise reporting info --------
  running_report running_rept = {0, 0, false};
  first_running_report(dims[0], running_rept, quiet);
  size_t n_av = 0, n_omega = 0;

//------------------Main loops here-------------------------
//We have deep nested loops. Move ANYTHING that can be as far up tree as possible

  for(size_t v_ind = 0; v_ind < dims[0]; v_ind++){
    //particle parallel velocity or momentum
    mod_v = get_axis_element(0, v_ind);
    gamma_particle = gamma_rel(mod_v);
    if(std::abs(mod_v) < 1.0) continue;//Skip if velocity is tiny
    //Print the progress info
    do_running_report(v_ind, mod_v, report.n_min, report.n_max, running_rept);

    for(size_t part_pitch_ind = 0; part_pitch_ind < dims[1]; part_pitch_ind++){
      //particle pitch angle
      alpha = get_axis_element_ang(part_pitch_ind);
      cos_alpha = std::cos(alpha);
      s2alpha = std::pow(std::sin(alpha), 2);
      if(single_n){
        //Use only the single selected resonance
        n_min = n_used;
        n_max = n_used;
      }else{
        //Get limits on n for this velocity and angle
        n_min = get_min_n(mod_v, cos_alpha, k_thresh, om_ce_ref);
        n_max = get_max_n(mod_v, cos_alpha, k_thresh, om_ce_ref);
        n_av += n_max;//Track the average n_max over all iterations
        report.n_max = std::max(report.n_max, (size_t)n_max);//Track the extreme n_max and n_min
        report.n_min = std::max(report.n_min, (size_t) -n_min);
      }
      for(int wave_ang_ind = 0; wave_ang_ind < n_thetas; wave_ang_ind++){
      //theta loop for wave angle or x=tan theta
        tan_theta = x[wave_ang_ind];
        theta = atan(tan_theta);
        c2th = std::pow(cos(theta), 2);
        D_part_n_sum = 0.0;

        for(int n = n_min; n <= n_max; ++n){
          // n is resonant number
          omega_calc = plas.get_resonant_omega(theta, mod_v*cos_alpha, gamma_particle, (calc_type) n);
          omega_n = (calc_type) n * om_ce_ref;
          n_omega = omega_calc.size();
          for(size_t om_solution_num = 0; om_solution_num < n_omega; ++om_solution_num){
            //Loop over solutions
            /** \todo Check for double counting*/
            omega_solution = omega_calc[om_solution_num];
            my_mu = plas.get_high_dens_phi_mu_om(omega_solution, theta, alpha, n, gamma_particle);

            Eq6 = omega_solution/(omega_solution - omega_n)* my_mu.mu/(my_mu.mu + omega_solution * my_mu.dmudom);
            numerator = std::pow( -s2alpha + omega_n/omega_solution, 2);
            D_conversion_factor = sin(alpha)*cos(alpha)/(-s2alpha + omega_n/omega_solution);
            //Get the part of D summed over n. Note additional factors below
            //NB NB get_G_1 precancels the \Delta\omega and B_wave factors
            D_part_n_sum += numerator * my_mu.phi / std::abs(1.0 - Eq6) * get_G1(spect, omega_solution) * get_G2(spect, omega_solution, tan_theta);
            
            //Convert alpha_alpha to requested D type
            D_part_n_sum = (is_mixed ? D_part_n_sum*D_conversion_factor : D_part_n_sum);
            D_part_n_sum = (is_pp ? D_part_n_sum*D_conversion_factor*D_conversion_factor : D_part_n_sum);
          }
        }
        //Store into temporary theta array
        D_theta[wave_ang_ind] = D_part_n_sum * c2th * tan_theta;
      }
      //now integrate in x = tan theta, restore the velocity factor and save
      D_final_without_consts = integrator(D_theta, n_thetas, dx);
      D_final_without_consts /= mod_v*std::pow(cos_alpha, 3) * (gamma_particle*gamma_particle - 1.0);
      set_element(v_ind, part_pitch_ind, D_final_without_consts*D_consts);
    }
  }
//------------End main loops-----------------------------------

  free(dx);
  free(x);
  free(D_theta);
  
//----Assemble the final report -----------------------------
  report.error = false;
  report.n_av = n_av/dims[0]/dims[1];
  if(single_n){
    report.n_av = n_used;
    report.n_max = 0;
    report.n_min = 0;
    report.single_n = true;
  }

  return report;
  
}

int diffusion_coeff::get_min_n(calc_type mod_v, calc_type cos_alpha,  my_type k_thresh, calc_type om_ce){
/** \brief Limits on n 
*
* Uses the maximum k and the velocity to give min/max n to consider (note signs).
@param mod_v Abs. value of particle velocity (real value)
@param cos_alpha Cosine of particle pitch angle
@param k_thresh Maximum wavevector where spectrum has non-nelgible power
@param om_ce Reference cyclotron frequency
@return Maximum resonant number, bounded as < -1
*/

  calc_type gamma = gamma_rel(mod_v);
  int n = std::max(-(int)(gamma * k_thresh * std::abs(mod_v*cos_alpha / om_ce)), -n_n);
  return std::min(-1, n);
}

int diffusion_coeff::get_max_n(calc_type mod_v, calc_type cos_alpha, my_type k_thresh, calc_type om_ce){
/** \brief Limits on n 
*
* Uses the maximum k and the velocity to give min/max n to consider (note signs)
@param mod_v Abs. value of particle velocity (real value)
@param cos_alpha Cosine of particle pitch angle
@param k_thresh Maximum wavevector where spectrum has non-nelgible power
@param om_ce Reference cyclotron frequency
@return Maximum resonant number, always at least 1
*/

  calc_type gamma = gamma_rel(mod_v);

  int n = std::min((int)(gamma * k_thresh * std::abs(mod_v*cos_alpha / om_ce)), n_n);
  return std::max(1, n);
}

void diffusion_coeff::copy_ids( spectrum * src){
/** Copies ID fields from src array to this
*@param src Source array to copy from
*/

  strcpy(this->block_id, src->block_id);
  std::copy(src->time, src->time + 2, this->time);
  for(int i=0; i < 2; ++i) this->space[i] = src->space[i];
}

my_type diffusion_coeff::get_element_by_values(my_type p, my_type alpha){
/**  \brief Lookup D element using axis values
*
*
Get D element by values of p and alpha by looking up those values in the axes and returning the nearest value on grid.
@param p Momentum value
@param alpha Angle value
@return Diffusion coeff. value at location
\todo Q? try interpolating?
*/
  if(this->get_dims() != 2){
#ifdef DEBUG_DIMS
    my_error_print("Wrong dimensions, attempting 2 with "+mk_str(this->get_dims()));
#endif
    return 0.0;
  }
  size_t len_a = 0, len_p = 0;
  long alpha_ind, p_ind;

  const my_type * d_axis = this->get_axis(0, len_p);
  p_ind = where(d_axis, len_p, p);
  d_axis = this->get_axis(1, len_a);
  alpha_ind = where(d_axis, len_a, angle_to_stored_angle(alpha));

  if(p_ind >= 0 && alpha_ind >= 0 && (size_t)p_ind < len_p && (size_t)alpha_ind < len_a){
    return this->get_element(p_ind, alpha_ind);
  }else{
    return 0.0;
  }

}

my_type diffusion_coeff::get_axis_element_ang(size_t ind){
/** \brief Get angle axs element
*
*Returns an ANGLE value, rather than raw axis entry
@param ind Index of element wanted
@return Angle on axis at ind
*/
  return stored_angle_to_angle(this->get_axis_element(1, ind));

}

