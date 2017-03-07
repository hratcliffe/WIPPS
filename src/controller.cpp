//
//  controller.cpp
//  
//
//  Created by Heather Ratcliffe on 19/11/2015.
//
//

#include <stdio.h>
#include <cmath>
#include <mpi.h>

#include "controller.h"
#include "my_array.h"
#include "d_coeff.h"
#include "spectrum.h"
#include "my_array.h"

/********Basic setup functions ****/
controller::controller(std::string file_prefix){
/** \brief Setup
*
* Create plasma object and initialise. A controller without a plasma is not meaningful. Plasma guarantees to be valid after construction, but may contain defaults if the specified file was not found. 
*/

  my_plas = plasma(file_prefix);
  current_pair=0;
  current_special_d=0;
}

controller::~controller(){
/** \brief Delete 
*
*Destroys all the D and spectrum objects
*/
  for(size_t i=0; i<spect_D_list.size(); i++){
    if(spect_D_list[i].first) delete spect_D_list[i].first;
    if(spect_D_list[i].second) delete spect_D_list[i].second;
    spect_D_list[i].first = nullptr;
    spect_D_list[i].second = nullptr;
  }
}

/********Plasma, spectrum, D setup functions ****/
bool controller::add_spectrum(std::string file){
/** \brief Add spectrum from file
*
*Add spectrum read from a file, created using e.g. data_array::write_to_file or write_data in IDL
*/

  //Read spectrum from file
  spectrum * tmp_spect = new spectrum(file);
  if(tmp_spect->is_good()){
    //Assign controller and add to list
    tmp_spect->my_controller = this;
    spect_D_list.push_back(std::make_pair(tmp_spect, nullptr));
    current_pair = spect_D_list.size()-1;
    return 0;
  }else{
    my_error_print("Spectrum construction failed", mpi_info.rank);
    return 1;
  }
}

bool controller::add_spectrum(int nx, int n_ang, bool separable){
/** \brief Create and add spectrum
*
*Create a spectrum of the specified size and adds to list. See spectrum::spectrum(int n_om, int n_ang, bool separable) for details.
*/
  //Create the spectrum
  spectrum * tmp_spect;
  tmp_spect = new spectrum(nx, n_ang, separable);
  if(tmp_spect->is_good()){
    //Assign controller and add to list
    tmp_spect->my_controller = this;
    spect_D_list.push_back(std::make_pair(tmp_spect, nullptr));
    current_pair = spect_D_list.size()-1;
    return 0;

  }else{
    my_error_print("Spectrum construction failed", mpi_info.rank);
    return 1;
  }
}

bool controller::add_d(int nx, int n_angs){
/** \brief Create and add diffusion_coefficient
*
*Creates a diffusion coefficient of size nx x n_angs, paired with the last spectrum that was added.
*/
  //Create empty diffusion coefficient
  diffusion_coeff * tmp_d;
  tmp_d = new diffusion_coeff(nx, n_angs);

  if(!tmp_d->is_good()) return 1;

  //Create axes etc and add to list
  tmp_d->my_controller = this;
  tmp_d->make_velocity_axis();
  tmp_d->make_pitch_axis();
  spect_D_list[current_pair].second = tmp_d;
  return 0;
}

void controller::add_d_special(int nx, int n_angs){
/** \brief Create and add special diffusion_coefficient
*
*Creates a diffusion coefficient of size nx x n_angs. Special coefficients do not have a matching spectrum. We use them mainly to hold results of bounce-averaging.
*/

  diffusion_coeff * tmp_d;
  tmp_d = new diffusion_coeff(nx, n_angs);
  tmp_d->my_controller = this;
  tmp_d->make_velocity_axis();
  tmp_d->make_pitch_axis();
  
  d_specials_list.push_back(tmp_d);
  current_special_d = d_specials_list.size()-1;
}

void controller::delete_current_spectrum(){
/** \brief Delete the current spectrum object
*
*Deletes the current (last added) spectrum and any attached D_coeff
*/
  if(spect_D_list[current_pair].first) delete spect_D_list[current_pair].first;
  if(spect_D_list[current_pair].second) delete spect_D_list[current_pair].second;
  spect_D_list.resize(spect_D_list.size() -1);

}

/********Plasma, spectrum, D getters ****/
void controller::get_D_size(int dims[2]){
/** \brief Get the dimensions of D
*
*Returns D dimensions in the dims array. Note that all D ought to be created the same size, but this is not guaranteed.
*/
  dims[0] = spect_D_list[current_pair].second->get_dims(0);
  dims[1] = spect_D_list[current_pair].second->get_dims(1);

}

spectrum * controller::get_current_spectrum(){
/** \brief Return current spectrum
*
* Return a pointer to the current (i.e. the latest added) spectrum, or nullptr if list is empty.
*/
  if(!spect_D_list.empty()) return spect_D_list[current_pair].first;
  else return nullptr;

}

spectrum * controller::get_spectrum_by_num(size_t indx){
/** \brief Return spectrum by indx
*
* Return a pointer to the spectrum indx ago, or nullptr if list is empty or indx is invalid.
*/
  if(indx >= spect_D_list.size()) return nullptr;
  if(!spect_D_list.empty()) return spect_D_list[current_pair - indx].first;
  else return nullptr;

}


diffusion_coeff * controller::get_current_d(){
/** \brief Return current D
*
* Return a pointer to the current (i.e. the latest added) d_coeff, or nullptr if list is empty.
*/
  if(!spect_D_list.empty()) return spect_D_list[current_pair].second;
  else return nullptr;
}
diffusion_coeff * controller::get_special_d(){
/** \brief Return current special D
*
* Return a pointer to the current (i.e. the latest added) special d_coeff, or nullptr if list is empty.
*/
  if(!d_specials_list.empty()) return d_specials_list[current_special_d];
  else return nullptr;
}

/********Bounce averaging specials ****/
void controller::bounce_average(bounce_av_data bounce_dat){
/** \brief Bounce average D
*
*Assumes the list contains D in order across space and performs bounce average to create special D. We use bounce_data to inform the field shape etc etc. The end result on each processor should be something which just has to be plain-summed by the mpi part. Calls controller::handle_d_mpi() and  \todo Finish integrand \todo Move these to dedicated code?
*/
  if(spect_D_list.size() ==0) return;
  //Empty list, nothing to do.

  int dims[2];
  get_D_size(dims);
  
  add_d_special(dims[0], dims[1]);
  //Add new averaged D
  d_specials_list[current_special_d]->tag = BOUNCE_AV;
  strcpy(d_specials_list[current_special_d]->block_id, spect_D_list[current_pair].second->block_id);
  //Tag as bounce av'd and copy block id from block used

  //Flatten down each d onto this one with all integrand included...
  my_type val, current_lat = 0.0, d_lat = 0.0, D_val, lat_factor;
  size_t n_blocks = spect_D_list.size();
  my_type sin_theta, cos_theta_sq;
  //Latitude increment block to block
  d_lat = bounce_dat.max_latitude/(my_type) n_blocks;
  //Initial block-centred latitude
  for(int p_i = 0; p_i < dims[0]; p_i++){
    for(int ang_i = 0; ang_i < dims[1]; ang_i++){
      val = 0.0;
      current_lat = d_lat/2.0;
      for(size_t block_i = 0; block_i < n_blocks; block_i++){
        sin_theta = std::sin(pi* (90.0-current_lat)/180.0);
        cos_theta_sq = 1.0 - sin_theta*sin_theta;
        D_val = spect_D_list[block_i].second->get_element(p_i, ang_i); //* integrand!
        switch (bounce_dat.type) {
          case plain:
            lat_factor = bounce_dat.L_shell * R_E * std::sqrt(1.0+3.0*cos_theta_sq) * sin_theta * (pi*d_lat/180.0);//This is ds as in Schulz/Lanzerotti Eq 1.18.
            break;
          case alpha_alpha:
          /** \todo Fill other cases*/
            break;
          case alpha_p:
            break;
          case p_p:
            break;
        }
        val += D_val * lat_factor;
        //Sum up the blocks done on this processor
        current_lat += d_lat;
      }
      d_specials_list[current_special_d]->set_element(p_i, ang_i, val);
    }
  }
  handle_d_mpi();
  //Now mpi reduce the results to get a global average
  
}

void controller::handle_d_mpi(){
/** \brief MPI reduce diffusion coeffs
*
* Creates an MPI Summed d on the root node
*/

//Make a reduced d on root, but leaving the current_special_d counter unchanged.
  int dims[2];
  get_D_size(dims);
  int total = dims[0]*dims[1];

  int prev_current_d = current_special_d;
  if(mpi_info.rank ==0){
  //Add a new global value on root only
    add_d_special(dims[0], dims[1]);
    //Copy roots bounce average into this
    *(d_specials_list[current_special_d]) = *(d_specials_list[prev_current_d]);
    d_specials_list[current_special_d]->tag = GLOBAL;
  }
  //reset current_special_d (affects root only)
  int root_global_d = current_special_d;
  current_special_d = prev_current_d;

  //Reduce onto this new root d so we preserve all the local versions too
  MPI_Reduce(get_current_d()->data, d_specials_list[root_global_d]->data , total, MPI_MYTYPE, MPI_SUM, 0, MPI_COMM_WORLD);
  
}

/********File IO functions ****/
bool controller::save_spectra(std::string pref){
/** \brief Save spectra to files (one per chunk)
*
* Writes each spectrum object to a file, identified by space range and time. \todo Add logging?
*/

  std::fstream file;
  std::string filename, tmp;
  for(size_t i=0; i<spect_D_list.size(); ++i){
    tmp = spect_D_list[i].first->block_id;
    filename = pref+"spec_"+tmp +"_"+mk_str(spect_D_list[i].first->time[0]) + "_"+mk_str(spect_D_list[i].first->time[1])+"_"+mk_str(spect_D_list[i].first->space[0])+"_"+mk_str(spect_D_list[i].first->space[1])+".dat";
    file.open(filename.c_str(),std::ios::out|std::ios::binary);
    if(file.is_open()) spect_D_list[i].first->write_to_file(file);
    else return 1;
    file.close();
  }
  return 0;
}

bool controller::save_D(std::string pref){
/** \brief Save D's to files (one per chunk)
*
* Writes each spectrum object to a file, identified by space range and time. Note root will also write a bounce averaged file
*/

  std::fstream file;
  std::string filename, tmp;

  //Do local list
  for(size_t i=0; i<spect_D_list.size(); ++i){
    tmp = spect_D_list[i].second->block_id;
    //They might have different blocks
    filename = pref+"D_"+tmp +"_"+mk_str(spect_D_list[i].second->time[0]) + "_"+mk_str(spect_D_list[i].second->time[1])+"_"+mk_str(spect_D_list[i].second->space[0])+"_"+mk_str(spect_D_list[i].second->space[1])+".dat";
    //Don't bother saving these...
    
    file.open(filename.c_str(),std::ios::out|std::ios::binary|std::ios::in|std::ios::trunc);
    if(file.is_open()) spect_D_list[i].second->write_to_file(file);
    else return 1;
    file.close();
  
  }

  //Save global D from root only
  int d_global = d_specials_list.size()-1;
  if(mpi_info.rank == 0 ){
    filename = pref+"D_"+tmp +"_"+mk_str(d_specials_list[d_global]->time[0]) + "_"+mk_str(d_specials_list[d_global]->time[1])+"_global.dat";
    file.open(filename.c_str(),std::ios::out|std::ios::binary|std::ios::in|std::ios::trunc);
    if(file.is_open()) d_specials_list[d_global]->write_to_file(file);
    else return 1;
    file.close();
  }
  return 0;
}

