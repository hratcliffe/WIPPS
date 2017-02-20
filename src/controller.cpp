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

diffusion_coeff * controller::get_current_d(){
/** \brief Return current D
*
* Return a pointer to the current (i.e. the latest added) d_coeff, or nullptr if list is empty.
*/
  if(!spect_D_list.empty()) return spect_D_list[current_pair].second;
  else return nullptr;
}

/********Bounce averaging specials ****/
void controller::bounce_average(){
/** \brief Bounce average D
*
*Assumes the list contains D in order across space and performs bounce average to create special D. \todo Finish integrand \todo test case? \todo Move these to dedicated code?
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
  my_type val;
  for(int j=0; j< dims[0]; j++){
    for(int k=0; k<dims[1]; k++){
      val = 0;
      for(size_t i=0; i<spect_D_list.size() -1; i++){
        val += spect_D_list[i].second->get_element(j, k); //* integrand!
        //Sum up the blocks done on this processor
      }
      d_specials_list[current_special_d]->set_element(j, k, val);

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

