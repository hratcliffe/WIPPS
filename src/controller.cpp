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


extern deck_constants my_const;
extern const mpi_info_struc mpi_info;

controller::controller(std::string file_prefix){
/** \brief Setup
*
* Create plasma object and initialise
*/

  my_plas = plasma(file_prefix);
  current_spect=0;
  current_d=0;
};

controller::~controller(){
/** \brief Delete 
*
*Destroys all the D and spectrum objects
*/
  for(size_t i=0; i<my_list.size(); i++){
    if(my_list[i].first) delete my_list[i].first;
    if(my_list[i].second) delete my_list[i].second;
    my_list[i].first = nullptr;
    my_list[i].second = nullptr;
  }
};

void controller::add_spectrum(std::string file){
/** \brief Add spectrum from dump
*
*Add spectrum read from file dump
*/
  spectrum * tmp_spect = new spectrum(file);
  if(tmp_spect->is_good()){
    tmp_spect->my_controller = this;
    tmp_spect->init();
    my_list.push_back(std::make_pair(tmp_spect, nullptr));
    current_spect = my_list.size()-1;
  }else{
    my_print("Spectrum construction failed", mpi_info.rank);
  }
}

void controller::add_spectrum(int nx, int n_ang, bool separable){
/** \brief Create and add spectrum
*
*If separable is true it is assumed the angle distrib does not depend on omega, else it does
*/
  spectrum * tmp_spect;
  tmp_spect = new spectrum(nx, n_ang, separable);
  if(tmp_spect->is_good()){
    tmp_spect->my_controller = this;
    my_list.push_back(std::make_pair(tmp_spect, nullptr));
    current_spect = my_list.size()-1;

  }else{
    my_print("Spectrum construction failed", mpi_info.rank);
  }
}

void controller::add_d(int nx, int n_angs){
/** \brief Create and add diffusion_coefficient
*
*Creates a diffusion coefficient of size nx x n_angs, paired with the last spectrum that was added.
*/
  diffusion_coeff * tmp_d;
  tmp_d = new diffusion_coeff(nx, n_angs);
  tmp_d->my_controller = this;
  tmp_d->make_velocity_axis();
  tmp_d->make_pitch_axis();
  
  my_list[current_spect].second = tmp_d;
}
void controller::add_d_special(int nx, int n_angs){

  diffusion_coeff * tmp_d;
  tmp_d = new diffusion_coeff(nx, n_angs);
  tmp_d->my_controller = this;
  tmp_d->make_velocity_axis();
  tmp_d->make_pitch_axis();
  
  d_specials.push_back(tmp_d);
  current_d = d_specials.size()-1;
}
spectrum * controller::get_current_spectrum(){
/** \brief Return current spectrum
*
*Returns pointer because spectra vector may be empty
*/
  if(!my_list.empty()) return my_list[current_spect].first;
  else return nullptr;

}

diffusion_coeff * controller::get_current_d(){
/** \brief Return current D
*
*Returns pointer because D list may be empty
*/
  if(!my_list.empty()) return my_list[current_spect].second;
  else return nullptr;
}

void controller::bounce_average(){
/** \brief Bounce average D
*
*Assumes the list contains D in order across space and performs bounce average to create special D. \todo Finish integrand \todo test case?
*/
  if(my_list.size() ==0) return;
  //Empty list, nothing to do.

  int dims[2];
  get_size(dims);
  
  add_d_special(dims[0], dims[1]);
  //Add new averaged D
  d_specials[current_d]->tag = BOUNCE_AV;
  strcpy(d_specials[current_d]->block_id, my_list[current_spect].second->block_id);
  //Tag as bounce av'd and copy block id from block used

  //Flatten down each d onto this one with all integrand included...
  my_type val;
  for(int j=0; j< dims[0]; j++){
    for(int k=0; k<dims[1]; k++){
      val = 0;
      for(size_t i=0; i<my_list.size() -1; i++){
        val += my_list[i].second->get_element(j, k); //* integrand!
        //Sum up the blocks done on this processor
      }
      d_specials[current_d]->set_element(j, k, val);

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

//Make a reduced d on root, but leaving the current_d counter unchanged.
  int dims[2];
  get_size(dims);
  int total = dims[0]*dims[1];

  int prev_current_d = current_d;
  if(mpi_info.rank ==0){
  //Add a new global value on root only
    add_d_special(dims[0], dims[1]);
    //Copy roots bounce average into this
    *(d_specials[current_d]) = *(d_specials[prev_current_d]);
    d_specials[current_d]->tag = GLOBAL;
  }
  //reset current_d (affects root only)
  int root_global_d = current_d;
  current_d = prev_current_d;

  //Reduce onto this new root d so we preserve all the local versions too
  MPI_Reduce(get_current_d()->data, d_specials[root_global_d]->data , total, MPI_MYTYPE, MPI_SUM, 0, MPI_COMM_WORLD);
  
}

void controller::get_size(int dims[2]){

  dims[0] = my_list[current_spect].second->get_dims(0);
  dims[1] = my_list[current_spect].second->get_dims(1);

}

bool controller::save_spectra(std::string pref){
/** \brief Save spectra to files (one per chunk)
*
* Writes each spectrum object to a file, identified by space range and time. \todo Add logging?
*/

  std::fstream file;
  std::string filename, tmp;
  for(size_t i=0; i<my_list.size(); ++i){
    tmp = my_list[i].first->block_id;
    filename = pref+"spec_"+tmp +"_"+mk_str(my_list[i].first->time[0]) + "_"+mk_str(my_list[i].first->time[1])+"_"+mk_str(my_list[i].first->space[0])+"_"+mk_str(my_list[i].first->space[1])+".dat";
    file.open(filename.c_str(),std::ios::out|std::ios::binary);
    if(file.is_open()) my_list[i].first->write_to_file(file);
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
  for(size_t i=0; i<my_list.size(); ++i){
    tmp = my_list[i].second->block_id;
    //They might have different blocks
    filename = pref+"D_"+tmp +"_"+mk_str(my_list[i].second->time[0]) + "_"+mk_str(my_list[i].second->time[1])+"_"+mk_str(my_list[i].second->space[0])+"_"+mk_str(my_list[i].second->space[1])+".dat";
    //Don't bother saving these...
    
    file.open(filename.c_str(),std::ios::out|std::ios::binary|std::ios::in|std::ios::trunc);
    if(file.is_open()) my_list[i].second->write_to_file(file);
    else return 1;
    file.close();
  
  }

  //Save global D from root only
  int d_global = d_specials.size()-1;
  if(mpi_info.rank == 0 ){
    filename = pref+"D_"+tmp +"_"+mk_str(d_specials[d_global]->time[0]) + "_"+mk_str(d_specials[d_global]->time[1])+"_global.dat";
    file.open(filename.c_str(),std::ios::out|std::ios::binary|std::ios::in|std::ios::trunc);
    if(file.is_open()) d_specials[d_global]->write_to_file(file);
    else return 1;
    file.close();
  }
  return 0;
}

