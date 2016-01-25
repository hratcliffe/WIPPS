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
#include "support.h"
#include "controller.h"
#include "plasma.h"
#include "my_array.h"
#include "d_coeff.h"
#include "spectrum.h"

extern deck_constants my_const;
extern mpi_info_struc mpi_info;
controller::controller(std::string file_prefix){
/** \brief Setup
*
*\todo Plasma object should be setup from files or such
*/
  my_plas = new plasma(my_const.omega_ce * me/std::abs(q0), file_prefix);
  current_spect=0;
  current_d=0;

};

controller::~controller(){

  if(my_plas) delete my_plas;
  
  for(int i=0; i<my_spect.size(); i++){
    delete my_spect[i];
    my_spect[i] = nullptr;
  }
  for(int i=0; i<my_d.size(); i++){
    delete my_d[i];
    my_d[i] = nullptr;
  }

};

void controller::add_spectrum(int nx, int n_ang){
/** \brief Create and add spectrum
*
*rectangular spectrum, for when angular distribution depends on omega
*/
  spectrum * tmp_spect;
  tmp_spect = new spectrum(nx, n_ang);
  tmp_spect->my_controller = this;
  my_spect.push_back(tmp_spect);
  current_spect = my_spect.size()-1;

}

void controller::add_spectrum(int * row_lengths, int ny){
/** \brief Create and add spectrum
*
*Spectrum for when angular distribution does not depend on omega
*/
  spectrum * tmp_spect;
  tmp_spect = new spectrum(row_lengths, ny);
  tmp_spect->my_controller = this;
  my_spect.push_back(tmp_spect);
  current_spect = my_spect.size()-1;
  
}

void controller::add_d(int nx, int n_angs, int pos){
/** \brief Create and addd diffusion_coefficient and create it's basic axes etc
*
*
*/
  diffusion_coeff * tmp_d;
  tmp_d = new diffusion_coeff(nx, n_angs);
  tmp_d->my_controller = this;
  tmp_d->make_velocity_axis();
  tmp_d->make_pitch_axis();
  if(pos == -1 || pos == my_d.size()){
    my_d.push_back(tmp_d);
    current_d = my_d.size()-1;
  }
  else if(pos >=0 && pos< my_d.size()){
    //insert at position specified
    my_d.insert(my_d.begin() + pos, tmp_d);
    current_d = pos;
  }
}

spectrum * controller::get_current_spectrum(){

  if(!my_spect.empty()) return my_spect[current_spect];
  else return nullptr;

}

diffusion_coeff * controller::get_current_d(){
  
  if(!my_d.empty()) return my_d[current_d];
  else return nullptr;
}

void controller::bounce_average(){
  /** \todo finish this!!*/
  if(my_d.size() ==1) return;
  //No averaging to do!

  int dims[2];
  get_size(dims);
  add_d(dims[0], dims[1]);
  //Add new averaged D
  my_d[current_d]->tag = BOUNCE_AV;

    //Flatten down each d onto this one with all integrand included...
  my_type val;
  for(int j=0; j< dims[0]; j++){
    for(int k=0; k<dims[1]; k++){
      val = 0;
      for(int i=0; i<my_d.size() -1; i++){
        val += my_d[i]->get_element(j, k);
        
      }
// FAKENUMBERS      val *= integrand extras!!!
      my_d[current_d]->set_element(j, k, val);

    }
  }
  handle_d_mpi();
  //Now mpi reduce the results
  
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

  int current_d_keep = current_d;
  if(mpi_info.rank ==0){
    add_d(dims[0], dims[1], 0);
    my_d[current_d]->tag = GLOBAL;
  }
  current_d = current_d_keep;
  
  //Reduce onto this new root d so we preserve all the local versions too

  MPI_Reduce(my_d[0]->data, get_current_d()->data, total, MPI_MYTYPE, MPI_SUM, 0, MPI_COMM_WORLD);
  
}

void controller::get_size(int dims[2]){

  dims[0] = my_d[current_d]->get_dims(0);
  dims[1] = my_d[current_d]->get_dims(1);

}

bool controller::save_spectra(std::string pref){
/** \brief Save spectra to files (one per chunk)
*
* Writes each spectrum object to a file, identified by space range and time.
*/

  std::fstream file;
  std::string filename, tmp;
  tmp = my_spect[0]->block_id;
  for(int i=0; i<my_spect.size(); ++i){
    filename = pref+"spec_"+tmp +"_"+mk_str(my_spect[i]->time[0]) + "_"+mk_str(my_spect[i]->time[1])+"_"+mk_str(my_spect[i]->space[0])+"_"+mk_str(my_spect[i]->space[1])+".dat";
    file.open(filename.c_str(),std::ios::out|std::ios::binary);
    if(file.is_open()) my_spect[i]->write_to_file(file);
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
  tmp = my_d[0]->block_id;
  for(int i=0; i<my_d.size(); ++i){
    if(my_d[i]->tag == LOCAL) filename = pref+"D_"+tmp +"_"+mk_str(my_d[i]->time[0]) + "_"+mk_str(my_d[i]->time[1])+"_"+mk_str(my_d[i]->space[0])+"_"+mk_str(my_d[i]->space[1])+".dat";
    else if(my_d[i]->tag == BOUNCE_AV) continue;
    //Don't bother saving these...
    else if(my_d[i]->tag == GLOBAL) filename = pref+"D_"+tmp +"_"+mk_str(my_d[i]->time[0]) + "_"+mk_str(my_d[i]->time[1])+"_global.dat";
    
    file.open(filename.c_str(),std::ios::out|std::ios::binary);
    if(file.is_open()) my_d[i]->write_to_file(file);
    else return 1;
    file.close();
  
  }

  return 0;


}

