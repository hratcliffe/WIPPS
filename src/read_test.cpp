//
//  read_test.cpp
//  
//
//  Created by Heather Ratcliffe on 09/11/2015.
//
//

#include <stdio.h>
#include <fstream>
#include <iostream>
#include "sdf.h"
//SDF file libraries
#include <mpi.h>
#include "support.h"

using namespace std;
mpi_info_struc mpi_info;


int main(int argc, char *argv[]){
/**
*In theory, these classes and functions should be named well enough that the function here is largely clear. Remains to be seen, eh?
*
*/

  int ierr;
  
{
  int rank, n_procs;
  ierr = MPI_Init(&argc, &argv);
  //Note any other command line arg processing should account for this...
  ierr = MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  ierr = MPI_Comm_size(MPI_COMM_WORLD, &n_procs);

  mpi_info.rank = rank;
  mpi_info.n_procs = n_procs;
}


//  sdf_file_t *handle = sdf_open("/Users/heatherratcliffe/epoch_HR/epoch1d/Data/0000.sdf", MPI_COMM_WORLD, SDF_READ, 0);
  sdf_file_t *handle = sdf_open("./test0000.sdf", MPI_COMM_WORLD, SDF_READ, 0);

  sdf_block_t * block, * next;


  bool err=sdf_read_blocklist(handle);
  if(!err) cout<<"Yay! Blocks read"<<endl;
  else cout<<"Bum. Block read failed"<<endl;
/*
  block = sdf_find_block_by_id(handle, "ex_acc");

  if(block) cout<<"Success!!! Requested block found"<<endl;
  else{
    cout<<"Bleh. Requested block not found. Aborting"<<endl;
    return 1;
  }

  int n_dims = block->ndims;
  for(int i=0; i<n_dims; i++) cout<<block->dims[i]<<endl;
  cout<<block->next_block_location<<endl;


  block = sdf_find_block_by_id(handle, "grid_accum");

  if(block) cout<<"Success!!! Requested block found"<<endl;
  else{
    cout<<"Bleh. Requested block not found. Aborting"<<endl;
    return 1;
  }
  cout<<block->next_block_location<<endl;

  n_dims = block->ndims;
  for(int i=0; i<n_dims; i++) cout<<block->dims[i]<<endl;

  block = sdf_find_block_by_id(handle, "grid");

  if(block) cout<<"Success!!! Requested block found"<<endl;
  else{
    cout<<"Bleh. Requested block not found. Aborting"<<endl;
    return 1;
  }
  cout<<block->next_block_location<<endl;

  n_dims = block->ndims;
  for(int i=0; i<n_dims; i++) cout<<block->dims[i]<<endl;

  block = sdf_find_block_by_id(handle, "ex");

  if(block) cout<<"Success!!! Requested block found"<<endl;
  else{
    cout<<"Bleh. Requested block not found. Aborting"<<endl;
    return 1;
  }

  n_dims = block->ndims;
  for(int i=0; i<n_dims; i++) cout<<block->dims[i]<<endl;

  cout<<block->next_block_location<<endl;
*/

next = handle->current_block;
for(int i =0; i< handle->nblocks; i++){
  cout<<next->name<<" "<<next->id<<endl;
  cout<<next->blocktype<<endl;
    cout<<next->next_block_location<<endl;
    cout<<next->data_length<<endl;

  if(strcmp(next->name, "Electric Field/Ex")==0){
    cout<<"This one!"<<endl;
    cout<<"ndims "<<next->ndims<<endl;
    cout<<"dims are :"<<next->dims[0]<<" "<<next->dims[1]<<" "<<next->dims[2]<<endl;

  }

  //handle->current_block = next;
  next = next->next;
}

  
  sdf_close(handle);



}