//
//  reader.h
//  
//
//  Created by Heather Ratcliffe on 02/10/2015.
//
//

#ifndef _reader_h
#define _reader_h

#include <iostream>
#include <vector>


class data_array;

/** \brief Reads SDF files into data_array

*
*Takes file prefixes, block id (see SDF documentation) and time and space ranges and a data_array pointer to fill and does so. The data array is not stored here.
*\todo maybe we write a verify sdf which checks our files have the correct dimensionalities etc etc and contain needed blocks...
* @author Heather Ratcliffe @date 02/10/2015
*/
class reader{

  int n_z; /*<*Number of characters in filename's dump number*/
  int ref_file_num;/**< Number of any file which exists to test dimensions etc*/
public:

  std::string file_prefix;/**< Prefix of files before dump number*/

  int space_range[2];/**< Space range in x to extract*/
  int time_range[2];/**< Time range to extract*/

  char block_id[10];/**< Name of block to extract*/

  reader(std::string file_prefix_in,  char * block_id_in, int first=0);
  ~reader(){;}

  bool read_dims(int &n_dims, std::vector<int> &dims);

  bool read_data(data_array * my_data_in, int time_range[2], int space_range[2]);
  int get_file_size();
  std::string get_full_name(int num);

};

#endif
