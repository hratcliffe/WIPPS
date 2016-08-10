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
*Takes file prefixes, block id (see SDF documentation) and time and space ranges and a data_array to fill and does so. NOTE dumped arrays etc can be read directly, see e.g. data_array(std::string filename);
* @author Heather Ratcliffe @date 02/10/2015
*/
class reader{

  int n_z; /*<*Number of characters in filename's dump number*/
  int ref_file_num;/**< Number of any file which exists to test dimensions etc*/
  std::string get_full_name(int num);
  bool is_accum(std::string block_id);

public:

  std::string file_prefix;/**< Prefix of files before dump number*/

  int space_range[2];/**< Space range in x to extract*/
  int time_range[3];/**< Time range to extract*/

  char block_id[ID_SIZE];/**< Name of block to extract*/
  reader();
  reader(std::string file_prefix_in,  char * block_id_in, int first=0);
  ~reader(){;}

  bool read_dims(size_t &n_dims, std::vector<size_t> &dims);

  int read_data(data_array & my_data_in, int time_range[3], int space_range[2]);
  bool current_block_is_accum();
  int get_file_size();
  
};

#endif
