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
#include "sdf.h"


class data_array;

/** \brief Reads SDF files into data_array

*
*Takes file prefixes, block id (see SDF documentation) and time and space ranges and a data_array to fill and does so.
* @author Heather Ratcliffe @date 02/10/2015 \ingroup cls
*/
class reader{

  int n_chars; /**<Number of characters in filename's dump number*/
  int ref_file_num;/**< Number of a file which exists to test dimensions, get reference B value etc*/

/********Filename and file manipulators ****/
  int check_file_num(int file_num);
  int get_filename_n_chars(int file_num);
  std::string get_full_name(int file_num);

/********Block manipulators ****/
  bool is_accum(std::string block_id);
  bool is_field_block(std::string block_id);

/********File read routines ****/
  int pre_read(data_array& data_in, int ref_time, bool accumulated, int flatten_on, size_t &n_dims, size_t * &source_sizes);
  void read_acc_time(data_array & data_in, sdf_file_t * handle, size_t start_pos, size_t n_times);
  void read_plain_time(data_array& data_in, sdf_file_t * handle, size_t pos);
public:

  std::string file_prefix;/**< Prefix of files before dump number*/
  size_t space_range[2];/**< Space range in x to extract*/
  int time_range[3];/**< Time range to extract*/
  char block_id[ID_SIZE];/**< Name of block to extract*/
/********Basic setup and allocation functions ****/
  explicit reader();
  explicit reader(std::string file_prefix_in,  const std::string block_id_in="", int ref_file_num_in=0);
  
/********Filename and file manipulators ****/
  void update_ref_filenum(int num);
  int get_file_size();
  
/********Block manipulators ****/
  bool change_block_id(std::string new_id);
  std::vector<std::pair<std::string, std::string> > list_blocks();
  bool current_block_is_accum();
  bool has_accum_data();

/********File read routines ****/
  bool read_dims(size_t &n_dims, std::vector<size_t> &dims);
  bool read_dims(size_t &n_dims, std::vector<size_t> &dims, std::string b_id);

  int read_data(data_array & my_data_in, size_t time_range[3], size_t space_range[2], int flatten_on = -1);
  
  bool read_distrib(data_array & my_data_in, std::string dist_id,int dump_number);

};

#endif
