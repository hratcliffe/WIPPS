//
//  my_array.h
//  
//
//  Created by Heather Ratcliffe on 21/09/2015.
//
//

#ifndef ____my_array__
#define ____my_array__

#include <stdlib.h> //for malloc eeejit
#include <vector>

class my_array{

protected:
  int n_dims;/** Number of dimensions*/
  int * dims;/** Array dimensions*/
  bool ragged; /**Flag for if array has different lengths on each row...*/
  int * row_lengths;/** Row lengths if ragged*/
  int * cumulative_row_lengths; /**Offset to start of given row to iterate through*/
  my_type *data;/** The data */
  bool defined; /**flag to check memory allocation sucess*/
  my_type * get_ptr(int nx, int ny);/** Shouldn't need this. Should use get/set_element*/

public:

  my_array();
  virtual void construct();
  my_array(int nx, int ny);
  my_array(int * row_lengths, int ny);
  virtual ~my_array();

  virtual bool is_good(){return !defined;}/** Check memory allocation etc worked*/

  virtual int get_index(int nx, int ny);
  int get_total_elements();
  /** These two account for all details of internal layout in memory */

  int get_dims();
  int get_dims(int dim);
  int get_length(int dim);

  my_type get_element(int nx, int ny);
  bool set_element(int nx, int ny, my_type val);

  bool populate_data(my_type * dat_in, int n_tot);
  bool populate_row(void * dat_in, int nx, int y_row);
  /** Populate data

  */
  virtual bool write_to_file(std::fstream &file);
  bool read_from_file(std::fstream &file);

  //we can use [] to wrap get elements and have the args pushed into vector which we then work with to be generic

  std::string array_self_test();

};


class data_array : public my_array{
/**Data array plus axes to hold field data. Includes ability to window in space/time */
protected:

  bool ax_defined;
  my_type *axes;/**1-d array in sections, so can be arbitary length and dims*/

public:

  char block_id[10]; /**the field name id from SDF file*/

  float time[2];/**time range over data are taken*/
  int space[2];/**space range ditto*/

  virtual void construct();
  data_array(int nx, int ny);
  data_array(int * row_lengths, int ny);
  virtual ~data_array();

  virtual bool is_good(){return defined && ax_defined;}
  my_type * get_axis(int dim, int & length);

  int get_axis_index(int nx, int ny);

  my_type get_axis_element(int dim, int pt);
  bool set_axis_element(int dim, int pt, my_type val);
  bool populate_axis(int dim, my_type * dat_in, int n_tot);
  float get_res(int i);
  int get_total_axis_elements();

  void make_linear_axis(int dim, float res, int offset=0);

  bool write_to_file(std::fstream &file);
  /** \todo Make mpi safe.. */
  bool read_from_file(std::fstream &file);

  bool fft_me(data_array * data_out);
  void copy_ids( data_array * src);

};



#endif /* defined(____my_array__) */
