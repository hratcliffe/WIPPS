//
//  data_array.h
//  
//
//  Created by Heather Ratcliffe on 03/08/2016.
//
//

#ifndef _data_array_h
#define _data_array_h

#include <stdlib.h>
#include <vector>
#ifndef NO_FFT
#include <fftw3.h> //Include here so propagates to anywhere used
#endif
#include "my_array.h"

/** \brief Extended my_array class including axes
*
*Contains also data axes and various id's describing the data set from which data came. Also has ability to fft itself into a new instance. Move from file my_array by \author Heather Ratcliffe \date 3/08/2016
*/
class data_array : public my_array{

protected:

  my_type *axes;/**< Axes data*/
  
/********Basic setup and allocation functions ****/
  virtual void construct();
  void alloc_ax(const size_t els);

/********Indexers, getters and setters ****/
  size_t get_total_axis_elements();
  long get_axis_index(size_t dim, size_t pt)const;
  std::vector<size_t> get_bounds(std::vector<my_type> limits);

public:
  char block_id[ID_SIZE]; /**< ID describing data*/
  my_type time[2];/**< Time range over which data are taken*/
  size_t space[2];/**< Space range over which data are taken*/
  my_type B_ref;/**< Reference average B field by location*/

/********Basic setup and allocation functions ****/
  data_array();
  data_array(size_t nx, size_t ny=0, size_t nz=0, size_t nt=0);
  data_array(size_t n_dims, size_t * dims);
#ifdef DEFAULT_NOVERS
  data_array(std::string filename, bool no_version_check = true);
#else
  data_array(std::string filename, bool no_version_check = false);
#endif
  virtual ~data_array();
  virtual bool is_good()const{return my_array::is_good() && axes;}

/********Technical stuff making my_array a proper "object" ****/
  data_array(const data_array &src);
  data_array(data_array && src);
  data_array & operator=(const data_array& src);

/********Helpers for working with data_array ****/
  my_type * disown_axes();
  void clone_empty(const data_array &src);
  void copy_ids( const data_array & src);
  bool check_ids(const data_array & src);

/********Indexers, getters and setters ****/
  my_type get_axis_element(size_t dim, size_t pt)const;
  bool set_axis_element(size_t dim, size_t pt, my_type val);
  my_type * get_axis(size_t dim, size_t & length);
  float get_res(size_t i)const;

/********Data/axis fillers, file IO ****/
  bool populate_axis(size_t dim, my_type * dat_in, size_t n_tot);
  void make_linear_axis(size_t dim, float res, long offset=0);

  bool write_to_file(std::fstream &file, bool close_file=true);
#ifdef DEFAULT_NOVERS
  bool read_from_file(std::fstream &file, bool no_version_check=true);
#else
  bool read_from_file(std::fstream &file, bool no_version_check=false);
#endif
  bool write_section_to_file(std::fstream &file, std::vector<my_type> limits, bool close_file=true);
  bool write_raw_section_to_file(std::fstream &file, std::vector<size_t> limits, bool close_file=true);
  bool write_closer(std::fstream &file);
  
/********Helpful functions working on entire array as a thing ****/
  bool resize(size_t dim, size_t sz, bool verbose=0);
  bool shift(size_t dim, long n_els, bool axis=1);
  data_array total(size_t dim);
  data_array total(size_t dim, my_type min, my_type max);
  data_array total(size_t dim, size_t min, size_t max);
  data_array average(size_t dim, my_type min, my_type max);
  data_array average(size_t dim);
  data_array average(size_t dim, size_t min, size_t max);
  bool subtract(const data_array& rhs);

};

/********FFT functions and helpers ****/
bool fft_array(const data_array &data_in, data_array &data_out);
bool populate_mirror_fastest(data_array &data_out, my_type * result_in, size_t total_els);



#endif
