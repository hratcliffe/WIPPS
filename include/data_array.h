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
*Extends the my_array class to add axes and some data tags describing the data set derived from.

*Move from file my_array by \author Heather Ratcliffe \date 3/08/2016 \ingroup cls
*/
class data_array : public my_array{

protected:

  my_type *axes;/**< Axes data*/
  
/********Basic setup and allocation functions ****/
  virtual void construct();
  void alloc_ax(const size_t els);

/********Indexers, getters and setters ****/
  size_t get_total_axis_elements()const;
  long get_axis_index(size_t dim, size_t pt)const;
  std::vector<size_t> get_bounds(std::vector<my_type> limits);

public:
  char block_id[ID_SIZE]; /**< ID describing data*/
  my_type time[2];/**< Time range over which data are taken*/
  size_t space[2];/**< Space range over which data are taken*/
  my_type B_ref;/**< Reference average B field by location*/

/********Basic setup and allocation functions ****/
  explicit data_array();
  explicit data_array(size_t nx, size_t ny=0, size_t nz=0, size_t nt=0);
  explicit data_array(size_t n_dims, size_t * dims);
  explicit data_array(std::string filename);
  virtual ~data_array();
  virtual bool is_good()const{return my_array::is_good() && axes;}/**< Whether array is useable @return Boolean true for good state, false for bad*/

/********Technical stuff making my_array a proper "object" ****/
  data_array(const data_array &src);
  data_array(data_array && src);
  data_array & operator=(const data_array& src);
  bool operator==(const data_array &rhs)const;
  bool operator!=(const data_array &rhs)const{return !(*this == rhs);}/**< See data_array::operator==()*/
  data_array(const my_array &src);
  data_array & operator=(const my_array & src);
  bool operator==(const my_array &rhs)const{return my_array::operator==(rhs);}/**< Equality (size and data only) with a my_array @param rhs Array to compare to*/
  bool operator!=(const my_array &rhs)const{return !(*this == rhs);}/**< See data_array::operator==()*/

/********Helpers for working with data_array ****/
  my_type * disown_axes();
  void clone_empty(const data_array &src);
  void copy_ids( const data_array & src);
  bool check_ids(const data_array & src)const;

/********Indexers, getters and setters ****/
  my_type get_axis_element(size_t dim, size_t pt)const;
  bool set_axis_element(size_t dim, size_t pt, my_type val);
  const my_type * get_axis(size_t dim, size_t & length);
  float get_res(size_t i)const;
  long get_axis_index_from_value(size_t dim, my_type value)const;

/********Data/axis fillers, file IO ****/
  bool populate_axis(size_t dim, my_type * dat_in, size_t n_tot);
  void make_linear_axis(size_t dim, float res, long offset=0);

  bool write_to_file(std::fstream &file, bool close_file=true);
  bool read_from_file(std::fstream &file);
  bool write_section_to_file(std::fstream &file, std::vector<my_type> limits, bool close_file=true);
  bool write_raw_section_to_file(std::fstream &file, std::vector<size_t> index_limits, bool close_file=true);
  bool write_closer(std::fstream &file);
  
/********Helpful functions working on entire array as a thing ****/
  bool resize(size_t dim, size_t sz, bool verbose=0);
  bool shift(size_t dim, long n_els, bool axis=1);
  data_array total(size_t dim);
  data_array total(size_t dim, my_type min, my_type max);
  data_array total(size_t dim, size_t min_ind, size_t max_ind);
  data_array average(size_t dim, my_type min, my_type max);
  data_array average(size_t dim);
  data_array average(size_t dim, size_t min, size_t max);
};

/********FFT functions and helpers ****/
bool fft_array(const data_array &data_in, data_array &data_out);
bool populate_mirror_fastest(data_array &data_out, my_type * result_in, size_t total_els);

/********Other non-member helpers ****/
my_type avval(const data_array & array_in);

#endif
