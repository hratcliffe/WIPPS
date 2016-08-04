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
#include "my_array.h"

/** \brief Extended my_array class including axes
*
*Contains also data axes and various id's describing the data set from which data came. Also has ability to fft itself into a new instance. Move from file my_array by \author Heather Ratcliffe \date 3/08/2016
*/
class data_array : public my_array{

protected:

  bool ax_defined;/**< Flag showing whether axes are fully defined*/
  my_type *axes;/**< 1-d array in sections, so can be arbitary length and dims*/

  virtual void construct();
  void alloc_ax(const size_t els);
  float get_res(size_t i);
  size_t get_total_axis_elements();
  long get_axis_index(size_t dim, size_t pt)const;

  std::vector<size_t> get_bounds(std::vector<my_type> limits);
  void copy_ids( const data_array & src);


public:
  //friend class spectrum;/**< \todo Can we remove this plz?*/
  char block_id[ID_SIZE]; /**< The field name id from SDF file*/

  float time[2];/**< Time range over which data are taken*/
  size_t space[2];/**< Space range over which data are taken*/
  data_array();
  data_array(size_t nx, size_t ny=0, size_t nz=0, size_t nt=0);
  data_array(const data_array &src);
  data_array & operator=(const data_array& src);
  data_array(std::string filename, bool no_version_check = false);
  data_array(size_t n_dims, size_t * dims);
  virtual ~data_array();

  virtual bool is_good()const{return defined && ax_defined;}

  my_type get_axis_element(size_t dim, size_t pt)const;
  bool set_axis_element(size_t dim, size_t pt, my_type val);
  bool populate_axis(size_t dim, my_type * dat_in, size_t n_tot);
  my_type * get_axis(size_t dim, size_t & length);
  void make_linear_axis(size_t dim, float res, long offset=0);

  bool write_to_file(std::fstream &file, bool close_file=true);
  bool read_from_file(std::fstream &file, bool no_version_check=0);
  bool write_section_to_file(std::fstream &file, std::vector<my_type> limits, bool close_file=true);
  
  bool fft_me(data_array * data_out);
  bool populate_mirror_fastest(my_type * result_in, size_t total_els);
  bool check_ids(const data_array & src);
  bool resize(size_t dim, size_t sz);
  bool shift(size_t dim, long n_els, bool axis=1);
};



#endif
