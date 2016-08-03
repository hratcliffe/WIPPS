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

/** \brief A basic array class
*
*Contains dimension information and data. Can be rectangular of any n_dims or ragged of 2 (rows of different lengths). Get_index and get_total_elements account for all details of internal layout in memory. For 1-4 dims individual getter/setter functions are given. For larger arrays one must construc the array of indexes. NOTE the backing memory is old style with Fortran style internal ordering (for ease of SDF interfacing). But contigous memory and pointer arithmetic give major speed advantage and we very rarely change size on the fly. However nothing outside this class should need to do anything except access by index and populate by element, slice or entire. Internal ordering is Fortran style (for ease of SDF interfacing). \todo Break out data array to file \todo Check a vector works
 \author Heather Ratcliffe \date 21/09/2015
*/

class my_array{

protected:
  size_t n_dims;/**< Number of dimensions*/
  size_t * dims;/**< Array dimensions*/
  my_type *data;/**< The data */
  bool defined; /**< Flag to check memory allocation sucess*/

  virtual void construct();
  virtual void alloc_all(const size_t n_dims, const size_t * const dims);
  
  virtual std::vector<size_t> get_index_from_offset(size_t offset);
  virtual long get_index(size_t n_dims, size_t * dim);
  long get_index(size_t nx);
  long get_index(size_t nx, size_t ny);
  long get_index(size_t nx, size_t ny, size_t nz);
  long get_index(size_t nx, size_t ny, size_t nz, size_t nt);

  size_t get_total_elements();

public:

  my_array();
  my_array(size_t nx, size_t ny=0, size_t nz=0, size_t nt=0);
//  my_array(size_t * row_lengths, size_t ny);
  my_array(size_t n_dims, size_t * dims);
  virtual ~my_array();
  my_array(const my_array &src);
  my_array & operator=(const my_array& src);
  virtual bool is_good()const {return !defined;}/**< Check memory allocation etc worked*/

  size_t get_dims();
  size_t get_dims(size_t dim);
  size_t get_length(size_t dim);
  my_type get_element(size_t nx);
  my_type get_element(size_t nx, size_t ny);
  my_type get_element(size_t nx, size_t ny, size_t nz);
  my_type get_element(size_t nx, size_t ny, size_t nz, size_t nt);
  my_type get_element(size_t n_dims, size_t * dim);
  bool set_element(size_t nx, my_type val);
  bool set_element(size_t nx, size_t ny, my_type val);
  bool set_element(size_t nx, size_t ny, size_t nz, my_type val);
  bool set_element(size_t nx, size_t ny, size_t nz, size_t nt, my_type val);
  bool set_element(size_t n_dims, size_t * dim, my_type val);

  bool populate_data(my_type * dat_in, size_t n_tot);
  bool populate_slice(my_type * dat_in, size_t n_dims, size_t * offsets);
  bool populate_complex_slice(my_type * dat_in, size_t n_dims, size_t * offsets, size_t* sizes);

  virtual bool write_to_file(std::fstream &file);
  virtual bool read_from_file(std::fstream &file, bool no_version_check=0);
  virtual bool write_section_to_file(std::fstream &file, std::vector<size_t> bounds);
  virtual bool resize(size_t dim, size_t sz);
  virtual bool shift(size_t dim, long n_els);
  
  my_type minval(size_t offset=0);
  my_type maxval(size_t offset=0);
  my_type minval(std::vector<size_t> &ind, size_t offset=0);
  my_type maxval(std::vector<size_t> &ind, size_t offset=0);

};

/** \brief Extended my_array class including axes
*
*Contains also data axes and various id's describing the data set from which data came. Also has ability to fft itself into a new instance. \author Heather Ratcliffe \date 21/09/2015
*/
class data_array : public my_array{

protected:

  bool ax_defined;/**< Flag showing whether axes are fully defined*/
  my_type *axes;/**< 1-d array in sections, so can be arbitary length and dims*/

  virtual void construct();
  void alloc_ax(const size_t els);
  float get_res(size_t i);
  size_t get_total_axis_elements();
  long get_axis_index(size_t dim, size_t pt);

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

  virtual bool is_good(){return defined && ax_defined;}

  my_type get_axis_element(size_t dim, size_t pt);
  bool set_axis_element(size_t dim, size_t pt, my_type val);
  bool populate_axis(size_t dim, my_type * dat_in, size_t n_tot);
  my_type * get_axis(size_t dim, size_t & length);
  void make_linear_axis(size_t dim, float res, long offset=0);

  virtual bool write_to_file(std::fstream &file);
  virtual bool read_from_file(std::fstream &file, bool no_version_check=0);
  virtual bool write_section_to_file(std::fstream &file, std::vector<my_type> limits);
  
  bool fft_me(data_array * data_out);
  bool populate_mirror_fastest(my_type * result_in, size_t total_els);
  bool check_ids(const data_array & src);
  virtual bool resize(size_t dim, size_t sz);
  virtual bool shift(size_t dim, long n_els, bool axis=1);

};



#endif /* defined(____my_array__) */
