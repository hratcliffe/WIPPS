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
*Contains dimension information and data. Can be rectangular of any n_dims or ragged of 2 (rows of different lengths). Get_index and get_total_elements account for all details of internal layout in memory. For 1-4 dims individual getter/setter functions are given. For larger arrays one must construc the array of indexes. NOTE the backing memory is old style with Fortran style internal ordering (for ease of SDF interfacing). But contigous memory and pointer arithmetic give major speed advantage and we very rarely change size on the fly. However nothing outside this class should need to do anything except access by index and populate by element, slice or entire. Internal ordering is Fortran style (for ease of SDF interfacing).
 \author Heather Ratcliffe \date 21/09/2015
*/

class my_array{

protected:
  size_t n_dims;/**< Number of dimensions*/
  size_t * dims;/**< Array dimensions*/
  my_type *data;/**< The data */

  virtual void construct();
  virtual void alloc_all(const size_t n_dims, const size_t * const dims);
  
  virtual std::vector<size_t> get_index_from_offset(size_t offset)const;
  virtual long get_index(size_t n_dims, size_t * dim)const;
  long get_index(size_t nx)const;
  long get_index(size_t nx, size_t ny)const;
  long get_index(size_t nx, size_t ny, size_t nz)const;
  long get_index(size_t nx, size_t ny, size_t nz, size_t nt)const;
  my_type get_element_from_index(size_t ind)const;

public:
  size_t get_total_elements()const;

  my_array();
  my_array(size_t nx, size_t ny=0, size_t nz=0, size_t nt=0);
  my_array(size_t n_dims, size_t * dims);
  virtual ~my_array();
  my_array(const my_array &src);
  my_array(my_array && src);
  my_array & operator=(const my_array& src);

  template<typename T> void divide( T val){
    for(size_t i=0; i< this->get_total_elements(); i++) *(this->data+i) /= val;
  }//Works for any type T where my_type/T is defined

  virtual bool is_good()const {return (data && dims);}/**< Check memory allocation etc worked*/
  void clone_empty(const my_array &src);
  size_t get_dims()const;
  size_t get_dims(size_t dim)const;
  my_type get_element(size_t nx)const;
  my_type get_element(size_t nx, size_t ny)const;
  my_type get_element(size_t nx, size_t ny, size_t nz)const;
  my_type get_element(size_t nx, size_t ny, size_t nz, size_t nt)const;
  my_type get_element(size_t n_dims, size_t * dim)const;
  bool set_element(size_t nx, my_type val);
  bool set_element(size_t nx, size_t ny, my_type val);
  bool set_element(size_t nx, size_t ny, size_t nz, my_type val);
  bool set_element(size_t nx, size_t ny, size_t nz, size_t nt, my_type val);
  bool set_element(size_t n_dims, size_t * dim, my_type val);

  bool populate_data(my_type * dat_in, size_t n_tot, bool convert = 0);
  bool populate_slice(my_type * dat_in, size_t n_dims, size_t * offsets);
  bool populate_complex_slice(my_type * dat_in, size_t n_dims, size_t * offsets, size_t* sizes);

  bool write_to_file(std::fstream &file);
  bool read_from_file(std::fstream &file, bool no_version_check=0);
  std::vector<size_t> read_dims_from_file(std::fstream &file, bool no_version_check=0);
  
  my_type * disown_data();
  
  bool write_section_to_file(std::fstream &file, std::vector<size_t> bounds);
  bool resize(size_t dim, size_t sz, bool verbose=0);
  bool shift(size_t dim, long n_els);
  
  my_type minval(size_t offset=0);
  my_type maxval(size_t offset=0);
  my_type minval(std::vector<size_t> &ind, size_t offset=0);
  my_type maxval(std::vector<size_t> &ind, size_t offset=0);
  my_type avval();
  my_type partial_maxval(std::vector<std::pair<size_t, size_t> > ranges,std::vector<size_t> &ind);
  void smooth_1d(int n_pts);
};

#endif /* defined(____my_array__) */
