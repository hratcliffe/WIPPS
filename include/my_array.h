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
#include <functional>

//Set debug flag if we've set DEBUG_ALL
#ifdef DEBUG_ALL
#define DEBUG_DIMS
#define DEBUG_RANGE
#endif

class data_array;

/** \brief A basic array class
*
*Contains dimension information and data. Can be rectangular of any n_dims or ragged of 2 (rows of different lengths). Get_index and get_total_elements account for all details of internal layout in memory. For 1-4 dims individual getter/setter functions are given. For larger arrays one must construc the array of indexes. NOTE the backing memory is old style with Fortran style internal ordering (for ease of SDF interfacing). But contigous memory and pointer arithmetic give major speed advantage and we very rarely change size on the fly. However nothing outside this class should need to do anything except access by index and populate by element, slice or entire. Internal ordering is Fortran style (for ease of SDF interfacing).
 \author Heather Ratcliffe \date 21/09/2015 \ingroup cls
*/
class my_array{

protected:
  size_t n_dims;/**< Number of dimensions*/
  size_t * dims;/**< Array dimensions*/
  my_type *data;/**< The data */

/********Basic setup and allocation functions ****/
  virtual void construct();
  virtual void alloc_all(const size_t n_dims, const size_t * const dims);
  
/********Indexers, getters and setters ****/
  virtual std::vector<size_t> get_indices_from_offset(size_t offset)const;
  virtual long get_index(size_t n_dims, size_t * dim)const;
  long get_index(size_t nx)const;
  long get_index(size_t nx, size_t ny)const;
  long get_index(size_t nx, size_t ny, size_t nz)const;
  long get_index(size_t nx, size_t ny, size_t nz, size_t nt)const;
  my_type get_element_from_index(size_t ind)const;

public:
/********Basic setup and allocation functions ****/
  //Explicit constructors, can only be called directly, not used in conversions etc
  explicit my_array();
  explicit my_array(size_t nx, size_t ny=0, size_t nz=0, size_t nt=0);
  explicit my_array(size_t n_dims, size_t * dims);
  virtual ~my_array();
  virtual bool is_good()const {return (data && dims);}/**< Check memory allocation etc worked @return True if good, false else*/

/********Technical stuff making my_array a proper "object" ****/
  my_array(const my_array &src);
  my_array(my_array && src);
  my_array & operator=(const my_array& src);
  bool operator==(const my_array &rhs)const;
  bool operator!=(const my_array &rhs)const{return !(*this == rhs);}/**< See my_array::operator==()*/

/********Helpers for working with my_array ****/
  my_type * disown_data();
  void clone_empty(const my_array &src);
  bool copy_data(my_type * destination)const;
  void zero_data();
/********Indexers, getters and setters ****/
  size_t get_dims()const;
  size_t get_dims(size_t dim)const;
  
  my_type get_element(size_t nx)const;
  my_type get_element(size_t nx, size_t ny)const;
  my_type get_element(size_t nx, size_t ny, size_t nz)const;
  my_type get_element(size_t nx, size_t ny, size_t nz, size_t nt)const;
  my_type get_element(size_t n_dims, size_t * dim)const;

  size_t get_total_elements()const;
  
  bool set_element(size_t nx, my_type val);
  bool set_element(size_t nx, size_t ny, my_type val);
  bool set_element(size_t nx, size_t ny, size_t nz, my_type val);
  bool set_element(size_t nx, size_t ny, size_t nz, size_t nt, my_type val);
  bool set_element(size_t n_dims, size_t * dim, my_type val);

/********Data fillers, file IO ****/
  //Template populate for any data type
  template<typename T> bool populate_data(T dat_in, size_t n_tot);
  
  bool populate_slice(my_type * dat_in, size_t n_dims, size_t * offsets);
  bool populate_complex_slice(my_type * dat_in, size_t n_dims, size_t * offsets, size_t* sizes);

  bool write_to_file(std::fstream &file);
  bool read_from_file(std::fstream &file);
  std::vector<size_t> read_dims_from_file(std::fstream &file);
  
  bool write_section_to_file(std::fstream &file, std::vector<size_t> bounds);

/********Helpful functions working on entire array as a thing ****/
  bool resize(size_t dim, size_t sz, bool verbose=0);
  bool shift(size_t dim, long n_els);
  
  my_type minval(size_t offset=0);
  my_type maxval(size_t offset=0);
  my_type minval(std::vector<size_t> &ind, size_t offset=0);
  my_type maxval(std::vector<size_t> &ind, size_t offset=0);
  my_type partial_maxval(std::vector<std::pair<size_t, size_t> > ranges,std::vector<size_t> &ind);
  void smooth_1d(int n_pts);

  /** Apply a function to each element of array. func must take and return a my_type or type convertible to this 
  @param func Function to apply */
  void apply(std::function<my_type(my_type arg)> func){
    for(size_t i=0; i< this->get_total_elements(); i++) *(this->data+i) = func(*(this->data+i));
  }
  /** Apply a function to each element of array. func must take two my_type and return one my_type or type convertible to this
  @param func Function to apply 
  @param arg Second argument to function */
  void apply(std::function<my_type(my_type arg, my_type arg2)> func, my_type arg){
    for(size_t i=0; i< this->get_total_elements(); i++) *(this->data+i) = func(*(this->data+i), arg);
  }
  /** Fill one array from another mapping elements using function. func must take and return a my_type or type convertible to this. Each element of the rhs is tranformed with func and placed into this array 
  @param func Transform function to apply to elements 
  @param rhs Array to transform*/
  void apply(std::function<my_type(my_type arg)> func, const my_array & rhs){
    if(this->n_dims != rhs.get_dims() || this->get_total_elements() != rhs.get_total_elements()) return;
    for(size_t i=0; i< this->get_total_elements(); i++) *(this->data+i) = func(*(rhs.data+i));
  }
  /** Transform one array using another. func must take a pair of my_type and return a my_type or types convertible to this. Func is called with each element of this array and each element of rhs and the result placed back into this. @param func Function to apply @param rhs Second parameters to func*/
  void apply(std::function<my_type(my_type, my_type) > func, const my_array & rhs){
    if(this->n_dims != rhs.get_dims() || this->get_total_elements() != rhs.get_total_elements()) return;
    for(size_t i=0; i< this->get_total_elements(); i++) *(this->data+i) = func(*(this->data+i), *(rhs.data+i));
  }

friend bool populate_mirror_fastest(data_array &data_out,my_type * result_in, size_t total_els);

};

#endif /* defined(____my_array__) */
