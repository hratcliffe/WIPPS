//
//  data_array.cpp
//  
//
//  Created by Heather Ratcliffe on 03/08/2016.
//
//

#include <stdio.h>
#include <fstream>
#include <iostream>
#include <iomanip>
#include <cstdlib>
#include <algorithm>
#include <complex.h>
#include <limits.h>
#include <cmath>
#include <cstring>
#include "support.h"
#include "my_array.h"
#include "data_array.h"

/********Basic setup and allocation functions ****/

void data_array::construct(){
/** \brief Common constructor logic
*
*Sets fields for empty, dimensionless array
*/
  axes=nullptr;
  time[0]=0; time[1]=1;
  space[0]=0; space[1]=1;
  memset((void *) block_id, 0, ID_SIZE*sizeof(char));
  B_ref = 0;
  
}

void data_array::alloc_ax(const size_t els){
/** \brief Allocate axis memory
*
* Alocate memory for axes
@param els Number of elements to allocate
*/
  if(els > 0 && els <= this->n_dims*MAX_SIZE){
    axes=(my_type*)calloc((els),sizeof(my_type));
  }else{
    if(els != 0) my_error_print("Array size exceeds max. Axes alloc failed", mpi_info.rank);
    else my_error_print("Array size cannot be 0. Axes alloc failed", mpi_info.rank);
  }
}

data_array::data_array() : my_array(){
/** \brief Default constructor
*
*Create an empty, dimensionless array
*/
  construct();
}

data_array::data_array(size_t nx, size_t ny, size_t nz, size_t nt) : my_array(nx, ny, nz, nt){
/**\brief Construct a 1-4 D array
*
*Adds axes to a normal rectangular my_array of correct size
@param nx Size of x dimension
@param ny Size of y dimension, default 0
@param nz Size of z dimension, default 0
@param nt Size of t dimension, default 0
*/
  construct();
  size_t els= this->get_total_axis_elements();
  alloc_ax(els);
}

data_array::data_array(size_t n_dims, size_t * dims) : my_array(n_dims, dims){
/** \brief Construct arbitrary dimension array
*
*Adds axes to a normal rectangular my_array of correct size
@param n_dims Rank of array to create
@param dims Array of dimensions of array to create
*/
  construct();
  size_t els= this->get_total_axis_elements();
  alloc_ax(els);
}

data_array::data_array(std::string filename){
/**\brief Create data array from file
*
* Create a data array by reading from the named file. If the file does not exist nothing is done and this will be an empty array. Otherwise it reads the dimensions, sets up sizes and populates data and info 
@param filename Full path of file to read from
*/

  //Initialise empty array
  construct();

  std::fstream infile;
  infile.open(filename, std::ios::in|std::ios::binary);
  if(!infile.is_open()){
    my_error_print("File open or access error", mpi_info.rank);
    return;
  }

  //Read dims from file, unpack and allocate correct memory
  std::vector<size_t> dims_vec = my_array::read_dims_from_file(infile);
  size_t n_dims_in = dims_vec.size();
  size_t * dims_in;

  if(n_dims_in >0){
    dims_in = (size_t*)malloc(n_dims_in*sizeof(size_t));
    //Unpack vector to array
    for(size_t i=0;i<n_dims_in;i++){
      dims_in[i] = dims_vec[i];
    }
    //Allocate data space etc
    my_array::alloc_all(n_dims_in, dims_in);
    //Now we can use get_total_axis_elements to get the number of axis elements we need and then allocate
    alloc_ax(get_total_axis_elements());

    //Finally, rewind file and read in data and axes using normal routines
    infile.seekg(0, std::ios::beg);
    bool err= this->read_from_file(infile);
    if(err){
      my_error_print("IO error, could not read", mpi_info.rank);
      //Clean up and reset sizes
      if(axes) delete axes;
      if(data) delete data;
      construct();
      my_array::construct();
    }
  
  }else{
    my_error_print("Invalid dimensionality in input file", mpi_info.rank);
  }
}

data_array::~data_array(){
/**\brief Destructor 
*
*Free axis memory
**/
  if(axes) free(axes);
  axes = nullptr; // technically unnecessary as desctructor deletes members.
}

/********Technical stuff making my_array a proper "object" ****/

data_array & data_array::operator=(const data_array& src){
/** \brief Copy assignment
*
*Set this array equal to src by (deep) copying src including data
@param src Array to copy
@return Copy of array
*/
  //Trap self-assigning or bad copy before destructing
  if(&src == this || !src.is_good()) return *this;
  //Clear any existing axes and reset sizes etc
  if(this->axes) free(axes);
  this->construct();

  //Call the equivalent constructor for base class my_array
  my_array::operator=(src);

  if(this->dims){
    //Allocate axis memory
    size_t els= this->get_total_axis_elements();
    alloc_ax(els);
    //Copy axes and ids from src
    if(axes) std::copy(src.axes, src.axes+els, this->axes);
    copy_ids(src);
  }
  return *this;

}

data_array::data_array(data_array && src) : my_array(src){
/** \brief Move constructor
*
*Move src to new location. Copies data pointers but does not reallocate memory. Src is left empty
@param src Array to copy
*/
  this->axes = src.axes;
  src.axes = nullptr;
  copy_ids(src);
}

data_array::data_array(const data_array &src) : my_array(src){
/** \brief Copy constructor
*
*Copy src to a new instance, making a duplicate of data
@param src Array to copy
*/
  construct();
  //Basic construction of additionals, already called base class copy constructor
  //Allocate axes, copy data and ids
  size_t els= this->get_total_axis_elements();
  alloc_ax(els);
  if(axes) std::copy(src.axes, src.axes+els, this->axes);
  copy_ids(src);
}

bool data_array::operator==(const data_array &rhs)const{
/** \brief Equality operator
*
* Check this is equal to rhs. Since copies are always deep, we check values, not data pointers
@param rhs Array to compare to
@return Boolean true if equal, false else
*/
  if(!my_array::operator==(rhs)) return false;

  //Check axes elementwise
  for(size_t i=0; i< this->get_total_axis_elements(); i++) if(*(this->axes + i) != *(rhs.axes + i)) return false;
  if(this->check_ids(rhs)) return false;

  return true;
}

data_array & data_array::operator=(const my_array& src){
/** \brief Conversion equality operator
*
*Convert a my_array into a data array. Data is deep copied, and axes are added
@param src My_array object to copy data and sizes from
@return A data array containing copy of data in src plus empty axes
*/

  //Trap self-assigning or bad copy before destructing
  if(&src == this || !src.is_good()) return *this;
  //Clear any existing axes and reset sizes etc
  if(this->axes) free(axes);
  this->construct();

  //Call the equivalent constructor for base class my_array
  my_array::operator=(src);

  if(this->dims){
    //Allocate axis memory
    size_t els= this->get_total_axis_elements();
    alloc_ax(els);
    //Axes and ids are present but empty because of construct and alloc_ax calls
  }
  return *this;

}

data_array::data_array(const my_array & src) : my_array(src){
/** \brief Conversion operator
*
*Convert a my_array into a data array. Data is deep copied, and axes are added
@param src My_array object to copy data and sizes from
*/
  construct();
  //Basic construction of additionals, already called base class copy constructor
  //Allocate axes, copy data and ids
  size_t els= this->get_total_axis_elements();
  alloc_ax(els);
}

/********Helpers for working with data_array ****/

my_type * data_array::disown_axes(){
/** \brief Disown and return axes pointer
*
*Surrenders ownership of memory pointed to by axes NB if this pointer is not kept an manually freed, memory will leak. To aquire ownership of both the axes and data, use disown_axes and my_array::disown_data in any order, but note that after the latter, this will not be a valid array and getter/setter for data and axes will return nothing.
@return Pointer to the axis memory
*/
  my_type * data_pointer = this->axes;
  this->axes = nullptr;
  return data_pointer;
}

void data_array::clone_empty(const data_array &src){
/** \brief Initialise this to match sizes of src
*
*This will be a valid empty array of size matching src.
@param src Array to copy from
*/
  if(axes) free(axes);
  construct();
  my_array::clone_empty(src);
  size_t els= this->get_total_axis_elements();
  alloc_ax(els);

}

void data_array::copy_ids( const data_array &src){
/** Copies ID fields from src array to this 
@param src Array to copy from
*/

  //Strncpy in case the src field does not actually contain a null-terminated string
  strncpy(this->block_id, src.block_id, ID_SIZE);
  for(size_t i=0; i < 2; ++i) this->time[i] = src.time[i];
  for(size_t i=0; i < 2; ++i) this->space[i] = src.space[i];
  this->B_ref = src.B_ref;
}

bool data_array::check_ids( const data_array & src)const{
/** Checks ID fields match src 
@param src Array to check against
@return False for non-matching ids, true for complete match
*/

  bool err = false;
  if(strcmp(this->block_id, src.block_id) != 0) err = true;
  for(size_t i=0; i< 2; i++) if(src.time[i] != this->time[i]) err=true;
  for(size_t i=0; i < 2; ++i) if(this->space[i] != src.space[i]) err=true;
  if(this->B_ref != src.B_ref) err = true;
  return err;
}

/********Indexers, getters and setters ****/

size_t data_array::get_total_axis_elements()const{
/** \brief Return total axes length
*
*
@return The number of elements in all axes combined
*/
  size_t tot_els=0;
  for(size_t i=0; i<n_dims;++i) tot_els +=dims[i];
  return tot_els;

}

long data_array::get_axis_index(size_t dim, size_t pt)const{
/** \brief Get index of axis element
*
*Get the position in the backing data store of point pt on dimension dim. Takes care of all bounds checking and disposition in memory
@param dim Dimension of axis
@param pt Index into required axis
@return Index into 1-D array of pt, or -1 for any sort of out-of-range issue
*/

  if(dim >=n_dims || pt >=get_dims(dim)) return -1;
  //Out of range error
  
  long offset = 0;
  //Axes are stored sequentially in 1-D layout
  for(size_t i=0; i< dim; i++) offset +=dims[i];
  return offset + pt;

}

float data_array::get_res(size_t i)const{
/**\brief Get (linear) axis resolution
*
*Return resolution of axis on dimension i. Takes the total length and divides by the number of elements to avoid rounding errors. If axis is not linear, this is meaningless. 
@param i Dimension of axis
@return 0.0 if the axis is undefined or 0 or 1 in length, otherwise the calculated linear resolution
*/
  if(dims[i] > 1) return std::abs(get_axis_element(i, dims[i]-1) - get_axis_element(i, 0))/(float)(dims[i]-1);
  else return 0.0;

}

std::vector<size_t> data_array::get_bounds(std::vector<my_type> limits){
/** \brief Convert axis values to indices
*
* For a vector of bounds, 2 per dimension, convert the required axis bounds into index bounds
@param limits The real physical axis values
@return The corresponding index values
*/

  std::vector<size_t> index_limits;

  if(limits.size() != 2*n_dims){
    //Wrong size, report and return empty vector
    my_error_print("Limits vector size does not match array!", mpi_info.rank);
    return index_limits;
  }

  index_limits.resize(n_dims*2);

  size_t len;
  const my_type * ax_start;
  long where_val;
  //Lookup the bounds one dimension at a time.
  for(size_t i=0; i< n_dims; i++){
    //Set defaults for if bounds are out of range
    index_limits[i*2] = 0;
    index_limits[i*2 + 1] = this->dims[i]-1;
    //Now attempt lookups and set values if in range
    ax_start = get_axis(i, len);
    where_val = where(ax_start, len, limits[2*i]);
    if(where_val != -1) index_limits[i*2] = where_val;
    where_val = where(ax_start, len, limits[2*i + 1]);
    if(where_val != -1) index_limits[i*2 + 1] = where_val;
  }

  return index_limits;
}

my_type data_array::get_axis_element(size_t dim, size_t pt)const{
/** \brief Get axis value
*
* 
@param dim Dimension
@param pt Point to get
@return value at pt in dimension dim if in range, else 0.0
*/

  long ind = get_axis_index(dim, pt);
  if(ind  != -1){
    return axes[ind];
  }else{
    return 0.0;
  }

}

bool data_array::set_axis_element(size_t dim, size_t pt, my_type val){
/** \brief Sets axis element
*
*Sets elements at pt on dimension dim, 
@param dim Dimension
@param pt Point to set
@param val Value to set
@return 1 if out of range, 0 else.
*/

  long index = get_axis_index(dim, pt);
  if(index >= 0){
    axes[index] = val;
    return 0;
  }else{
    return 1;
  }

}

const my_type * data_array::get_axis(size_t dim, size_t & length){
/**  \brief Get pointer to axis
*
*Returns pointer to given axis and its length. NB do not muck with this pointer. It's provided for ease of using where but is a const my_type * for a reason
@param dim Dimension of axis to get
@param[out] length Returns the length of axis
@return Pointer to start of axis, or nullptr if axis doesn't exist or requested dim is out of range
*/

  if(!axes || (dim >= n_dims)) return nullptr;

  long index = get_axis_index(dim, 0);
  //Get index of 0th element
  length = get_dims(dim);
  if(index != -1) return axes + index;
  else return nullptr;

}

long data_array::get_axis_index_from_value(size_t dim, my_type value)const{
/** Get the index on axis dim of value value
*
@param dim Dimension for lookup
@param value Value to find
@return Index of value in axis for dimension dim
*/
  size_t len = dims[dim];
  long index = get_axis_index(dim, 0);
  if(index != -1){
    long where_val = where(axes+index, len, value);
    return where_val;
  }else{
    return -1;
  }
}

/********Data/axis fillers, file IO ****/

bool data_array::populate_axis(size_t dim, my_type * dat_in, size_t n_tot){
/** \brief Fill axis from dat_in
*
*Populates axis from dat_in. Number of elements copied will be the smaller of n_tot and size of dimension dim.
@param dim Dimension of axis to fill
@param dat_in Pointer to start of data to copy
@param n_tot Size of input data array
@return 0 (success) 1 (error)
*/

  if(dim >=n_dims) return 1;
  //Out of range error

  //Min of n_tot and tot_els
  size_t tot_els = dims[dim];
  if(n_tot < tot_els) tot_els = n_tot;
  size_t len = 0;
  my_type * ax_ptr = const_cast<my_type *>(this->get_axis(dim, len));
  std::copy(dat_in, dat_in + tot_els, ax_ptr);

  return 0;

}

void data_array::make_linear_axis(size_t dim, float res, long offset){
/**\brief Make an axis
*
*Generates a linear axis for dimension dim, with resolution res, starting at value of  - offset*res That allows one to guarantee that 0 appears regardless of the number of cells. 
@param dim Dimension to build axis for 
@param res Axis resolution 
@param offset Number of grid cells to shift downwards (leftwards) by
*/

  size_t len = get_dims(dim);
  for(size_t i = 0; i < len; i++){
    set_axis_element(dim, i, ((float) ((long)i - offset)) * res);
  }
}


//Dummy to document file format
/**
    * @class dummy_data_file_format
    * After the layout in my_array::write_to_file is added
    * Next_block axes
    *"Footer:"
    * Next_block time[2] space[2] B_ref 
    * Next_block Block_id
    * If file is to be closed then we finish with Footer_start, i.e. the position of the start of "Footer" section
    *If multiple arrays are written, we write each without closing and the close the file with a final Block_id tag

*/

bool data_array::write_to_file(std::fstream &file, bool close_file){
/** \brief Write data array to file
* \copydoc dummy_data_file_format  
@param file Filestram to write to
@param close_file Whether to write file final data
@return 0 (success), 1 (error)
* \todo Test read/write with double, also IDL
*/

  if(!file.is_open()) return 1;
  bool write_err=0;

  //Call base class method to write the my_array data
  my_array::write_to_file(file);

  //First write the axes
  //Get position of next section
  size_t tot_ax = get_total_axis_elements();
  size_t next_location = (size_t) file.tellg() + sizeof(size_t) + tot_ax*sizeof(my_type);
  file.write((char*) & next_location, sizeof(size_t));
  //Write axis data
  file.write((char *) axes , sizeof(my_type)*tot_ax);
  //Check position is where we expect
  if((size_t)file.tellg() != next_location) write_err=1;

  //Grab the start position of the "footer" block
  size_t ftr_start = next_location;

  //Write the "footer" data in two sections
  next_location += sizeof(my_type)*3 + sizeof(size_t)*2 + sizeof(size_t);
  file.write((char*) & next_location, sizeof(size_t));

  file.write((char *)time, sizeof(my_type)*2);
  file.write((char *)space, sizeof(size_t)*2);
  file.write((char *) &B_ref, sizeof(my_type));

  next_location += sizeof(char)*ID_SIZE + sizeof(size_t);
  file.write((char*) & next_location, sizeof(size_t));
  file.write(block_id, sizeof(char)*ID_SIZE);

  if((size_t)file.tellg() != next_location) write_err=1;
  
  //Finish with position of start of footer if requested
  if(close_file) file.write((char*) & ftr_start, sizeof(size_t));

  //Report on any write errors
  if(write_err) my_error_print("Error writing offset positions", mpi_info.rank);
  return 0;

}

bool data_array::read_from_file(std::fstream &file){
/** \brief Read data array from file
*
*Reads data from file. This array should have already been created in the correct shape, otherwise we return an error. \copydoc dummy_data_file_format
@param file Filestram to read from
@return 0 (success), 1 (error)
*/

  bool err=false;

  //Call my_array for initial reads
  if(file.good()) err = my_array::read_from_file(file);
  if(err){
    my_error_print("File read failed", mpi_info.rank);
    return err;
  }

  //Read the axis elements
  size_t next_block = 0;
  file.read((char*) &next_block, sizeof(size_t));
  file.read((char *) this->axes , sizeof(my_type)*get_total_axis_elements());


  //Read the footer data
  file.read((char*) &next_block, sizeof(size_t));
  file.read((char *) time, sizeof(my_type)*2);
  file.read((char *) space, sizeof(size_t)*2);
  file.read((char *) &B_ref, sizeof(my_type));

  file.read((char*) &next_block, sizeof(size_t));
  char id_in[ID_SIZE];
  if(file){
    file.read(id_in, sizeof(char)*ID_SIZE);
    strncpy(this->block_id, id_in, ID_SIZE);
  }
  return err;

}

bool data_array::write_section_to_file(std::fstream &file, std::vector<my_type> limits, bool close_file){
/** \brief Write array section to file
*
*Write section between given AXIS values to file. To use one dimension entire supply values less/greater than min and max axis values.
@param file Filestram to write to
@param limits Vector of min and max physical axis values for section. Must be 2 per dimension
@param close_file Whether to write file final data
@return 0 (success), 1 (error)
*/
  //Identify limits of segment from axes
  std::vector<size_t> index_limits = this->get_bounds(limits);
  return write_raw_section_to_file(file, index_limits, close_file);

}

bool data_array::write_raw_section_to_file(std::fstream &file, std::vector<size_t> index_limits, bool close_file){
/** \brief Write section of array to file
*
*Write section defined by the limits to supplied file. If limits are out of range then the entire dimension is used. See data_array::write_to_file for file details
@param file Filestram to write to
@param index_limits Vector of min and max axis indices for section. Must be 2 per dimension
@param close_file Whether to write file final data
@return 0 (success), 1 (error)
*/

  if(!file.is_open()) return 1;
  bool write_err = 0;

  if(index_limits.size() != 2*n_dims){
    my_error_print("Limits vector size does not match array!", mpi_info.rank);
    return 1;
  }

  for(size_t i=0; i< index_limits.size(); i+=2){
    //Make sure limits are within ranges, adjust if necessary
    if(index_limits[i] >= dims[i/2]) index_limits[i] = 0;
    if(index_limits[i+1] >= dims[i/2]) index_limits[i+1] = dims[i/2];
  }

  //Write initial data using my_array function
  my_array::write_section_to_file(file, index_limits);

  //Calculate total number of axis elements to be written
  size_t tot_ax = 0;
  for(size_t i=0; i< n_dims; i++){
    tot_ax +=(index_limits[2*i +1]-index_limits[2*i]);
  }
  
  //Get and write position of next section
  size_t next_location = (size_t) file.tellg() + sizeof(size_t) + tot_ax*sizeof(my_type);
  file.write((char*) & next_location, sizeof(size_t));

  //Write axis data
  size_t len;
  for(size_t i=0; i< n_dims; i++){
    file.write((char *) (get_axis(i, len)+index_limits[2*i]), sizeof(my_type)*(index_limits[2*i +1]-index_limits[2*i]));
  }
  if((size_t)file.tellg() != next_location) write_err=1;

  //Write footer as in write_to_file
  size_t ftr_start = next_location;

  next_location += sizeof(my_type)*3 + sizeof(size_t)*2 + sizeof(size_t);
  file.write((char*) & next_location, sizeof(size_t));
  file.write((char *) time, sizeof(my_type)*2);
  file.write((char *) space, sizeof(size_t)*2);
  file.write((char *) &B_ref, sizeof(my_type));

  next_location += sizeof(char)*ID_SIZE + sizeof(size_t);
  file.write((char*) & next_location, sizeof(size_t));
  file.write(block_id, sizeof(char)*ID_SIZE);

  if((size_t)file.tellg() != next_location) write_err=1;

  //Close file if requested
  if(close_file) file.write((char*) & ftr_start, sizeof(size_t));

  if(write_err) my_error_print("Error writing offset positions", mpi_info.rank);

  return 0;
}

bool data_array::write_closer(std::fstream &file){
/** \brief Write close tag of file
*
*Writes the final footer into a file 
@param file Filestream to write to
@return 0 (success), 1 (error)
*/

  if(!file.is_open()) return 1;
  bool write_err = 0;

  size_t ftr_start = (size_t) file.tellg();
  //Start of ftr means where to start reading footer block, i.e. location of the next_location tag
  size_t next_location = ftr_start + sizeof(my_type)*3 + sizeof(size_t)*2 + sizeof(size_t);
  file.write((char*) & next_location, sizeof(size_t));
  file.write((char *)time, sizeof(my_type)*2);
  file.write((char *)space, sizeof(size_t)*2);
  file.write((char *) &B_ref, sizeof(my_type));

  next_location += sizeof(char)*ID_SIZE +2*sizeof(size_t);
  //Position of next section, i.e. of file end val
  file.write((char*) & next_location, sizeof(size_t));

  file.write(block_id, sizeof(char)*ID_SIZE);
  file.write((char*) & ftr_start, sizeof(size_t));

  //Check for error and report
  if((size_t)file.tellg() != next_location) write_err=1;
  if(write_err) my_error_print("Error writing offset positions", mpi_info.rank);

  return write_err;
}

/********Helpful functions working on entire array as a thing ****/

bool data_array::resize(size_t dim, size_t sz, bool verbose){
/** \brief Resize my_array on the fly
*
*If sz < dims[dim] the first sz rows will be kept and the rest deleted. If sz > dims[dim] the new elements will be added zero initialised. Similarly for axis elements. See my_array::resize() for more.
@param dim The dimension to resize
@param sz The new size
@param verbose Whether to print some info
@return 0 (success), 1 (error)
*/

  size_t old_sz = this->get_dims(dim);
  //Call my_array::resize to resize data portion
  bool err = my_array::resize(dim, sz, verbose);

  if(!err){
    //Success! Do axes!
    my_type * new_ax;
    size_t part_sz = 0;

    if(dim == n_dims-1){
      //Special case as we can shrink and maybe grow without copy
      //New size is sum of all other dims + sz
      for(size_t i=0; i<dim; ++i) part_sz+= dims[i];
      new_ax = (my_type *) realloc((void*) this->axes, (part_sz+sz)*sizeof(my_type));
      if(!new_ax){
        //Failure. Axes now don't match data array!
        //This is bad and will do weird out of range stuff if left
        my_error_print("Failed to reallocate axes", mpi_info.rank);
        this->axes = nullptr;
        return 1;
      }
      //Zero the newly allocated elements
      long new_els = (sz - dims[dim]);
      if(new_els > 0) memset((void*)(new_ax + part_sz), 0, new_els*sizeof(my_type));
      //Set axes to point to new block
      axes = new_ax;
    }else{
      //We have to allocate a new block and copy across
      //New size is sum of all other dims + sz
      for(size_t i=0; i<dim; ++i) part_sz += dims[i];
      for(size_t i=dim+1; i<n_dims; ++i) part_sz += dims[i];

      //Allocate and copy
      new_ax=(my_type*)calloc(part_sz+sz,sizeof(my_type));
      size_t els_to_copy = 0, old_starts=0, new_starts=0;
      for(size_t i=0; i<n_dims; ++i){
        els_to_copy = dims[i];
        if(i == dim && sz < old_sz) els_to_copy = sz;
        std::copy(axes+old_starts, axes+old_starts+els_to_copy, new_ax+new_starts);
        //Advance the offsets by correct amounts for old and new copies
        if(i!= dim) old_starts += dims[i];
        else old_starts +=old_sz;
        new_starts +=els_to_copy;
      
      }
      //Free old block and point axes to new block
      free(axes);
      axes = new_ax;

    }
    return 0;
  }else{
    return 1;
  }

}

bool data_array::shift(size_t dim, long n_els, bool axis){
/** Shift array on dim dim by n_els
*
*Shift is cyclical 
@param dim Dimension to shift on 
@param n_els Number of elements to shift by 
@param axis Whether to shift the corresponding axis 
@return 0 (success), 1 (error)
*/
  if(dim >= n_dims) return 0;

  //my_array takes care of data shifting
  bool err = my_array::shift(dim, n_els);

  if(axis && !err){
    //Rotate the specified axis if we succeeded
    //Adjust shift amount to be positive. This is always possible for cyclical shift
    long sign_n = (n_els<0 ? -1: 1);
    n_els = sign_n >0? (std::abs(n_els)%dims[dim]) : dims[dim]-(std::abs(n_els)%dims[dim]) ;
    //By now n_els guaranteed >= 0

    size_t len=0;
    const my_type * ax = get_axis(dim, len);
    my_type * new_axis = (my_type*)malloc(len*sizeof(my_type));
    if(!new_axis){
      //Since the data shift succeeded, this should never happen in practise
      my_error_print("Error allocating spare memory for shift", mpi_info.rank);
      return 1;
    }
    long actual_shift = 0;
    //Copy out axis, rotate, and copy back
    std::copy(ax,ax+len, new_axis);
    for(size_t j=0; j< len; j++){
      //This is inefficient but is minimmaly changed from shift for whole array
      actual_shift = (j+n_els >= dims[dim])? (j+n_els-dims[dim]): j+n_els;
      actual_shift += (actual_shift < 0 ? dims[dim]:0);
      *(new_axis + j) = *(ax + actual_shift);
//      std::copy(new_axis+j, new_axis+(j+1), ax+actual_shift);
    }
    //Clean up
    if(new_axis) free(new_axis);
  }
  return err;

}

data_array data_array::total(size_t dim){
/** \brief Total array on dim dim
*
* Wraps function data_array::total(size_t dim, size_t min, size_t max).
@param dim Dimension to total on
@return A new array of rank n_dims -1, containing data summed over entire range of dimension dim. If dim is out of range, empty array is returned
*/
  //No dimension to total, return empty array
  if(dim >= this->n_dims) return data_array();

  //Uses total with ranging, passing end values to capture entire range
  return this->total(dim, (size_t) 0, this->get_dims(dim));

}

data_array data_array::total(size_t dim, my_type min, my_type max){
/** \brief Total array on dim dim
*
*Wraps function data_array::total(size_t dim, size_t min, size_t max).
@param dim Dimension to total on
@param min Minimum physical axis value of slice to total. Note this expects a monotonic axis.
@param max Maximum physical axis value of slice to total.
@return A new array of rank n_dims -1, containing data summed over entire range of dimension dim. If dim is out of range, empty array is returned
*/

  //No dimension to total, return empty array
  if(dim >= this->n_dims) return data_array();

  size_t len;
  const my_type * ax_start = get_axis(dim, len);
  if(ax_start){
    //Start with indices set to ends of axis
    size_t min_ind = 0, max_ind = this->get_dims(dim);
    //Lookup provided range and update indices if within axis range
    int where_val = -1;
    where_val = where(ax_start, len, min);
    if(where_val != -1) min_ind = where_val;
    where_val = where(ax_start, len, max);
    if(where_val != -1) max_ind = where_val;
    //Call total function with indices
    return total(dim, min_ind, max_ind);
  }else{
    return data_array();
  }

}

data_array data_array::total(size_t dim, size_t min_ind, size_t max_ind){

/** \brief Total along dim dim
*
*
@param dim Dimension to total on
@param min_ind Minimum axis index of slice to total.
@param max_ind Maximum axis index of slice to total.
@return A new array of rank n_dims -1, containing data summed over entire range of dimension dim. If dim is out of range, empty array is returned. Note totalling a 1-d array gives a 1-element array
*/

  //No dimension to total, return empty array
  if(dim >= this->n_dims) return data_array();
  
  //Adjust indices to be within range
  //If min is big we've probably got an overflow...
  if(min_ind >= get_dims(dim)) min_ind = 0;
  if(max_ind >= get_dims(dim)) max_ind = get_dims(dim);
  
  //Create array of the sizes after totalling
  size_t * new_dims;
  new_dims = (size_t *) calloc((this->n_dims-1), sizeof(size_t));
  for(size_t i=0, i2=0; i< this->n_dims; i++, i2++){
    if(i == dim){
      i2--;
      continue;
    }
    new_dims[i2] = dims[i];
  }

  //Create the new, smaller, array
  size_t new_n_dims = 0;
  data_array new_array;
  if(this->n_dims > 1){
    new_n_dims = this->n_dims - 1;
    new_array = data_array(new_n_dims, new_dims);
  }
  else{
    //Return 1-element array rather than actual scalar
    new_dims = (size_t *) malloc(sizeof(size_t));
    new_dims[0] = 1;
    new_n_dims = 1;
    new_array = data_array(new_n_dims, new_dims);
  }
  //Flatten the actual data
  flatten_fortran_slice(this->data, new_array.data, this->n_dims, this->dims, dim, min_ind, max_ind);
  //Copy axes and ids
  size_t len;
  //One element array has no axes
  if(this->n_dims > 1){
    for(size_t i = 0, i2 = 0; i2 < new_n_dims; i++, i2++ ){
      if(i == dim) i++;
      len = new_array.get_dims(i2);
      for(size_t el = 0; el < len; el ++) new_array.set_axis_element(i2, el, get_axis_element(i, el));
    }
  }
  new_array.copy_ids(*this);

  free(new_dims);
  return new_array;
}

data_array data_array::average(size_t dim){
/** \brief Average array over dim dim
*
*We guarantee to match total's behaviour if dim is out of range.
@param dim Dimension to average over
@return New array of rank n_dims-1 filled with the average values over the given dim.
*/
  size_t old_dim = this->get_dims(dim);
  data_array new_arr = this->total(dim);
  //If dim doesn't exist we match total's behaviour, but don't want a divide by 0 error
  if(old_dim > 0) new_arr.apply(divide, old_dim);
  return new_arr;
}

data_array data_array::average(size_t dim, my_type min, my_type max){
/** \brief Average array over dim dim
*
*We guarantee to match total's behaviour if dim is out of range.
@param dim Dimension to average over
@param min Minimum axis value emcompassing desired slice
@param max Maximum axis value emcompassing desired slice
@return New array of rank n_dims-1 filled with the average values over the given dim.
*/

  //Guarantee to match behaviour of total by calling it right now
  if(dim >= this->n_dims) return total(dim, min, max);

  size_t len, min_ind=0, max_ind = this->get_dims(dim)-1;
  int where_val;
  const my_type * ax_start = get_axis(dim, len);
  where_val = where(ax_start, len, min);
  if(where_val != -1) min_ind = where_val;
  where_val = where(ax_start, len, max);
  if(where_val != -1) max_ind = where_val;

  return average(dim, min_ind, max_ind);
}

data_array data_array::average(size_t dim, size_t min_ind, size_t max_ind){
/** \brief Average array over dim dim
*
*We guarantee to match total's behaviour if dim is out of range.
@param dim Dimension to average over
@param min Minimum axis index emcompassing desired slice
@param max Maximum axis index emcompassing desired slice
@return New array of rank n_dims-1 filled with the average values over the given dim.
*/

  //Guarantee to match behaviour of total by calling it right now
  if(dim >= this->n_dims) return total(dim, min_ind, max_ind);

  //If dim is OK, then re-range the indices
  if(min_ind >= this->get_dims(dim)) min_ind = 0;
  if(max_ind >= this->get_dims(dim)) max_ind = this->get_dims(dim)-1;

  data_array new_arr = total(dim, min_ind, max_ind);
  if(max_ind != min_ind) new_arr.apply(divide, max_ind - min_ind);
  return new_arr;
}

/********FFT functions and helpers ****/
bool fft_array(const data_array &data_in, data_array & data_out){
#ifndef NO_FFT
/** \brief FFT data_array
*
* Data and axes in this object are FFT'd using FFTW and stored into the instance pointed to by data_out. Data_out must be created with correct dimensions first, but we check and return error (1) if it is not so. The returned values are the abs-square of the FFT.
*/

  //Check output array is defined and has correct dimensions
  if(!data_out.is_good()){
    my_error_print("Output array for FFT undefined", mpi_info.rank);
    return 1;
  }
  if(data_out.get_dims() != data_in.get_dims()){
    my_error_print("Wrong output dimensions for FFT", mpi_info.rank);
    return 1;
  }
  for(size_t i=0; i<data_in.get_dims();++i){
    if(data_out.get_dims(i) != data_in.get_dims(i)){
      my_error_print("Wrong output dimensions for FFT", mpi_info.rank);
      return 1;
    }
  }

  data_out.copy_ids(data_in);
  size_t total_size = 1; /* Total number of elements in array*/
  for(size_t i=0; i<data_in.get_dims();++i) total_size *= data_in.get_dims(i);

  //Do the FFTs according to FFTWs requirements
  ADD_FFTW(plan) p;
  cplx_type *out;
  my_type * in, *result;
  size_t dim_0 = data_in.get_dims(0);

  size_t output_size=data_in.get_total_elements()/dim_0*(dim_0/2+1);
  //Size of r2c transform output
  //This has been checked manually with valgrind and is CORRECT

  //Allocate memory for transforms
  //ADD_FFTW is a macro taking care of proper typing
  in = (my_type*) ADD_FFTW(malloc)(sizeof(my_type) * total_size);
  out = (cplx_type *) ADD_FFTW(malloc)(sizeof(cplx_type) * output_size);

  result = (my_type*) ADD_FFTW(malloc)(sizeof(my_type) * output_size);

  //Construct dims array in reverse order as we're using the opposite majority to FFTW
  int *rev_dims;
  rev_dims = (int *) malloc(data_in.get_dims()*sizeof(size_t));
  for(size_t i=0; i< data_in.get_dims(); i++) rev_dims[i] =data_in.get_dims(data_in.get_dims()-1-i);

  //Create the "plan"
  p = ADD_FFTW(plan_dft_r2c)((int)data_in.get_dims(), rev_dims, in, out, FFTW_ESTIMATE);
  free(rev_dims);

  //Copy data into in. Because the plan creation changes in, so we don't want to feed our actual data array in, and it's safer to run the plan with the memory block it was created with
  data_in.copy_data(in);
  
  //Execute the plan
  ADD_FFTW(execute)(p);

  cplx_type * addr;
  addr = out;
  //because double indirection is messy and cplx type is currently a 2-element array of floats

  for(size_t i=0; i< output_size; i++){
    *(result+i) = (my_type)(((*addr)[0])*((*addr)[0]) + ((*addr)[1])*((*addr)[1]));
    addr++;
  }
  //Absolute square of out array to produce final result of type my_type

  ADD_FFTW(free)(in);
  ADD_FFTW(free)(out);
  //Free as early as possible...
  
  bool err = false;
  err = populate_mirror_fastest(data_out, result, total_size);
  //Copy result into out array, restoring the lacking dimensions due to FFTW optimising

  if(err){
    my_error_print("Error populating result array", mpi_info.rank);
    return 1;
  }

  ADD_FFTW(destroy_plan)(p);
  ADD_FFTW(free)(result);
  //Destroy stuff we don't need

  float N2, res;
  for(size_t i = 0; i < data_out.get_dims(); ++i){
  //Loop over the dimensions and construct each axis in place. We KNOW now that data_out has correct dimensions. We checked before getting here.
    size_t dim_i = data_out.get_dims(i);
    N2 = ((float) dim_i)/2.0;
    res = data_in.get_res(i);
    for(size_t j = 0; j < dim_i; j++) data_out.set_axis_element(i, j, pi * ((float)j -N2)/N2/res);
  }

  return 0;

#endif
return 1;
//Return for noFFT version
}

bool populate_mirror_fastest(data_array &data_out,my_type * result_in, size_t total_els){
/** \brief Copy FFTW data into array

*
*For real data an FFT has a redundant half so FFTW returns array od size (dims[0]/2 +1)*dims[...]*dims[n-1] (Note that our first dim is FFTWs last). Since we want all k we have to mirror this to obey H(-f, -g) = H(f, g). For simplicity we enforce this on k_x, omega pair only. Data is assumed unshifted. Should work for 1-3 dimensions. 1st dimension will be returned shifted to 0-in-centre \todo Which side is negative k? \todo Can we replace reverse_copy with slice fillers and remove friendship?
*/
  size_t dim_0 = data_out.get_dims(0);
  size_t last_size = floor(dim_0/2) + 1;
  size_t num_strides = total_els/dim_0;//Always integer

  size_t num_seconds = data_out.get_dims() >= 2? data_out.get_dims(1): 1;//Size of "second dimension"
  size_t slice_sz = num_seconds*dim_0;

  //Zero data_out in case it was already populated
  data_out.zero_data();

  for(size_t i=0; i< num_strides; i++){
    //Reverse the result into place (left side of a 2-d)
    std::reverse_copy(result_in+ i*last_size, result_in+(i+1)*last_size-1, data_out.data+i*dim_0+1);
  }

  //Now fill the redundant 2-d side. Zero stays as zero, we flip the rest. This is to give us "correct" 2-D result when we take only +ve omega
  for(size_t j=0; j< num_strides/num_seconds; j++){
    std::reverse_copy(data_out.data+1+j*slice_sz, data_out.data+last_size+j*slice_sz, data_out.data+last_size-1+j*slice_sz);
    for(size_t i=1; i< num_seconds; i++){
      std::reverse_copy(data_out.data+i*dim_0+1+j*slice_sz, data_out.data+i*dim_0+last_size+j*slice_sz, data_out.data+j*slice_sz+(num_seconds-i)*dim_0+last_size-1);
    }
  }
  
  //do shifts for all but 0th dim, this is already done by steps above
  size_t shft = 0;
  for(size_t i=1; i< data_out.get_dims(); i++){
    shft = data_out.get_dims(i)/2;
    data_out.shift(i, shft, 0);
  }

  
  if(data_out.get_dims() >=3){
    //Reverse omega vs k_x on -ve k_x side only.
    size_t dim_1 = data_out.get_dims(1), dim_2=data_out.get_dims(2);
    for(size_t i=0; i<dim_0; i++){
      for(size_t j=0; j<dim_1; j++){
        for(size_t k=0; k<dim_2/2; k++){
        
          *(data_out.data+ i + j*dim_0+(dim_2-1-k)*dim_1*dim_0) = *(data_out.data+ i + j*dim_0+ k*dim_1*dim_0);
        }
      }
    }
  
  }
  
  return 0;
}

/********Other non-member helpers ****/
my_type avval(const data_array & array_in){
/** Find average value of array data*/
  size_t n_dims = array_in.get_dims();
  data_array new_array = array_in;
  for(size_t i = n_dims; i!=0; i--){
    new_array = new_array.average(i-1);
  }
  //Averaged on all dims so is now 1-d, but check to be sure
  if(new_array.get_dims() != 1){
    my_error_print("Error calculating average value", mpi_info.rank);
  }
  return new_array.get_element((size_t) 0);
}


