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
*Contains dimension information and data. Can be rectangular of 2,3,or 4 dimensions or ragged of 2 (rows of different lengths). Get_index and get_total_elements account for all details of internal layout in memory. Everything else should go via these as much as possible. NOTE this is completely C style. But contigous memory and pointer arithmetic give major speed advantage and we very rarely change size on the fly. However nothing outside this class should need to do anything except access by index and populate by element, row or entire.
 \author Heather Ratcliffe \date 21/09/2015
*/

class my_array{

protected:
  int n_dims;/**< Number of dimensions*/
  int * dims;/**< Array dimensions*/
  bool ragged; /**< Flag for if array has different lengths on each row...*/
  int * row_lengths;/**< Row lengths if ragged*/
  int * cumulative_row_lengths; /**< Offset to start of given row to iterate through*/
  my_type *data;/**< The data */
  bool defined; /**< Flag to check memory allocation sucess*/

public:

  my_array();
  virtual void construct();
  my_array(int nx, int ny, int nz=0, int nt=0);
  my_array(int * row_lengths, int ny);
  virtual ~my_array();

  virtual bool is_good(){return !defined;}/**< Check memory allocation etc worked*/
  virtual int get_index(int nx);
  virtual int get_index(int nx, int ny);
  virtual int get_index(int nx, int ny, int nz);
  virtual int get_index(int nx, int ny, int nz, int nt);
  int get_total_elements();

  int get_dims();
  int get_dims(int dim);
  int get_length(int dim);
  my_type get_element(int nx);
  my_type get_element(int nx, int ny);
  my_type get_element(int nx, int ny, int nz);
  my_type get_element(int nx, int ny, int nz, int nt);
  bool set_element(int nx, int ny, my_type val);
  bool set_element(int nx, int ny, int nz, my_type val);
  bool set_element(int nx, int ny, int nz, int nt, my_type val);

  bool populate_data(my_type * dat_in, int n_tot);
  bool populate_row(void * dat_in, int nx, int y_row);

  virtual bool write_to_file(std::fstream &file);
  bool read_from_file(std::fstream &file, bool no_version_check=0);
  bool write_section_to_file(std::fstream &file, std::vector<int> bounds);
  virtual bool resize(int dim, int sz);
  
  //void tmp_function();
  
};

/** \brief Extended my_array class including axes
*
*Contains also data axes and various id's describing the data set from which data came. Also has ability to fft itself into a new instance. \author Heather Ratcliffe \date 21/09/2015
*/
class data_array : public my_array{

protected:

  bool ax_defined;/**< Flag showing whether axes are fully defined*/
  my_type *axes;/**< 1-d array in sections, so can be arbitary length and dims*/

public:

  char block_id[ID_SIZE]; /**< The field name id from SDF file*/

  float time[2];/**< Time range over which data are taken*/
  int space[2];/**< Space range over which data are taken*/

  virtual void construct();
  data_array(int nx, int ny);
  data_array(int nx, int ny, int nz);
  data_array(int nx, int ny, int nz, int nt);
  data_array(int * row_lengths, int ny);
  data_array(std::string filename, bool no_version_check = 0);
  virtual ~data_array();

  virtual bool is_good(){return defined && ax_defined;}
  my_type * get_axis(int dim, int & length);

  int get_axis_index(int dim, int pt);

  my_type get_axis_element(int dim, int pt);
  bool set_axis_element(int dim, int pt, my_type val);
  bool populate_axis(int dim, my_type * dat_in, int n_tot);
  float get_res(int i);
  int get_total_axis_elements();

  void make_linear_axis(int dim, float res, int offset=0);

  bool write_to_file(std::fstream &file);
  bool read_from_file(std::fstream &file, bool no_version_check=0);
  bool write_section_to_file(std::fstream &file, std::vector<my_type> limits);

  bool fft_me(data_array * data_out);
  void copy_ids( data_array * src);
  virtual bool resize(int dim, int sz);

};



#endif /* defined(____my_array__) */
