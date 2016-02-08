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

#define my_type float /**< Input data type*/
const long MAX_SIZE = 100000;/**< Maximum array size allowed (per processor if MPI in use) */

class my_array{

protected:

public:

  int n_dims;/**< Number of dimensions*/
  long * dims;/**< Array dimensions*/
  my_type *data;/**< The data */
  bool defined; /**< Flag to check memory allocation sucess*/


  my_array();
  void construct();
  my_array(long nx, long ny);
  ~my_array();

  bool is_good(){return defined;}/**< Check memory allocation etc worked*/
  int get_index(int nx);
  long get_index(long nx, long ny);
  int get_total_elements();

  int get_dims();
  int get_dims(int dim);
  int get_length(int dim);
  my_type get_element(int nx, int ny);
  bool set_element(long nx, long ny, my_type val);


};
#endif /* defined(____my_array__) */
