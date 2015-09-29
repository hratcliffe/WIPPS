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

//Shortcut so we can change array internal type later if needed. And matching SDF type

template<class T> class arrayl{

T* array_values;
T getvalue(int ix,int iy);

};

//using namespace std;
class my_array{
protected:
int n_dims;
std::vector<int> dims;
bool defined;
//flag to check memory allocation sucess

public:

my_type *data;

my_array(int nx, int ny);
~my_array();

bool is_good(){return !defined;}
my_type * get_ptr(int nx, int ny);

int get_index(int nx, int ny);
int get_dims();
int get_dims(int dim);

my_type get_element(int nx, int ny);
bool set_element(int nx, int ny, int val);

bool populate_data(my_type * dat_in, int nx, int ny);
bool populate_row(void * dat_in, int nx, int y_row);
/** Populate data

*/
bool write_to_file(std::fstream &file);
//we can use [] to wrap get elements and have the args pushed into vector which we then work with to be generic

std::string array_self_test();

};


class data_array : public my_array{
//Adds axes for data

bool ax_defined;

public:

my_type *axes;
//This will be 1-d array in sections, so can be arbitary length and dims
char block_id[10]; //the field name id form SDF file

data_array(int nx, int ny);
~data_array();

bool is_good(){return defined && ax_defined;}
my_type * get_axis(int dim, int & length);

void make_linear_axis(int dim, float res, int offset=0);
//Generates a linear axis for dimension dim, with resolution res, starting at value of offset*res

my_type * get_chunk();
//To get say a set of rows?

bool write_to_file(std::fstream &file);

};



#endif /* defined(____my_array__) */
