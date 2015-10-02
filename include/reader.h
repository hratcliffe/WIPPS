//
//  reader.h
//  
//
//  Created by Heather Ratcliffe on 02/10/2015.
//
//

#ifndef _reader_h
#define _reader_h

#include <iostream>


class data_array;

class reader{

int n_z;
//Number of characters in filenames
public:

std::string file_prefix;

int space_range[2];
//This will work in x only, because my stuff only needs blocking in x...
int time_range[2];

char block_id[10];

reader(std::string file_prefix_in,  char * block_id_in);
~reader(){;}

bool read_data(data_array * my_data_in, int time_range[2], int space_range[2]);

};

#endif
