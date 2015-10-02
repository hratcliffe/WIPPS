//
//  reader.cpp
//  
//
//  Created by Heather Ratcliffe on 02/10/2015.
//
//

#include <stdio.h>
#include <fstream>
#include <iostream>
#include "support.h"
#include "reader.h"


/** This is going to take care of openeing and reading from SDF. We'll give it a pointer to a storage class instance, and it can then use that to obtain pointer to the bit of memory to fill. We'll also give it a numerical range, a file prefix, and a block name to grab. Perhaps we also have option to restrict on the ranges, so we can easily block in space.
*
*It will be sort of safe. It wont store the data array to write to, but will need to be fed it. So we can;t overwrite the memory.
*
*/

reader::reader(std::string file_prefix_in,  char * block_id_in){
/**We'll also have to work out how many zeros to use at some point. For nw we'll do it when we make the reader eh, as that's where we assign the prefix. We'll assume a 0th file exists, and default and minimum is 4. We'll also try 5-7.
*/
  strcpy(this->block_id, block_id_in);
  this->file_prefix = file_prefix_in;
  n_z=0;
  
  std::ifstream file;
  std::string name;
  name = file_prefix + "000.sdf";

  file.open(name);
  if(!file.is_open()){
    //wrong number, try with increasing numbers of zeros
    
    for(int i=4; i<8;i++){
      name.insert(file_prefix.size()+1, "0");
      std::cout<<name<<std::endl;

      file.open(name);
      if(file.is_open()){
        file.close();
        n_z = i;
        break;
      }
      
    }
    
  }else{
    n_z =4;
    file.close();
  }
}

bool reader::read_data(data_array * my_data_in, int time_range[2], int space_range[2]){
/** This will open the files dictated by time range seuqentially, and populate them into the data_array. It'll stop when the end of rnage is reached, or it goes beyond the size available. Space range entry of -1 is taken as respective limit.
*/


//We use sprintf to make the file name strings...

return 0;
}

