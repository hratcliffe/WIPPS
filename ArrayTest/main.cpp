

#include "stripped_array.h"
#include <stdio.h>
#include <cstdlib>
#include <iostream>


int main(int argc, char *argv[]){

  my_array *testy = new my_array(20000, 20000); // 
  
  if(!testy->is_good()) printf("%s", "Cannot allocate memory");

  for(long i=0; i<20000;i++){
    for(long j=0; j<20000; j++){
      testy->set_element(j, i, 1.0);
    
    }  
  }
  std::cout<<testy->get_element(10, 10)<<std::endl;
  
  
  /*my_type * datas = (my_type *) malloc(20000*20000*sizeof(my_type));
  
  long index ;

  for(long i=0; i<20000;i++){
    for(long j=0; j<20000; j++){
      index = (((j < 20000) && (i<20000)) ?  i*20000 + j : -1);
      if(index >= 0) datas[index] = 1.0;
    }
  }
  std::cout<<datas[100]<<std::endl;
  */
  return 0;
}

