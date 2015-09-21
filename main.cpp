//
//  main.cpp
//  
//
//  Created by Heather Ratcliffe on 18/09/2015.
//
//

#include "main.h"
#include "my_array.h"
//#include <math>
#include <boost/math/special_functions.hpp>
#include <fstream>

using namespace std;

int main(int argc, char *argv[]){

//Array class from std C++11 or boost
//own multi dim version
//version with axes


//Test use and syntax of FFTW

my_array dat = my_array(5, 5);
bool err;
for(int i=0; i<5; i++){
  for(int j =0; j<5; j++){
    err=dat.set_element(i, j, i*j);
    if(err) cout<<"Uh oh!";
  }
}

for(int i=0; i<5; i++){
  for(int j =0; j<5; j++){
    cout<<dat.get_element(i,j)<<" ";
  }
  cout<<endl;
}



}
