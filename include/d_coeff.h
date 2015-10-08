//
//  d_coeff.h
//  
//
//  Created by Heather Ratcliffe on 23/09/2015.
//
//

#ifndef ____d_coeff__
#define ____d_coeff__

#include <stdio.h>


class data_array;

class diffusion_coeff: public data_array{

public:

int wave_id;
//ID for which wave cutout we're going for...

void construct();
diffusion_coeff(int nx, int n_ang);

void set_ids(float time1, float time2, int space1, int space2, int wave_id, char block_id[10], int function_type=FUNCTION_DELTA){;}

bool write_to_file(std::fstream &file){return 0;}
};



#endif /* defined(____d_coeff__) */
