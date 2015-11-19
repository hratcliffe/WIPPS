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
class spectrum;
class plasma;

class diffusion_coeff: public data_array{

int n_thetas;
int n_n;


public:

int wave_id;
//ID for which wave cutout we're going for...

//calc_type * D_raw;
calc_type * D_bounceav;

diffusion_coeff(int nx, int n_ang);
virtual ~diffusion_coeff(){;}

void set_ids(float time1, float time2, int space1, int space2, int wave_id, char block_id[10]);

bool write_to_file(std::fstream &file);

void calculate(spectrum * spect, plasma * plas);

void bounce_averag(){;}

};



#endif /* defined(____d_coeff__) */
