//
//  spectrum.h
//  
//
//  Created by Heather Ratcliffe on 24/09/2015.
//
//

#ifndef _spectrum_h
#define _spectrum_h

const int WAVE_WHISTLER = 1;
const int WAVE_PLASMA = 2;
const int WAVE_O = 3;

const int FUNCTION_DELTA = 1;
const int FUNCTION_GAUSS = 2;
const int FUNCTION_ISO = 3;

//Lets have a spectrum class then
//Contains data, axis, sizes, ids (field, time range, space range)
//For now, let's stick to the 1-d, then there's only delta say in angle...
//We can impose say a narrow Gaussian later.
//NB Note in 1-d also omega and k_pllel and k are freely interchangeable
//Alright. Lets just add an extra row for the B^2, and fill the rest of the array with g_omega. It's never going to be a massive amount of data
//But if the angle profile is say a function, we don't need to...

class data_array;

class spectrum : public data_array{
//specialised data array with extra ID fields to hold additional data

public:

float time[2];
//time range over which spectrum was derived
int space[2];
//space range ditto

int wave_id;
//ID for which wave cutout we're going for...

bool angle_is_function;
//Says either we're in 1-D, and we impose functional form for anguar profile, e.g delta function or Gaussian, or we're in 2-D but still wish to impose.

int function_type;

int n_angs;


//We can't hold parent ID, as we don't know when parent might be destroyed...

spectrum(int nx);
spectrum(int nx, int n_ang);
spectrum(int nx, data_array* parent);

void set_ids(float time1, float time2, int space1, int space2, int wave_id, char block_id[10], int function_type=FUNCTION_DELTA);

bool generate_spectrum(data_array * parent);

float get_dispersion(my_type k, int wave_type);

my_type * get_angle_distrib(my_type ang, my_type omega);

int where(my_type * ax_ptr, int len, my_type target);

};


#endif
