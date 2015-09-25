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


//Lets have a spectrum class then
//Contains data, axis, sizes, ids (field, time range, space range)
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

//We can't hold parent ID, as we don't know when parent might be destroyed...

spectrum(int nx);
spectrum(int nx, data_array* parent);

void set_ids(float time1, float time2, int space1, int space2, int wave_id, char block_id[10]);

bool generate_spectrum(data_array * parent);

float get_dispersion(my_type k, int wave_type);

};


#endif
