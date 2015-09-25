//
//  spectrum.h
//  
//
//  Created by Heather Ratcliffe on 24/09/2015.
//
//

#ifndef _spectrum_h
#define _spectrum_h


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

spectrum(int nx, data_array* parent);

};


#endif
