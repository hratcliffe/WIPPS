//
//  d_coeff.cpp
//  
//
/** \file d_coeff.cpp This file will calculate the diffusion coefficients from a spectrum. We'll try first just going for it. Does it make sense to have this be any sort of object? Possibly, in that we want to keep info on how, where, etc, etc, and the axes. So possibly it is again descended from a data_array. For now, implement some basic calcs
*
* @author Heather Ratcliffe @date 23/09/2015.
*/

#include <math.h>
#include <boost/math/special_functions.hpp>
#include "support.h"
#include "my_array.h"
#include "spectrum.h"
#include "d_coeff.h"

extern mpi_info_struc mpi_info;

/**First we want to generate a particle velocity axis, and a pitch angle axis. Possibly the latter matches the one from spectrum?
*We need the refractive index. 
*Probably numerical integration/differentiation routines


*/


