//
//  support.h
//  
//
//  Created by Heather Ratcliffe on 24/09/2015.
//
//

#ifndef _support_h
#define _support_h


#define cplx_type fftw_complex
#define VERSION 1.0

#define my_type float
#define my_sdf_type 3

const float pi = 3.14159;
const float v0 = 2.997924e8; 

//q0 = 1.602176565d-19 ; C
//m0 = 9.10938291d-31  ; kg
//v0 = 2.99792458d8    ; m/s^2
//kb = 1.3806488d-23   ; J/K
//mu0 = 4.0d-7 * !dpi  ; N/A^2
//epsilon0 = 8.8541878176203899d-12 ; F/m
//h_planck = 6.62606957d-34 ; J s

struct deck_constants{
  float v_t;
  float omega_pe;
  float omega_ce;
  int ppc;
  

};

const int DEFAULT_N_ANG = 100;


#endif
