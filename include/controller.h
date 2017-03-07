//
//  controller.h
//  
//
//  Created by Heather Ratcliffe on 19/11/2015.
//
//

#ifndef _controller_h
#define _controller_h

#include "support.h"
#include "plasma.h"
#include "data_array.h"

class spectrum;
class diffusion_coeff;
/*Circular dependencies, don't include headers*/

typedef enum bounce_av_type_specs {plain, alpha_alpha, alpha_p,  p_p} bounce_av_type;
//Struct holding the needed stuff for bounce-averaging
//Controller should only access B_x via this accessor, because we might later allow 2-d B_x or such
class bounce_av_data{
  public:
    bounce_av_type type {plain};
    data_array Bx {data_array()};
    size_t len_x {0};
    my_type L_shell {4.0};
    my_type max_latitude {90.0};/**< Maximum latitude of "field line"*/
    bounce_av_data(){;}
    //{B_x = data_array(); len_x = 0; max_latitude = 90.0;}
    my_type get_Bx_at(size_t x_pos){return Bx.get_element(x_pos);}
};

typedef std::pair<spectrum*, diffusion_coeff*> spect_D_pair;/**< Link spectrum and D_coeff as pair*/

/** \brief Controls plasma, spectrum and d_coeff objects and their connections
*
*This is the public facing class controlling plasma, spectrum and d_coeff objects. It makes sure there is a plasma to provide needed functions for the latters and keeps each D_coeff attached to the spectrum used to generate it. Should also be responsible for supplying D's in order to whatever does the bounce-averaging.
*Because a spectrum is meaningless without a plasma, and  a diffusion coefficient meaningless without both a plasma and a spectrum, the controller class is the only thing allowed to create or destroy spectra and diffusion coefficients. Plasma's have other purposes so are not restricted in this way. When a new spectrum is created, the get_current_spectrum is set to refer to it. Any subsequent add_d operation will update the D linked to this spectrum.
* @author Heather Ratcliffe @date 19/11/2015 \ingroup cls
*/
class controller{
  plasma my_plas; /**< Plasma object*/
  std::vector<spect_D_pair> spect_D_list;/**< pairs of spectrum and corresponding D*/
  std::vector<diffusion_coeff *> d_specials_list;/**< Special D's for bounce averaging etc*/
  size_t current_pair;/**<Index of current spectrum-D pair. Usually this is the last one added*/
  size_t current_special_d;/**<Index of current special D object*/

/********Plasma, spectrum, D getters ****/
  void get_D_size(int dims[2]);

public:

/********Basic setup functions ****/
  explicit controller(std::string file_prefix);
  ~controller();
  bool is_good(){return my_plas.is_good();}/**< Whether controller is fully setup*/

/********Plasma, spectrum, D setup functions ****/
  void set_plasma_B0(my_type Bx_ref){my_plas.set_B0(Bx_ref);}/**<Set the reference B field used by plasma, and thus the local om_ce value*/
  bool add_spectrum(std::string file);
  bool add_spectrum(int nx, int n_ang,bool separable);
  bool add_d(int nx, int n_angs);
  void add_d_special(int nx, int n_angs);
  void delete_current_spectrum();
/********Plasma, spectrum, D getters ****/
  spectrum * get_current_spectrum();
  spectrum * get_spectrum_by_num(size_t indx);
  diffusion_coeff * get_current_d();
  diffusion_coeff * get_special_d();
  
  const plasma& get_plasma(){return my_plas;};/**<Get reference to the plasma object to use*/

/********Bounce averaging specials ****/
  void bounce_average(bounce_av_data bounce_dat);
  void handle_d_mpi();

/********File IO functions ****/
  bool save_spectra(std::string pref);
  bool save_D(std::string pref);
};


#endif
