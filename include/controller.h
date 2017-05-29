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

/** 
\ingroup halp
*\defgroup bounce Bounce-averaging helpers
*@{ 
\brief Helpers for bounce-averaging process
*
* Functions to solve for mirror latitude, bounce period etc and various helpers for these. We assume the usual dipole magnetic field.
*/

/**\brief Specifies "parameters" of D to determine how to bounce average*/
typedef enum bounce_av_type_specs {plain, alpha_alpha, alpha_p, p_alpha, p_p} bounce_av_type;
/**\brief Bounce-averaging data
*
*Holds the needed stuff for bounce-averaging. Controller should only access B_x via the accessor, which returns the value from the B array. NB NB we don't currently use this B_x, rather an assumed dipole. \ext Actually use this B_x rather than assuming a dipole */
class bounce_av_data{
  private:
    data_array Bx {data_array(1)};/**< Array holding the ambient B_x field. Note this MUST be 1-D, see accessor functions get_Bx_at and set_Bx*/
    size_t len_x {0};/**< Length of the B field array*/
  public:
    bounce_av_type type {plain};/**< The specification for how to bounce average, a value from bounce_av_type */
    my_type L_shell {4.0};/**< The L shell we're at */
    my_type max_latitude {90.0};/**< Maximum latitude of "field line" in degrees*/

    /** \brief Set size of Bx
    *
    *Resizes the Bx array and updates the stored length
    @param len The new length */
    void set_Bx_size(size_t len){this->Bx.resize(0, len);len_x = Bx.get_dims(0);}

    /** \brief Set value of Bx
    *
    *Sets the Bx array to (a copy of) Bx_in
    @param Bx_in Input array, must match current size of Bx
    @return 0 for success, 1 if problem */
    bool set_Bx(data_array Bx_in){if(Bx_in.is_good() && Bx_in.get_dims() == 1 && Bx_in.get_dims(0) == len_x){Bx = Bx_in; return true;}else{return false;}}
  
    /** \brief Get Bx value
    *
    * Get the value of Bx
    @param lat The latitude to get value at in radians
    @return The value at given latitude */
    inline my_type get_Bx_at(my_type lat){return Bx.get_element(Bx.get_axis_index_from_value(0, lat));}
};

my_type solve_mirror_latitude(my_type alpha_eq, bool print_iters=false);

/** Polynomial describing the mirror latitude
@param L cos^2 lambda_mirror
@param s4alpha sin^4 alpha for alpha particle pitch angle
@return Value of mirror polynomial
*/
inline my_type mirror_poly(my_type L, my_type s4alpha){
  return std::pow(L, 6) +  (3.0*L - 4.0)*s4alpha;
}

/** Derivative of mirror_poly
@param L cos^2 lambda_mirror
@param s4alpha sin^4 alpha for alpha particle pitch angle
@return Value of the derivative of the mirror polynomial
*/
inline my_type d_mirror_poly(my_type L, my_type s4alpha){
  return 6.0*std::pow(L, 5) +  3.0*s4alpha;
}

/**Iterate solution of mirror polynomial using Newton-Raphson
@param last_guess Previous guess for N-R iteration
@param s4alpha sin^4 alpha for alpha particle pitch angle
@return Next guess for mirror poly root
*/
inline my_type Newton_Raphson_iteration(my_type last_guess, my_type s4alpha){
  return last_guess - mirror_poly(last_guess, s4alpha)/d_mirror_poly(last_guess, s4alpha);
}

/**Bounce time alpha factor from Summers Et Al \cite SummersEtAl2007 Eq 29 or Glauert and Horne \cite GlauertHorne2005 Eq 27
@param alpha_eq Equatorial particle pitch angle
@return The alpha factor
*/
inline my_type bounce_period_approx(my_type alpha_eq){
  return 1.30 - 0.56 * sin(alpha_eq);
}

/** Calculate f as in Summers Et Al \cite SummersEtAl2007 Eq 20.
@param lat Latitude in RADIANS
@return Value of f(lat)
*/
inline my_type f_latitude(my_type lat){
  return std::sqrt(1.0 + 3.0*std::pow(sin(lat), 2))/std::pow(cos(lat), 6);
}

/** Calculate alpha at given latitude from alpha_eq Summers Et Al\cite SummersEtAl2007 Eq 22.
@param alpha_eq Equatorial pitch angle in radians
@param lat Latitude in radians
@return Pitch angle at this latitude
*/
inline my_type alpha_from_alpha_eq(my_type alpha_eq, my_type lat){
  return asin(sin(alpha_eq)*std::sqrt(f_latitude(lat)));
}
/**@}*/
/**@}*/

typedef std::pair<spectrum*, diffusion_coeff*> spect_D_pair;/**< Link spectrum and D_coeff as pair*/

/** \brief Controls plasma, spectrum and d_coeff objects and their connections
*
*This is the public facing class controlling plasma, spectrum and d_coeff objects. It makes sure there is a plasma to provide needed functions for the latters and keeps each D_coeff attached to the spectrum used to generate it. Should also be responsible for supplying D's in order to whatever does the bounce-averaging.
*Because a spectrum is meaningless without a plasma, and  a diffusion coefficient meaningless without both a plasma and a spectrum, the controller class is the only thing allowed to create or destroy spectra and diffusion coefficients. Plasma's have other purposes so are not restricted in this way. When a new spectrum is created, the get_current_spectrum is set to refer to it. Any subsequent add_d operation will update the D linked to this spectrum. Note: currently controllers aren't fully implemented as objects, so work with them as pointers.
* @author Heather Ratcliffe @date 19/11/2015 \ingroup cls
\ext Consider adding copy and move constructors etc
*/
class controller{
  plasma my_plas; /**< Plasma object*/
  std::vector<spect_D_pair> spect_D_list;/**< pairs of spectrum and corresponding D*/
  std::vector<diffusion_coeff *> d_specials_list;/**< Special D's for bounce averaging etc*/
  size_t current_pair;/**<Index of current spectrum-D pair. Usually this is the last one added*/
  size_t current_special_d;/**<Index of current special D object*/

/********Plasma, spectrum, D getters ****/
  void get_D_size(size_t dims[2]);

public:

/********Basic setup functions ****/
  explicit controller(std::string file_prefix);
  /* The following disallows copying and moving, because we need a proper deep copy of the contents of spect_D_list and d_specials_list to do so, and right now there's no need*/
  controller(const controller &src) = delete;
  controller(const controller &&src) = delete;
  void clear_all();
  ~controller();
  bool is_good(){return my_plas.is_good();}/**< Whether controller is fully setup @return Boolean true for good state, false else*/

/********Plasma, spectrum, D setup functions ****/
  void set_plasma_B0(my_type Bx_ref){my_plas.set_B0(Bx_ref);}/**<Set the reference B field used by plasma, and thus the local om_ce value @param Bx_ref The value to set*/
  bool add_spectrum(std::string file);
  bool add_spectrum(int n_om, int n_ang,bool separable);
  bool add_d(int n_v, int n_angs);
  void add_d_special(int n_v, int n_angs);
  void delete_current_spectrum();
/********Plasma, spectrum, D getters ****/
  spectrum * get_current_spectrum();
  spectrum * get_spectrum_by_num(size_t indx);
  diffusion_coeff * get_current_d();
  diffusion_coeff * get_d_by_num(size_t indx);
  diffusion_coeff * get_special_d();
  
  const plasma& get_plasma(){return my_plas;}/**<Get reference to the plasma object to use @return Reference to the plasma object*/

/********Bounce averaging specials ****/
  void bounce_average(bounce_av_data bounce_dat);
  void handle_d_mpi();

/********File IO functions ****/
  bool save_spectra(std::string pref);
  bool save_D(std::string pref);
};


#endif
