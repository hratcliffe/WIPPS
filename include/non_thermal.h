

//  Created by Heather Ratcliffe on 11/02/2016.


#ifndef _non_thermal_h
#define _non_thermal_h
#include<map>
/** \brief Nonthermal electron description
*
*Very small class to hold a non-thermal electron distribution we can operate on. Mainly contains routine to parse from deck.status assuming particular input deck format. We use a couple of std::functions which allows us to bind any number of constants to values and then invoke the function using single interface. That is f_p_par(p) can map to exp(-p^2/const) or equally (1+1/const*p^kappa)^1/kappa.*/

class non_thermal{

private:
  bool configure_from_file(std::string file_prefix);
  std::vector<std::function<calc_type(calc_type p_par, calc_type p_perp)> > f_p_private;/**< Bound function for f*/
  std::map<std::string, calc_type> parameters;
  calc_type norm;
  size_t ncomps;
  calc_type total_dens;
  calc_type dp;
  calc_type ref_dens;/**<Reference density (background)*/
  calc_type ref_B;/**<Reference B field*/
  calc_type fraction;/**<Non-thermal fraction*/
  calc_type v_par;/**<Parallel velocity*/
  calc_type v_perp;/**Perpendicular velocity \todo Clean up redundant params*/

public:

  non_thermal(std::string file_prefix);
  ~non_thermal(){};
  void write(std::ofstream &outfile);
  void dump(std::fstream &outfile);
  calc_type f_p(calc_type p_par, calc_type p_perp);
  calc_type d_f_p(calc_type p_par, calc_type p_perp, bool parallel);
  size_t get_n_pars(){return parameters.size();}
  void set_dp(calc_type dp){this->dp = dp;}/**<Set the dp used to get numerical derivative*/
};

calc_type bikappa(calc_type p, calc_type p2, calc_type kappa, calc_type v_k, calc_type v_k2, calc_type A);
calc_type max(calc_type p, calc_type v_th, calc_type A);
calc_type double_max(calc_type p, calc_type p_th, calc_type p_th2, calc_type A, calc_type A2);
calc_type bimax(calc_type p, calc_type p2, calc_type p_th, calc_type p_th2, calc_type A);




#endif