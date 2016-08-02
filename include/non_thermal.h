

//  Created by Heather Ratcliffe on 11/02/2016.


#ifndef _non_thermal_h
#define _non_thermal_h

/** \brief Nonthermal electron description
*
*Very small class to hold a non-thermal electron distribution we can operate on. Mainly contains routine to parse from deck.status assuming particular input deck format.*/

class non_thermal{

private:
  bool configure_from_file(std::string file_prefix);

public:

  non_thermal(std::string file_prefix);
  ~non_thermal();
  calc_type ref_dens;/**<Reference density (background)*/
  calc_type ref_B;/**<Reference B field*/
  calc_type fraction;/**<Non-thermal fraction*/
  calc_type v_par;/**<Parallel velocity*/
  calc_type v_perp;/**Perpendicular velocity*/

  void write(std::ofstream &outfile);
};





#endif