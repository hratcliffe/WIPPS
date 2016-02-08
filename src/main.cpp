/** \file main.cpp \brief Main program
*
* This should: open a sequence of SDF files using the SDF library and read E and B field data. Fourier transform it. Extract frequency/wavenumber and angular spectra (if possible). Calculate the resulting particle diffusion coefficients using Lyons 1974 a, b, Albert 2005 and such.
* Depends on the SDF file libraries, the FFTW library, and boost's math for special functions. A set of test arguments is supplied. Call using ./main `<test_pars` to use these.
  \author Heather Ratcliffe \date 18/09/2015.
*/


#include <math.h>
#include <cmath>
#include <boost/math/special_functions.hpp>
//Provides Bessel functions, erf, and many more
#include <fstream>
#include <iostream>
#include <stdio.h>
#include "sdf.h"
//SDF file libraries
#include <mpi.h>
#include <complex.h>
#include <fftw3.h>
//FFTW3 Fourier transform libs

#include "main.h"
#include "support.h"
#include "reader.h"
#include "controller.h"
#include "plasma.h"
#include "my_array.h"
#include "d_coeff.h"
#include "spectrum.h"
#include "tests.h"

deck_constants my_const;/**< Physical constants*/
mpi_info_struc mpi_info;/**< MPI data */

#ifdef RUN_TESTS_AND_EXIT
tests* test_bed;/**<Test bed for testing */
#endif
//We wrap in ifdef for nice Doxygen docs

void get_deck_constants(std::string file_prefix);
int local_MPI_setup(int argc, char *argv[]);
setup_args process_command_line(int argc, char *argv[]);
void share_consts();
void print_help();
void divide_domain(std::vector<int>, int space[2], int per_proc, int block_num);

int main(int argc, char *argv[]){
/**
*In theory, these classes and functions should be named well enough that the function here is largely clear. Remains to be seen, eh?
*
*/


  my_array *testy = new my_array(10000,10000);

  my_type * row = (my_type *) calloc(10000, sizeof(my_type));
//  for(int i=0; i<10000;i++) row[i] = (float) i;

  for(int i=0; i<10000;i++){
//    for(int j=0; j<10000; j++){
  //    testy->set_element(i, j, 1.0);
  for(int j=0; j<10000;j++) row[j] = (float)(i+ j);
    
    //}
    testy->populate_row(row, 10000, i);
  
  
  }
  
  
  return 0;



  int err;
  
  int ierr = local_MPI_setup(argc, argv);
  if(ierr){
    std::cout<< "Error initialising MPI. ABORTING!";
    return 1;
  }

  my_print(std::string("Code Version: ")+ VERSION, mpi_info.rank);
  my_print("Code is running on "+mk_str(mpi_info.n_procs)+" processing elements.", mpi_info.rank);

  MPI_Barrier(MPI_COMM_WORLD);

  setup_args cmd_line_args = process_command_line(argc, argv);
  if(mpi_info.rank == 0) get_deck_constants(cmd_line_args.file_prefix);
  share_consts();
  /** Get constants from deck and share to other procs*/

#ifdef RUN_TESTS_AND_EXIT
  my_print("Running basic tests", mpi_info.rank);

  if(mpi_info.rank == 0) get_deck_constants("./files/test");
  share_consts();

  test_bed = new tests();
  test_bed->set_verbosity(2);
  test_bed->run_tests();
  delete test_bed;

  MPI_Finalize();

  return 0;
#else

  //Actually do the code...
  my_print("Processing "+mk_str(cmd_line_args.per_proc)+" blocks per core", mpi_info.rank);

  char block_id[ID_SIZE];
  strcpy(block_id, cmd_line_args.block.c_str());

  reader * my_reader = new reader(cmd_line_args.file_prefix, block_id);

  int n_tims = std::max(cmd_line_args.time[1]-cmd_line_args.time[0], 1);

  int my_space[2];
  my_space[0] = cmd_line_args.space[0];
  my_space[1] = cmd_line_args.space[1];

  int n_dims;
  std::vector<int> dims;
  err = my_reader->read_dims(n_dims, dims);
  if(err) safe_exit();
  int space_dim = dims[0];
  
  if(n_dims !=1) return 1;
  /**for now abort if data file wrong size... \todo FIX*/

  controller * contr;
  contr = new controller(cmd_line_args.file_prefix);


  //---------------- Now we loop over blocks per proc-------
  for(int block_num = 0; block_num<cmd_line_args.per_proc; block_num++){

    divide_domain(dims, my_space, cmd_line_args.per_proc, block_num);
    space_dim = my_space[1]-my_space[0];
    
    MPI_Barrier(MPI_COMM_WORLD);
    //--------------THIS will slightly slow down some cores to match the slowest. But it makes output easier. Consider removing if many blocks

    data_array  * dat = new data_array(space_dim, n_tims);

    if(!dat->is_good()){
      my_print("Bugger, data array allocation failed. Aborting.", mpi_info.rank);
      return 0;
    }

    err = my_reader->read_data(dat, cmd_line_args.time, my_space);
    if(err == 1) safe_exit();

    if(err == 2) n_tims = dat->get_dims(1);
    //Check if we had to truncate data array...
    data_array * dat_fft = new data_array(space_dim, n_tims);

    if(!dat_fft->is_good()){
      my_print("Bugger, data array allocation failed. Aborting.", mpi_info.rank);
      return 0;
    }

    err = dat->fft_me(dat_fft);
    
    if(mpi_info.rank ==0) MPI_Reduce(MPI_IN_PLACE, &err, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);
    else MPI_Reduce(&err, NULL, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);

    my_print("FFT returned err_state " + mk_str(err), mpi_info.rank);

    int row_lengths[2];
    row_lengths[0] = space_dim;
    row_lengths[1] = DEFAULT_N_ANG;
    
    contr->add_spectrum(row_lengths, 2);
    contr->get_current_spectrum()->make_test_spectrum(cmd_line_args.time, my_space);

    //Now we have some test spectral data we can work with...

    contr->add_d(cmd_line_args.d[0], cmd_line_args.d[1]);
    contr->get_current_d()->calculate();

    delete dat;
    delete dat_fft;

  }
  //-----------------end of per_proc loop---- Now controller holds one spectrum and d per block
  MPI_Barrier(MPI_COMM_WORLD);
  
  contr->bounce_average();
  
  contr->save_spectra(cmd_line_args.file_prefix);
  contr->save_D(cmd_line_args.file_prefix);

  //Cleanup objects etc
  delete my_reader;
  delete contr;

  std::cout<<"Grep for FAKENUMBERS !!!!"<<std::endl;

  ADD_FFTW(cleanup());
  MPI_Finalize();
  //call these last...
#endif

  exit(0);
}

/** \defgroup main Main Helper Functions
*@{ */

int local_MPI_setup(int argc, char *argv[]){
/** \brief Do the MPI init
*
* Calls MPI init and sets up communicator. Stores number of processors and each rank into mpi_info.
*/
  int ierr, rank, n_procs;

  ierr = MPI_Init(&argc, &argv);
  //Note any other command line arg processing should account for this...

  ierr = MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  ierr = MPI_Comm_size(MPI_COMM_WORLD, &n_procs);
  mpi_info.rank = rank;
  mpi_info.n_procs = n_procs;
  return ierr;
}

void safe_exit(){
/** \brief Exit program
*
* Does minimal cleanup and exits
*/

  MPI_Finalize();

  exit(0);
}

void share_consts(){
/** \brief MPI Share deck constants
*
*Share the deck constants read on root to all procs
*/

  const int count = 5;
  int lens[count] ={1, 1, 1, 1, 1};
  MPI_Aint disps[count];
  MPI_Datatype types[count] = {MPI_FLOAT, MPI_FLOAT, MPI_FLOAT, MPI_FLOAT, MPI_INT};
  MPI_Datatype deckConstType;
  
  disps[0] = offsetof(deck_constants, v_t);
  disps[1] = offsetof(deck_constants, omega_pe);
  disps[2] = offsetof(deck_constants, omega_ce);
  disps[3] = offsetof(deck_constants, omega_ci);
  disps[4] = offsetof(deck_constants, ppc);
  
  MPI_Type_create_struct(count, lens, disps, types, &deckConstType);
  MPI_Type_commit(&deckConstType);

  MPI_Bcast(&my_const, 1, deckConstType, 0, MPI_COMM_WORLD);

}

setup_args process_command_line(int argc, char *argv[]){
/** \brief Set basic parameters
*
*Sets defaults or those given via command line (see help.txt). Forces an constant integer number of space blocks on each core.
*/

  setup_args values;
  values.n_space = -1;
  values.space[0] = -1;
  values.space[1] = -1;

  values.time[0] = 0;
  values.time[1] = 1;
  values.time[2] = 0;
  values.file_prefix = "./files/";
  values.block = "ex";
  values.d[0] = 10;
  values.d[1] = 10;

  for(int i=1; i< argc; i++){
    if(strcmp(argv[i], "-h")==0) print_help();
    if(strcmp(argv[i], "-f")==0 && i < argc-1) values.file_prefix = argv[i+1];
    if(strcmp(argv[i], "-start")==0 && i < argc-1) values.time[0] = atoi(argv[i+1]);
    if(strcmp(argv[i], "-end")==0 && i < argc-1) values.time[1] = atoi(argv[i+1]);
    if(strcmp(argv[i], "-rows")==0 && i < argc-1) values.time[2] = atoi(argv[i+1]);

    if(strcmp(argv[i], "-block")==0 && i < argc-1) values.block = argv[i+1];
    if(strcmp(argv[i], "-n")==0 && i < argc-1) values.n_space= atoi(argv[i+1]);
    if(strcmp(argv[i], "-space")==0 && i < argc-2){
      values.space[0] = atoi(argv[i+1]);
      values.space[1] = atoi(argv[i+2]);
    }
    if(strcmp(argv[i], "-d")==0 && i < argc-2){
      values.d[0] = atoi(argv[i+1]);
      values.d[1] = atoi(argv[i+2]);
    }
    
  }

  if(values.space[0] == -1 && values.space[1] == -1 && values.n_space == -1){
    values.n_space = mpi_info.n_procs;
    values.per_proc = 1;
  }
  //By default do one block per processor
  else{
    values.per_proc = values.n_space / mpi_info.n_procs;
    //integer division!
    values.n_space = values.per_proc * mpi_info.n_procs;
  //Exactly divide
  }

  if(values.time[0]< 0 ) values.time[0] = 0;
  if(values.time[1]< 0 ) values.time[1] = 0;
  if(values.time[1]< 0 ) values.time[1] = 0;
  if(values.time[1] < values.time[0]) values.time[1] = values.time[0] + 1;
  if(values.d[0] < 0) values.d[0] = 0;
  if(values.d[1] < 0) values.d[1] = 0;
  if(values.d[0] >MAX_SIZE){
    values.d[0] = MAX_SIZE;
    my_print("WARNING: Requested size exceeds MAXSIZE", mpi_info.rank);
  }
  if(values.d[1] >MAX_SIZE){
    values.d[1] = MAX_SIZE;
    my_print("WARNING: Requested size exceeds MAXSIZE", mpi_info.rank);
  }
  //Protect from invalid user input
  
  return values;
}

void print_help(){
/** \brief Print command line help
*
*Prints contents of halp_file from rank zero and calls safe exit.
*/
  std::ifstream halp;
  
  halp.open(halp_file);
  if(mpi_info.rank == 0){
    std::cout<<"Command line options: "<<std::endl;
    std::cout<<halp.rdbuf();
  }
  safe_exit();
}

void divide_domain(std::vector<int> dims, int space[2], int per_proc, int block_num){
/** \brief Divide dims evenly between procs
*
*Uses the number of space blocks from args (if specified) and the domain size from dims to ensure perfect subdivision and set current proc's bounds. We can ignore incoming space vals as they should be -1
*/

    int end, block_start, block_end, block_len;
    float per_proc_size;
    end = dims[0];
    per_proc_size = std::ceil( (float) end / (float) mpi_info.n_procs);
    //Force overlap not missing
    block_start = mpi_info.rank * (int) std::floor((float) end/(float) (mpi_info.n_procs));
    block_end = block_start + (int) per_proc_size;
    block_len = block_end - block_start;
    space[0] = block_start + (block_num) * block_len / per_proc;
    space[1] = space[0] + block_len / per_proc;

}

int where(my_type * ax_ptr, int len, my_type target){
/** \brief Find where ax_ptr exceeds target
*
*Checks bounds and calls locally scoped recursive whereb to do the search. Special case if target equals bottom end
*/
  int whereb(my_type * ax_ptr, int len, my_type target, int &cut,int sign); //Recursive function to do the finding

  int sign = 1.0;
  int cut = 0;
  if(ax_ptr[0] == target) return 0;
  if(ax_ptr[0] > target || ax_ptr[len-1] < target) return -1;
  return whereb(ax_ptr, len, target, cut, sign);

}

int whereb(my_type * ax_ptr, int len, my_type target,int &cut, int sign){
/**\brief Recursive binary bisection find. 
*
*First index where ax_ptr exceeds target. If sign ==-1 it should do < but as yet untested....
*/
  
  if(len==1){
     cut++;
     return cut;
  }
  else if(sign*ax_ptr[len/2] == sign*target){
    //Exact equality. Save recursing deeper
    cut+= len/2;
    return cut;
  }
  else if(sign*ax_ptr[len/2] > sign*target){
    whereb(ax_ptr, len/2, target, cut, sign);

    return cut;
  }
  else if(sign*ax_ptr[len/2] < sign*target){
    cut+= len/2;
    whereb(ax_ptr+len/2, len-len/2, target, cut, sign);
    return cut;
  }

  return -1;
  //Shouldn't ever reach this case, but.
}

void get_deck_constants(std::string file_prefix){
/** \brief Setup run specific constants
*
*Reads deck.status and parses values for user defined constants etc. It will rely on using the specific deck, because it has to look for names. Any changes to deck may need updating here. Tag names are set as const strings in support.h
*/


//open file deck.status from whatever path
//Find the line " Constant block values after"
//Read into vector up to "Deck state" or eof
//parse out the values we want given by name in support.h

  std::ifstream infile;
  infile.open(file_prefix+"deck.status");
  std::string header_row, line;
  infile>> header_row;
  std::vector<std::string> lines;
  bool found=false;
  if(!infile.is_open()){
    my_print("No deck.status file found, aborting", mpi_info.rank);
    safe_exit();
  
  }
  while(infile){
    getline(infile, line);
    if(!(strcmp(line.substr(0, CONSTANTS.size()).c_str(), CONSTANTS.c_str()))){
      found = true;
      break;
    }
  }
  if(!found){
    my_print("No constants dump found in deck.status, aborting", mpi_info.rank);
    safe_exit();
  }

  while(infile){
    getline(infile, line);
    if(!(strcmp(line.substr(0, CONSTANTS_END.size()).c_str(), CONSTANTS_END.c_str()))) break;
    lines.push_back(line);
  }
  
  std::string name, val;
  bool parse_err;
  float val_f;
  for(size_t i=0; i< lines.size(); i++){

    parse_err = parse_name_val(lines[i], name, val);
    if(parse_err) continue;
    val_f = atof(val.c_str());
  
    if(name == OMEGA_CE) my_const.omega_ce = val_f;
    else if(name == OMEGA_PE) my_const.omega_pe = val_f;
    //else if(name == DENS_RAT) my_const.dens_rat = val_f;
    //else if(name == PPC) my_const.ppc = val_f;
    
  }

//my_const.omega_ce = 17588.1111;

//my_const.omega_pe =2.0*35176.401757;

  my_const.omega_ci = my_const.omega_ce * me/mp;
  //assumes same charge magnitude, I.e. H plasma

}

/** @} */


bool parse_name_val(std::string in, std::string &name, std::string &val){
/** \brief Parse x=y strings
*
* Basic line parser. Takes a string and if it contains an '=' splits into the left and right segments, stripping leading and trailing spaces. Returns 0 if success, 1 if no equals sign. Standard comment character is # as first non-whitespace
*/
  if(in[in.find_first_not_of(" \t\n")] == '#') return 1;
  size_t pos = in.find("=");
  if(pos == std::string::npos){
    //is not a x = y line
    name = "";
    val="";
    return 1;

  }else{
    name = in.substr(0, pos);
    trim_string(name);
    val = in.substr(pos+1, in.size());
    trim_string(val);
    return 0;
  }
}

void trim_string(std::string &str, char ch){
  std::string tmp;
  if(str.find_first_not_of(ch) !=std::string::npos) tmp = str.substr(str.find_first_not_of(ch), str.size());
  str=tmp;
  if(str.find_first_not_of(ch) !=std::string::npos) tmp = str.substr(0, str.find_last_not_of(ch)+1);
  str=tmp;

}

void my_print(std::string text, int rank, int rank_to_write){
/** \brief Write output
*
* Currently dump to term. Perhaps also to log file. Accomodates MPI also. Set rank_to_write to -1 to dump from all. Default value is 0
*/
  if(rank == rank_to_write || rank_to_write == -1){
  
    std::cout<< text<<std::endl;
  }

}
void my_print(std::fstream * handle, std::string text, int rank, int rank_to_write){
/** \brief Write output
*
* Currently dump to term. Perhaps also to log file. Accomodates MPI also. Set rank_to_write to -1 to dump from all. Default value is 0
*/
  if((rank == rank_to_write || rank_to_write == -1) && handle!=nullptr){
    *handle<<text<<std::endl;
  }else if(rank == rank_to_write || rank_to_write == -1){
    std::cout<<text<<std::endl;

  }

}

std::string mk_str(int i){

  char buffer[25];
  std::sprintf(buffer, "%i", i);
  std::string ret = buffer;
  return ret;
  
}

std::string mk_str(double i, bool noexp){

  char buffer[25];
  if(noexp) std::snprintf(buffer, 25, "%f", i);
  else std::snprintf(buffer, 25, "%e", i);
  std::string ret = buffer;
  return ret;
  
}

std::string mk_str(bool b){

  if(b) return "1";
  else return "0";

}

std::string mk_str(long double i, bool noexp){return mk_str((double) i, noexp);};
std::string mk_str(float i, bool noexp){return mk_str((double) i, noexp);};

template<typename T> T integrator(T * start, int len, T * increment){
/** \brief Basic numerical integrator
*
*Uses trapezium rule. WARNING this is working with contiguous memory. Not very C++ but faster.
*/

  T value=0.0;
  
  for(int i=0; i<len-1; i++){
  
    value += 0.5*(start[i] + start[i+1]) * increment[i];
    
  }
//  value += start[len-1]*increment[len-1];
  //top bnd we assume flat

 return value;

}

template float integrator<float>(float *, int, float *);
template double integrator<double>(double *, int, double *);
//We need both float and double versions

calc_type square_integrator(calc_type * start, int len, calc_type * increment){
/** \brief Basic numerical integrator
*
*Uses trapezium rule. WARNING this is working with contiguous memory. Not very C++ but faster.
*/

  calc_type value=0.0;
  
  for(int i=0; i<len-1; i++){
  
    value += 0.5*(start[i]*start[i] + start[i+1]*start[i+1]) * increment[i];
    
  }
//  value += start[len-1]*increment[len-1];
  //top bnd we assume flat

 return value;

}

template<typename T> void inplace_boxcar_smooth(T * start, int len, int width, bool periodic){
/** \brief Boxcar smoothing of specified width
*
*Smooths the array given by start and len using specified width. If periodic is set the ends wrap around. Otherwise they one-side
*/

  if(width > len) my_print("Really? How am I meant to smooth that?", mpi_info.rank);
  calc_type result = 0.0;
  int edge = width/2;
  for(int i=0; i<width; ++i) result += start[i];
  for(int i=edge+1; i<len-edge; ++i){
    result -= start[i-edge-1];
    result += start[i+edge];
    //remove behind and add in front. Faster for width>2
    start[i] = result / (calc_type) width;
  }

  //Handle ends
  if(periodic){
    result = 0.0;
    for(int i=0; i<width-1; ++i) result += start[i];
    //first width-1
    int wrap = len -1;
    //and one wrapped back
    result +=start[wrap];
    for(int i=0; i<edge; ++i){
      result -= start[width-2 - i];
      result += start[wrap - i - 1];
      start[edge-i] = result / (calc_type) width;
    }
    for(int i=0; i<edge; ++i){
      result -= start[width-2 - i];
      result += start[wrap - edge - i - 1];
      start[len-i-1] = result / (calc_type) width;
    }


  }else{
    //we just truncate at the bottom
    
  
  }


}
template void inplace_boxcar_smooth(calc_type *, int, int, bool);

std::vector<calc_type> cubic_solve(calc_type an, calc_type bn, calc_type cn){
/** \brief Finds roots of cubic x^3 + an x^2 + bn x + cn = 0
*
* Uses Num. Rec. equations, which are optimised for precision. Note that if x >>1 precision errors may result. Returns real solutions only
*/

  calc_type Q, R, bigA, bigB, Q3, R2, bigTheta;
  std::vector<calc_type> ret_vec;

  Q = (std::pow(an, 2) - 3.0 * bn)/9.0;
  R = (2.0* std::pow(an, 3) - 9.0 * an *bn + 27.0*cn)/54.0;
  
  R2 = std::pow(R, 2);
  Q3 = std::pow(Q, 3);
  
  if( R2 < Q3){
    
    bigTheta = std::acos(R/sqrt(Q3));
    calc_type minus2sqrtQ = -2.0*std::sqrt(Q);
    
    ret_vec.push_back(minus2sqrtQ*std::cos(bigTheta/3.0) - an/3.0);
    ret_vec.push_back(minus2sqrtQ*std::cos((bigTheta + 2.0*pi)/3.0) - an/3.0);
    ret_vec.push_back(minus2sqrtQ*std::cos((bigTheta - 2.0*pi)/3.0) - an/3.0);

  }else{
    calc_type ret_root;
    bigA = - boost::math::sign(R)*std::pow((std::abs(R) + std::sqrt(R2 - Q3)), 1.0/3.0 );

    (bigA != 0.0) ? (bigB = Q / bigA) : (bigB = 0.0);
    ret_root = (bigA + bigB) - an/3.0;

    ret_vec.push_back(ret_root);
  }

/** Used to test when writing
  calc_type tmp;
  for(int i=0; i<ret_vec.size(); ++i){
    
    tmp = std::pow(ret_vec[i], 3) + an*std::pow(ret_vec[i], 2) + bn*ret_vec[i] + cn;
    std::cout<<"solution gives "<<tmp<<std::endl;
  
  }
*/
  return ret_vec;

}

template<typename T> T interpolate(T* axis, T* vals, T target, int pts){
/** Interpolate vals on axis to target value
*
*For pts=1 uses closest value, pts=2 uses 2 pt linear, \todo add more pts options
*/

  T ret = 0.0;
  if(pts ==1){
    //select closer value
    if(std::abs(target - axis[0]) <= std::abs(target - axis[1])) ret = vals[0];
    else ret = vals[1];
  
  }else if(pts ==2){
    ret = (std::abs(target - axis[1]) * vals[0] + std::abs(target - axis[0]) * vals[1])/(std::abs(axis[1] - axis[0]));
  
  }

  return ret;
}
template float interpolate(float*, float*, float, int);
template double interpolate(double*, double*, double, int);

