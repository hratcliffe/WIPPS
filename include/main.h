//
//  main.h
//  
//
//  Created by Heather Ratcliffe on 18/09/2015.
//
//

#ifndef ____main__
#define ____main__

//#include <stdio.h>


#endif /* defined(____main__) */

//The following is Doxygen mainpage docs
/** \mainpage notitle
 
 * \section Introduction
 * This program post-processes EPOCH PIC code output files to explore Whistler wave-particle interactions. It can create wave spectra, and calculate theory based growth rates and particle diffusion coefficients.
 *
 * Function is divided across a set of utilities. `make list_utils` will list their names and each has an entry in \ref utils. For example, generate_ffts processes input data and outputs trimmed FFTs, or calculate_growth calculates theoretical growth rates of whistlers.
 *
 * The \htmlonly <a href="modules.html">Modules</a> tab \endhtmlonly \latexonly Modules section \endlatexonly shows the various parts of the code. Broadly, we have classes to store data, access data from files, represent physical entities such as plasma, and non-thermal electrons, and perform required calculations, and collections of helpers for common actions.
 *
 * \subsection setup Code Setup
 * Setup and use of the code is via a combination of config files and command line arguments. These specify the domain, working directory etc. For each program, help on the command line arguments is available with  `name -h`
 * Larger pieces of configuration are read from files. We use the deck.status file generated by EPOCH for user defined constants, a plasma.conf file for background plasma configuration and a non-thermal.conf file for non-thermal electron distributions (required for growth rate calculations). Examples of all the conf files are in the files subdirectory.
 *
 * Some of the programs can use multiple cores, by parallelising over space.
 *
 * \subsection Data Arrays
 * Data arrays are a specialised class containing data, axes and information. Get/set_element methods are provided for accessing specific parts. Spectrum and diffusion_coeff's are specialised data arrays representing physical entities.
 * \subsection FileIO SDF IO
 * SDF file reading uses the SDF C libraries. A copy of these is provided with the code, and an installer script install given to build them.
 * \subsection FFT Fourier Transforms and Special Functions
 * FFTs are handled by the FFTW routines in suitable precision (float or double). Special functions are provided by Boost.
 * \subsection contr Spectrum and D generation
 * Some objects are deeply connected. In particular, spectrums are connected to plasmas, and diffusion_coeffs are deeply connected to spectrums. Keeping these connected is the job on the controller class, which creates and stores them.

 *\section prereqs Prerequisites
 * As well as the code, an install of boost is needed (exists on OSX and most Linux systems) as well as the correct version of FFTW libraries (float for float data, double for double). A copy of EPOCH's SDF file library is included with the code. To generate the docs Doxygen and pdftex are used. See \ref build for information on how to build prereqs.

 * \section docs This Documentation
 * These docs describe all classes and methods under the \htmlonly <a href="annotated.html">Classes</a> \endhtmlonly \latexonly Classes \endlatexonly section. Helper functions, constants etc are grouped under \htmlonly <a href="modules.html">Modules</a>\endhtmlonly \latexonly Modules \endlatexonly. Assumptions and caveats within the code are collected under \ref caveats.
 *
 * Full and User docs are available. The former includes all private entities, function references, and a full source code listing. The latter does not. Change between these modes with `make DOCS=full` and `make DOCS=user`

 * \section build Building the Code
 * An install script is provided to build SDF, install fftw etc. For details run `./install --help`
 * Build using make. If input data is type double, use `make TYPE=double` to build with correct FFT etc libraries.
*
* Once the code has been built, `make tar_built` to produce a tarball Runnable.tgz of utilities and necessary files that can be copied elsewhere and run.
*\section idl IDL routines
 *Some IDL helpers are provided for reading the output files, reading the deck.status file etc. The simple way to use these is to copy the enclosed .idlstartup file to the directory you run analysis from and set this as a startup file like, for example, `pref_set, "IDL_STARTUP", "/path/to/.idlstartup",/commit` Then set a `WIPPS_PATH` environment variable giving the location of the wipps code. Start idl, and you should see a series of "% Compiled module:" lines. Try `IDL> v_to_kev(0.5)
       79.051798621061963` to test

 * \section mod Modifying the Code
  *Any changes to code include files or addition of files will change dependencies. In this case run `make echo_deps` before a clean build to regenerate the makefile listing of included headers.
 *
 * test, profile and debug modes are also available, and are invoked with `make MODE=[test,profile,debug]` after a clean. More details on testing are below.
 *
 * See the \ref todo and \ref exts pages for tasks and extensions.
 \subsection Versioning
 *Version numbers follow the usual major-minor numbering. So in general any change to file-io, data normalisation etc which may break compatibility should bump the major number, and any substantive change to how things are done bumps the minor. For example, adding a new field to data files, or changing the FFT normalisation by 2!pi is a major change. Tweaking spectrum extraction to be different but equally correct is a minor. Internal-only changes change neither number. 
 *
 *Version number is derived using a git tag, so to set one do something like `git tag -a vx.y -m "Messsage here"` before commit and make sure to push tags
 *
 *Generally only some utilities will be broken by changes, so we don't necessarily quit on mismatch. It's probably easiest to write a small file-converter to update old files and solve errors that way. 
 *\section test Integrated Testing
 * All significant parts of the code should be covered by inbuilt tests. These are defined in tests.cpp, tests_basic_and_code.cpp and tests_data_and_calc.cpp and cover a mixture of unit testing, science testing and library integration tests. To run the tests, clean build with `make clean && make MODE=test` and run `./main`. Errors and outcomes are written to stderr, information to stdout. Thus `./main 1> /dev/null` will show only the former, etc. A logfile is produced, by default named tests.log. Temporary testing files are written into tests::tests_tmp_dir.
* Some tests are rather heavy so are only run if flags are set. See \ref tests in the Test Docs for these. 
* By default, all compiles are at optimisation O3. To compile with O0 for testing or debugging, pass NO_OPT = 1 to the make command.

*Consider adding tests for any significant additions or changes, and running the existing ones in this case. To include tests in this documentation, run
 `make MODE=test docs`. Similarly, to omit them, make without test mode.
 
\section cpp Minimally wrong description of C++ as used here
*\subsection types Types and typedefs
*There's two ways to do a character string, old C-style array or characters, or a std::string (see below for meaning of std::). Some older functions expect the former, so the .c_str() conversion is used. I use C-strings for some things to match the libraries being used.
*The special type size_t is defined as an unsigned integer type large enough to count "anything". Unsigned so cannot be negative, but is used for counts, sizes etc.
* I have also added a typedef for the type in the data files and that to do the calculations in, my_type and calc_type respectively, defined in support.h. I decided to do all calculations as double, but as my data are float I use float versions of the FFTW libraries. 
*\subsection class Classes
*For this code, classes are basically structs containing data, with special methods (member functions). These always know the contents of the class and can access bits that might be hidden from the outside (private or protected).
* One class can extend another, as data_array does to my_array. The former holds just a chunk of data, the latter adds axes and additional functions.
*Once you have a class instance, i.e. a variable containing a thing of that class, you can call methods on it using the '.'
*Constructors are functions used to set up a new instance. For instance, if I want to make a 10x10 array, we set the parameters recording the dimension to define a 2-d array with sizes 10 and 10, and we grab some memory to store 10x10=100 data values etc. These special functions look like class_name(parameter list). We also have a destructor, which clean up when the variable is destroyed, and some special things which let us make copies, set one thing equal another and so on. These can be safely ignored.
*Quick example:
\code data_array dat = data_array(10, 10); //Make a new 10x10 array
 dat.set_element(5, 5, 2.0); //Set element 5, 5 to 2.0
 my_type element = dat.get_element(5, 5); //element = 2.0 \endcode
*
*\subsection default Default arguments
* Function definitions can set default values for arguments which will be used if nothing is explicitly given. For example, 
\code int add(int a, int b = 1){return a+b;}
cout<< add(1, 2)<<endl; //Prints 3
cout<< add(5)<<endl; //Prints 6
\endcode
With two arguments this prints their sum, but if only one if given, b is set to a value of 1.
*
*\subsection pointer Pointers * and ->
*Some classes get rather large, for instance if they hold a lot of data internally. In this case, you might want to pass them about not by value (copying everything) but just by getting a pointer to where they are. Pointer variables are defined like
class * my_pointer with an asterix. This variable holds only the address. Conversions between pointer and instance are:
\code class * my_pointer = new class();
class my_instance = class();
my_instance = *(my_pointer); 'dereference' pointer to get value it points to
my_pointer = &(my_instance); take 'address of' instance to get pointer. \endcode
The special operator '->' is used to apply a method to a pointer:
\code my_instance.set_element(5, 5, 2.0);
my_pointer->set_element(5, 5, 2.0);\endcode
Some of the core code uses these, none of the stuff in use should.
*The special id 'this' inside a class method refers to the instance the method was called on and is a pointer, so uses this->
*& in a function means that you may pass an instance, but a copy will not be made, the function will just be given the address. For instance the function to read data into an array dat is \code my_reader.read_data(dat, time, space);\endcode But a copy of dat is not made.

*\subsection colon The double colons ::
* The :: appears either with something like std:: or boost::math:: or with a class name, and means that this refers to the function X in that library, class etc (i.e. the function in the given namespace). So there might be a function abs() in the standard std library and in a math library and one must distinguish between them. Or both my_array and data_array have a function called is_good() and the definitions must state which they refer to.

\defgroup cls Main Classes
*@{
\brief Major classes
*
*The major classes we use, dealing with data, physics etc
@}
*/

