#H Ratcliffe, University of Reading 2015-2017

OMPI_MPICC=clang++
CC = mpic++
SDFPATH = ./SDF/C
#Path to the SDF libraries.

#This works for me on OSX and Ubuntu. It will NOT work on windows.
#Looking for /usr/???/include/boost/ and /usr/???/lib/libboost, libfftw etc
#Default to assuming just usr/include
USR = /usr
UNAME_S := $(shell uname -s)
ifeq ($(UNAME_S),Darwin)
  USR = /usr/local
endif
ifeq ($(UNAME_S),Linux)
  USR = /usr
endif

GIT_VERSION := $(shell git describe --always --tags | cut -c1-14)
# This cleverness encodes the git commit version into the source so we can write version number into our data files

SRCDIR = src
OBJDIR = obj
INCLUDE = -I $(USR)/include/ -I $(SDFPATH)/include/ -I ./include/
LIBSDF = $(SDFPATH)/lib/libsdfc.a
LIB := $(LIBSDF)
LIB += -L $(USR)/lib/
LIBBOOST = -lboost_filesystem -lboost_system

#========Edit these for optimisation, debug options etc===============
CFLAGS = -c $(INCLUDE) -DVERSION=\"$(GIT_VERSION)\" -std=c++0x -pedantic -DPARALLEL
CFLAGS += -g
OPTIMISE = -O3
#Optimiser level for non-debug builds (debug includes test)

#NO_FFT = 1
#Uncomment to not use FFT libraries. Note FFT routines will be unavailable if set, as will certain utilities

ifdef NO_FFT
  CFLAGS += -DNO_FFT
endif

#Path for FFTW libraries. Leave empty to use system install
FFTW_PATH = ~/FFTW_testdir/fftw-3.3.4/

DEBUG = -O0 -g -W -Wall -pedantic -D_GLIBCXX_DEBUG -Wextra
#Comment/uncomment these to hide specific errors...
#DEBUG+= -Wno-sign-compare -Wno-unused-parameter
DEBUG+= -DDEBUG_ALL

#Flags just for modern gcc which is very pedantic
#DEBUG += -Wnounused-but-set-variable

PROFILE = -g
LFLAGS = -g
#Flags for profile or debug MODEs
#=====================================================================

#======================Add files here=================================

#The main program
MAINSOURCE := main.cpp
#Alternative main's for utility programs
UTILSSOURCE := cutout_fft.cpp fft_to_spectrum.cpp compress_distributions.cpp calculate_growth.cpp calculate_diffusion.cpp extract_field.cpp
ifndef NO_FFT
 UTILSSOURCE += generate_ffts.cpp
endif
#All of these are included in tarball. The former is built as default target. The latter are _ALL_ built under utils target. E.g. cutout.cpp -> cutout utility

SOURCE = my_array.cpp data_array.cpp d_coeff.cpp spectrum.cpp plasma.cpp reader.cpp controller.cpp non_thermal.cpp support.cpp resonance_poly.cpp tests.cpp tests_basic_and_code.cpp tests_data_and_calc.cpp
#list of all other cpp files. These are assumed to have both a cpp and h pairing
EGSOURCE = example_singlecore.cpp example_multicore.cpp
#Example main programs

INCLS := $(SOURCE:.cpp=.h)
OBJS := $(SOURCE:.cpp=.o)
#make lists of include and object files from SOURCE list

#INCLS += support.h
#Add files which are header only (no .cpp)
#SOURCE += blank.cpp
#Uncomment to add a source-only (no .h) file
#=====================================================================

##################Don't need to edit below here!######################
INVOKEDFILE := $(lastword $(MAKEFILE_LIST))

SCRPSDIR := ./files/scripts/

DEPSFLAGS = -DRUN_TESTS_AND_EXIT
#Dependency generation happens only once, so we want to include any code hidden in the RUN_TESTS... defined regions

#Add the FFTW_PATH to include path IFF it is non-empty
ifneq ($(strip $(FFTW_PATH)),)
  INCLUDE += -I $(FFTW_PATH)
endif
#$(info $(INCLUDE))

#---------------------------------------------------------------------
#The below sets the correct type for the FFTW etc libraries
#and also deals with make options, MODE, TEST and DOCS
#Note the spaces in search strings. DO NOT REMOVE
SED_STR = sed -i.bak 's/ USE_FLOAT/ NO_USE_FLOAT/' Doxyfile
SED_STR_Test = sed -i.bak 's/ RUN_TESTS_AND_EXIT/ NO_RUN_TESTS_AND_EXIT/' Doxyfile
#Sed strings to put the options chosen here into our Doxyfile so we document the right version

ifeq ($(strip $(TYPE)),float)
  ISFLOAT = 'Using float'
endif
ifndef TYPE
  ISFLOAT = 'Using float'
endif
ifeq ($(strip $(TYPE)),double)
  ifndef NO_FFT
    LIB += -lfftw3 -lm
  else
    LIB += -lm
  endif
else ifdef ISFLOAT
  ifndef NO_FFT
    LIB += -lfftw3f -lm
  else
    LIB += -lm
  endif
  CFLAGS += -DUSE_FLOAT
  SED_STR = sed -i.bak 's/ NO_USE_FLOAT/ USE_FLOAT/' Doxyfile
#Set Doxygen to document correct version
else ifdef TYPE
  $(error Unknown TYPE)
endif
#Check we have a valid TYPE selection and set the right flags and docs
ifeq ($(strip $(MODE)),debug)
  CFLAGS += $(DEBUG)
  #CFLAGS += -DDEBUG_DIMS
else ifeq ($(strip $(MODE)),test)
  CFLAGS += -DRUN_TESTS_AND_EXIT
  #CFLAGS += $(PROFILE)
  CFLAGS += $(DEBUG)
  ifeq ($(strip $(OPTIMISE)), 1)
    CFLAGS += $(OPTIMISE)
  endif
  #Static link boost filesystem
  LIB += $(LIBBOOST)
  SED_STR_Test = sed -i.bak 's/ NO_RUN_TESTS_AND_EXIT/ RUN_TESTS_AND_EXIT/' Doxyfile
else ifeq ($(strip $(MODE)),profile)
  CFLAGS += $(PROFILE)
  ifeq ($(strip $(OPTIMISE)), 1)
    CFLAGS += $(OPTIMISE)
  endif
  #LFLAGS += $(PROFILE)
else ifdef MODE
  $(error Unknown MODE)
else
  ifeq ($(strip $(OPTIMISE)), 1)
    CFLAGS += $(OPTIMISE)
  endif
endif
#Check we have a valid MODE selection and set the right flags and docs

#Linker flag for SDF libraries
LIB += -ldl

#Check we have a valid DOCS selection and edit Doxyfile accordingly
DOCS_SCRIPT := $(SCRPSDIR)"set_docs.sh"
ifeq ($(strip $(DOCS)),full)
  DOCS_ARG="full"
else ifeq ($(strip $(DOCS)),user)
  DOCS_ARG="user"
else ifdef DOCS
  $(error Unknown DOCS)
endif
#---------------------------------------------------------------------

#add directory prefixes
SOURCE := $(addprefix $(SRCDIR)/, $(SOURCE))
OBJS := $(addprefix $(OBJDIR)/, $(OBJS))
MAINOBJS := $(MAINSOURCE:.cpp=.o)
MAINOBJS := $(addprefix $(OBJDIR)/, $(MAINOBJS))
MAINSOURCE := $(addprefix $(SRCDIR)/, $(MAINSOURCE))
UTILS :=$(UTILSSOURCE:.cpp=)
UTILSOBJS := $(UTILSSOURCE:.cpp=.o)
UTILSOBJS := $(addprefix $(OBJDIR)/, $(UTILSOBJS))
UTILSSOURCE := $(addprefix $(SRCDIR)/, $(UTILSSOURCE))
EGSOURCE := $(addprefix $(SRCDIR)/, $(EGSOURCE))

#Add the main.h header explicitly
INCLS += "main.h"
INCLS := $(addprefix include/, $(INCLS))

WARN_STR = "**************Run make clean before changing MODE or TYPE. Run echo_deps if \#include's change or new files added***********"

#=========================Main rules=================================
#The main file, main.cpp can be compiled in normal or test mode. Target below for utils also compiles a series of other main programs.
default : echo_float echo_warning $(OBJS) $(OBJDIR)/main.o $(SDFPATH)/lib/libsdfc.a
	$(CC) $(LFLAGS) $(INCLUDE) $(OBJS) $(OBJDIR)/main.o $(LIB) -o main
	@echo $(WARN_STR)

#Builds the SDF library if absent
$(SDFPATH)/lib/libsdfc.a :
	@./build_sdf.sh $(SDFPATH)

.PHONY: update_docs echo_warning echo_float echo_usr

#Updates the Doxyfile according to selected options
#If no docs updates needed DOCS_ARG is undefined and DOCS_SCRIPT does nothing
update_docs:
	@$(SED_STR)
	@$(SED_STR_Test)
	@$(DOCS_SCRIPT) $(DOCS_ARG)

#Some reporting that should happen before building
echo_warning:
	@echo $(WARN_STR)
echo_float:
	@echo $(ISFLOAT)
#This is a test of the usr dir setting
echo_usr:
	@echo $(USR)

#Build the utilities. Each <FILE>.cpp in UTILSOBJ is assumed to create a main program which will be named <FILE>
utils : $(UTILSOBJS) $(OBJS)
	@for var in $(UTILSOBJS); do name=$$(basename $$(basename $$var .o )) && echo "Building" $$name "....." && $(CC) $(LFLAGS) $(INCLUDE) $(OBJS) $$var $(LIB) -o $$name;done

example_singlecore : $(OBJDIR)/example_singlecore.o $(OBJS) FORCE
	$(CC) $(LFLAGS) $(INCLUDE) $(OBJS) $(OBJDIR)/example_singlecore.o $(LIB) -o example_singlecore

example_multicore : $(OBJDIR)/example_multicore.o $(OBJS) FORCE
	$(CC) $(LFLAGS) $(INCLUDE) $(OBJS) $(OBJDIR)/example_multicore.o $(LIB) -o example_multicore
#====================================================================

.PHONY: FORCE
FORCE:

.PHONY: debug
#testing makefile commands ;)
debug :
	@echo $(SOURCE)
	@echo " "
	@echo $(OBJS)
	@echo " "
	@echo $(MAINOBJS)
	@echo " "
	@echo $(UTILSOBJS)
	@echo " "
	@echo $(INCLS)
	@echo " "
	@echo $(CFLAGS)
	@echo " "
	@echo $(LIB)

-include dependencies.log
#Include the dependencies file we generate in echo_deps

.PHONY : echo_deps

#Print warning if make run with dependencies.log not present
dependencies.log :
	@if [ ! -e dependencies.log ]; then echo "+++++++++++++++++Run make echo_deps to get correct file dependencies+++++++++"; fi

#Do some basic autodependency generation
echo_deps : $(SCRPSDIR)process_deps.sh
	@echo "Regenerating dependencies..."
	@touch dependencies.log
	@rm dependencies.log
  #touch so must exist before rm
	@for var in $(SOURCE) $(MAINSOURCE) $(UTILSSOURCE) $(EGSOURCE); do $(CC) $(INCLUDE) $(DEPSFLAGS) -MM $$var |fmt -1 >> dependencies.log 2>&1;\
    done
  #-M dumps dependencies to file, -MM excludes system headers;
  #Recursive dependencies are also resolved under gcc
	@cp dependencies.log dependencies.log.bak
	@$(SCRPSDIR)process_deps.sh
  #post processing to fix up lines and remove irrelevant deps
	@sed -i.bak 's,[a-zA-Z/_]*\.o,$(OBJDIR)\/&,' dependencies.log
  #prepend OBJDIR string


$(OBJS): | $(OBJDIR)
$(UTILSOBJS): | $(OBJDIR)
$(OBJDIR):
	@mkdir -p $(OBJDIR)
#Create the object directory before it is used,
#no error if not exists (order-only prereqs, will not rebuild objects if directory timestamp changes)

$(OBJDIR)/%.o:./$(SRCDIR)/%.cpp
	$(CC) $(CFLAGS)  $< -o $@
#General rule to build any obj from corresponding cpp

.PHONY : tar tar_built clean veryclean docs list cleandocs list_utils

list:
	@$(MAKE) -pRrq -f $(INVOKEDFILE) : 2>/dev/null | awk -v RS= -F: '/^# File/,/^# Finished Make data base/ {if ($$1 !~ "^[#.]") {print $$1}}' | sort | egrep -v -e '^[^[:alnum:]]' -e '^$@$$' -v -e $(OBJDIR)
#List targets.
#From http://stackoverflow.com/questions/4219255/how-do-you-get-the-list-of-targets-in-a-makefile and added exclusion for objdir and .o

list_utils:
	@echo "main" $(UTILS)

#Refresh dependencies before building the tarball
tar: dependencies.log
	tar --no-recursion -cvzf Source.tgz $(SOURCE) $(INCLS) $(MAINSOURCE) $(UTILSSOURCE) ./files/* Makefile redox.sh process_deps.sh dependencies.log install ./install_fftw.sh ./build_sdf.sh Doxyfile test_pars
	tar -cvzf SDF.tgz ./SDF

#Tar up the runnable code, excluding build details and omitting test files etc. Includes IDL scripts,
tar_built: utils
	tar --no-recursion -cvzf Runnable.tgz $(UTILS) ./main ./files/help/*.txt test_pars .idlstartup ./files/IDL/* ./SDF/IDL/* ./scripts/*

tar_docs:
	tar -cvzf Docs.tgz WIPPS.html ./html/* ./latex/refman.pdf Derivations.pdf

clean:
	@rm -f main $(UTILS) $(OBJS) $(MAINOBJS) $(UTILSOBJS)

#Clean up all generated docs. Remove html and latex dirs
cleandocs:
	@rm -rf ./html ./latex
	@rm -f ./files/tests_runtime_flags.txt
	@find . -maxdepth 1 -name "Derivations.*" ! -name "*.tex" ! -name "*.pdf" -delete

#Removes all executables, the dependencies file, the entire objdir and the compiled SDF library
veryclean: cleandocs
	@rm -f main $(UTILS) dependencies.log*
	@rm -rf $(OBJDIR)
	@rm $(SDFPATH)/lib/libsdfc.a

docs: update_docs
	@echo Running Doxygen...
	@ $(SCRPSDIR)make_docs.sh
