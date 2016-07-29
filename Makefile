#H Ratcliffe, University of Reading 2015-2016

OMPI_MPICC=clang++
CC = mpic++
SDFPATH = ./SDF
#Path to the SDF libraries.

GIT_VERSION := $(shell git describe --dirty --always --tags)
# This cleverness encodes the git commit version into the source so we can write version number with data
#-I ./matplotpp/
SRCDIR = src
OBJDIR = obj
INCLUDE = -I /usr/local/include/ -I $(SDFPATH)/C/include/ -I ./include/
LIBSDF = -L /usr/local/lib/ $(SDFPATH)/C/lib/libsdfc.a
LIB := $(LIBSDF)
#LIB += ./matplotpp/matplotpp.a -lglut
#Add the libraries for glut (openGL) and the matplot library

#========Edit these for optimisation, debug options etc===============
CFLAGS = -O0 -c $(INCLUDE) -DVERSION=\"$(GIT_VERSION)\" -std=c++11 -pedantic
CFLAGS += -g
DEBUG = -g -W -Wall -pedantic -D_GLIBCXX_DEBUG -Wextra
#DEBUG+= -Wno-sign-compare
#DEBUG+= -Wno-unused-parameter
#Comment/uncomment these to hide specific errors...
PROFILE = -g
LFLAGS = -g
#Flags for profile or debug MODEs
#=====================================================================

#======================Add files here=================================
INCLS = my_array.h d_coeff.h spectrum.h  plasma.h tests.h reader.h controller.h non_thermal.h main_support.h
#list of all files with both header and cpp pair.

SOURCE := $(INCLS:.h=.cpp)
OBJS := $(SOURCE:.cpp=.o)
#make lists of source and object files from INCLS list
MAINSOURCE := main.cpp main_growth.cpp
UTILSSOURCE := generate_ffts.cpp
#List of source files containing a main.
#Valid program contains one and only one of these!
#Add a rule to the Main rules section to build a different one
#If not a small variant, test version etc should use own Makefile
#Include here to include in tarball
INCLS += support.h
#Add files which are header only (no .cpp)
#=====================================================================

##################Don't need to edit below here!######################
INVOKEDFILE := $(lastword $(MAKEFILE_LIST))

DEPSFLAGS = -DRUN_TESTS_AND_EXIT
#These flags should be always on for dependency generation as we do this only globally

#---------------------------------------------------------------------
#The below sets the correct type for the FFTW etc libraries
#and also deals with make options like MODE=debug or MODE=test

SED_STR = sed -i.bak 's/ _USE_FLOAT/ _NO_USE_FLOAT/' Doxyfile
SED_STR_Test = sed -i.bak 's/ RUN_TESTS_AND_EXIT/ _NORUN_TESTS_AND_EXIT/' Doxyfile
#Sed strings to put the options chosen here into our Doxyfile so we document the right version

ifeq ($(strip $(TYPE)),float)
  ISFLOAT = 'Using float'
endif
ifndef TYPE
  ISFLOAT = 'Using float'
endif
ifeq ($(strip $(TYPE)),double)
  LIB += -lfftw3 -lm
else ifdef ISFLOAT
  LIB += -lfftw3f -lm
  CFLAGS += -D_USE_FLOAT
  SED_STR = sed -i.bak 's/_NO_USE_FLOAT/_USE_FLOAT/' Doxyfile
#Set Doxygen to document correct version
else ifdef TYPE
  $(error Unknown TYPE)
endif
#Check we have a valid TYPE selection and set the right flags and docs

ifeq ($(strip $(MODE)),debug)
  CFLAGS += $(DEBUG)
else ifeq ($(strip $(MODE)),test)
  CFLAGS += -DRUN_TESTS_AND_EXIT
  CFLAGS += $(PROFILE)
  SED_STR_Test = sed -i.bak 's/ NO_RUN_TESTS_AND_EXIT/ RUN_TESTS_AND_EXIT/' Doxyfile
else ifeq ($(strip $(MODE)),profile)
  CFLAGS += $(PROFILE)
  #LFLAGS += $(PROFILE)
else ifdef MODE
  $(error Unknown MODE)
endif
#Check we have a valid MODE selection and set the right flags and docs
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

#We need these for clean etc
INCLS := $(addprefix include/, $(INCLS))

WARN_STR = "**************Run make clean before changing MODE or TYPE. Run echo_deps if \#include's change or new files added***********"
#reminder about -D options

#=========================Main rules=================================
#Yes we have multiple "main"s we can compile. They aren't in OBJS but are in SOURCE (for the dependency gen, and tarball). Deal with it and add to these compile lines explicitly
main : echo_float echo_warning $(OBJS) $(OBJDIR)/main.o 
	$(CC) $(LFLAGS) $(INCLUDE) $(OBJS) $(OBJDIR)/main.o $(LIB) -o main
	@echo $(WARN_STR)
	@$(SED_STR)
	@$(SED_STR_Test)

echo_warning:
	@echo $(WARN_STR)
echo_float:
	@echo $(ISFLOAT)
#These are rules so we can ensure they happen before any dependencies are built

read_test : $(OBJDIR)/read_test.o
	$(CC) $(INCLUDE) $(SRCDIR)/read_test.cpp $(LIBSDF) -o read_test

growth : $(OBJDIR)/main_growth.o $(OBJS)
	$(CC) $(LFLAGS) $(INCLUDE) $(OBJS) $(OBJDIR)/main_growth.o $(LIB) -o growth

#UTILSOBJS := $(UTILSSOURCE:.cpp=.o)
#UTILSOBJS := $(addprefix $(OBJDIR)/, $(UTILSOBJS))
utils : $(UTILSOBJS) $(OBJS)
	@for var in $(UTILSOBJS); do name=$$(basename $$(basename $$var .o )) && echo "Building" $$name && $(CC) $(LFLAGS) $(INCLUDE) $(OBJS) $$var $(LIB) -o $$name;done
#	for var in $(UTILSOBJS); do name=$$(basename $${var%.*}) && echo $$name && $(CC) $(LFLAGS) $(INCLUDE) $(OBJS) $$var $(LIB) -o $$name;done

#====================================================================

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
#Include in the dependencies file we generate in echo_deps

.PHONY : echo_deps

dependencies.log :
	@if [ ! -e dependencies.log ]; then echo "+++++++++++++++++Run make echo_deps to get correct file dependencies+++++++++"; fi

echo_deps :
	@echo "Regenerating dependencies..."
	@touch dependencies.log
	@rm dependencies.log
  #touch so must exist before rm
	@for var in $(SOURCE); do $(CC) $(INCLUDE) $(DEPSFLAGS) -MM $$var |fmt -1 >> dependencies.log 2>&1;\
    done
  #-M dumps dependencies to file, -MM excludes system headers;
	@cp dependencies.log dependencies.log.bak
	@./process_deps.sh
  #post processing to fix up lines and remove irrelevant deps
	@sed -i .bak 's,[a-z/_]*\.o,$(OBJDIR)\/&,' dependencies.log
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

.PHONY : tar clean veryclean docs list

list:
	@$(MAKE) -pRrq -f $(INVOKEDFILE) : 2>/dev/null | awk -v RS= -F: '/^# File/,/^# Finished Make data base/ {if ($$1 !~ "^[#.]") {print $$1}}' | sort | egrep -v -e '^[^[:alnum:]]' -e '^$@$$' -v -e $(OBJDIR)
#List targets.
#From http://stackoverflow.com/questions/4219255/how-do-you-get-the-list-of-targets-in-a-makefile and added exclusion for objdir and .o

tar: dependencies.log
	tar --no-recursion -cvzf Source.tgz $(SOURCE) $(INCLS) $(MAINSOURCE) ./files/* Makefile redox.sh process_deps.sh dependencies.log

clean:
	@rm -f main $(UTILS) $(OBJS) $(MAINOBJS)

veryclean:
	@rm -f main $(UTILS) dependencies.log*
	@rm -rf $(OBJDIR)

docs:
	@echo Running Doxygen...
	@echo "Mode: " $(MODE) "Type: " $(TYPE)
	@doxygen Doxyfile &> Doxy.log
	@echo Processing Doxygen output...
	@./redox.sh
	@echo Running pdftex...
	@cd latex ; pdflatex --file-line-error --synctex=1 -interaction nonstopmode ./refman.tex &> ../docs.log; cd ..
#	@mv ./latex/docs.log .
	@echo "Docs built. See Doxy.log and docs.log for details"

