


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
CFLAGS = -O0 -c $(INCLUDE) -DVERSION=\"$(GIT_VERSION)\" -std=c++11 -pedantic
CFLAGS += -g
DEBUG = -g -W -Wall -pedantic -D_GLIBCXX_DEBUG -Wextra
PROFILE = -g
LFLAGS = -g
DEPSFLAGS = -DRUN_TESTS_AND_EXIT
#These should be always on for dependency generation as we want allllll the code considered.

#DEBUG+= -Wno-sign-compare
#DEBUG+= -Wno-unused-parameter
#Comment/uncomment these to hide specific errors...

SED_STR = sed -i.bak 's/ _USE_FLOAT/ _NO_USE_FLOAT/' Doxyfile
SED_STR_Test = sed -i.bak 's/ RUN_TESTS_AND_EXIT/ _NORUN_TESTS_AND_EXIT/' Doxyfile

ifeq ($(strip $(TYPE)),double)
  LIB += -lfftw3 -lm
else
  LIB += -lfftw3f -lm
  CFLAGS += -D_USE_FLOAT
  SED_STR = sed -i.bak 's/_NO_USE_FLOAT/_USE_FLOAT/' Doxyfile
#Set Doxygen to document correct version
endif

ifeq ($(strip $(MODE)),debug)
  CFLAGS += $(DEBUG)
endif

ifeq ($(strip $(MODE)),test)
  CFLAGS += -DRUN_TESTS_AND_EXIT
  CFLAGS += $(PROFILE)
  SED_STR_Test = sed -i.bak 's/ NO_RUN_TESTS_AND_EXIT/ RUN_TESTS_AND_EXIT/' Doxyfile

endif

ifeq ($(strip $(MODE)),profile)
  CFLAGS += $(PROFILE)
  #LFLAGS += $(PROFILE)
endif

#list of all header and cpp pairs. 
INCLS = my_array.h d_coeff.h spectrum.h  plasma.h tests.h reader.h controller.h non_thermal.h main_support.h

#make lists of source and object files, all headers plus main
SOURCE := $(INCLS:.h=.cpp)
OBJS := $(SOURCE:.cpp=.o)
MAINOBJS := main.o main_growth.o
#multiple mains, only onne can be included in linker calls
#header files only (no .cpp)
INCLS += support.h

#add directory prefixes
SOURCE := $(addprefix $(SRCDIR)/, $(SOURCE))
OBJS := $(addprefix $(OBJDIR)/, $(OBJS))
MAINOBJS := $(addprefix $(OBJDIR)/, $(MAINOBJS))
INCLS := $(addprefix include/, $(INCLS))

WARN_STR = "**************Run make clean before changing MODE or TYPE. Run echo_deps if \#include's change or new files added***********"
#reminder about -D options

#Yes we have multiple "main"s we can compile. They aren't in OBJS but are in SOURCE (for the dependency gen, and tarball). Deal with it and add to these compile lines explicitly
main : echo_warning $(OBJDIR)/main.o $(OBJS)
	$(CC) $(LFLAGS) $(INCLUDE) $(OBJS) $(OBJDIR)/main.o $(LIB) -o main
	@echo $(WARN_STR)
	@$(SED_STR)
	@$(SED_STR_Test)
echo_warning:
	@echo $(WARN_STR)

read_test : $(OBJDIR)/read_test.o
	$(CC) $(INCLUDE) $(SRCDIR)/read_test.cpp $(LIBSDF) -o read_test

growth : $(OBJDIR)/main_growth.o $(OBJS)
	$(CC) $(LFLAGS) $(INCLUDE) $(OBJS) $(OBJDIR)/main_growth.o $(LIB) -o growth

#testing makefile commands ;)
debug :
	@echo $(SOURCE)
	@echo " "
	@echo $(OBJS)
	@echo " "
	@echo $(MAINOBJS)
	@echo " "
	@echo $(INCLS)
	@echo " "
	@echo $(CFLAGS)
	@echo " "
	@echo $(LIB)

-include dependencies.log

echo_deps :
	@echo "Regenerating dependencies..."
	@touch dependencies.log
	@rm dependencies.log
 #touch so must exist before rm
	@for var in $(SOURCE); do $(CC) $(INCLUDE) $(DEPSFLAGS) -MM $$var |fmt -1 >> dependencies.log 2>&1;\
    done
	cp dependencies.log dependencies.log.bak
  #post processing to fix up lines and remove irrelevant deps
	@./process_deps.sh
  #prepend OBJDIR string
	@sed -i .bak 's,[a-z/_]*\.o,$(OBJDIR)\/&,' dependencies.log

 #-M dumps dependencies to file; we extract those we don't care about, such as stdio, boost and SDF libs. Each src file is appended together. To do this, first split each onto a new line, then remove the lines we don't want, append the object dir to the targets and add blank line

#Create the object directory before it is used, no error if not exists (order-only prereqs, will not rebuild objects if directory timestamp changes)
$(OBJS): | $(OBJDIR)
$(OBJDIR):
	@mkdir -p $(OBJDIR)

$(OBJDIR)/%.o:./$(SRCDIR)/%.cpp
	$(CC) $(CFLAGS)  $< -o $@

$(OBJDIR)/main.o:./$(SRCDIR)/main.cpp
	$(CC) $(CFLAGS)  $< -o $@


$(OBJDIR)/read_test.o:./$(SRCDIR)/read_test.cpp
	$(CC) $(CFLAGS)  $< -o $@

.PHONY : tar tartest clean veryclean docs

tar:
	tar -cvzf Source.tgz $(SOURCE) $(INCLS) ./files/* Makefile redox.sh process_deps.sh

clean:
	@rm main $(OBJS) $(MAINOBJS)

veryclean:
	@rm main dependencies.log*
	@rm -r $(OBJDIR)

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

