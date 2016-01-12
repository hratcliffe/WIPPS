


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
CFLAGS = -O0 -c $(INCLUDE) -DVERSION=\"$(GIT_VERSION)\" -std=c++11
CFLAGS += -g
DEBUG = -g -W -Wall -pedantic -D_GLIBCXX_DEBUG -Wextra
PROFILE = -g
LFLAGS = -g

#DEBUG+= -Wno-sign-compare
#DEBUG+= -Wno-unused-parameter
#Comment/uncomment these to hide specific errors...

SED_STR = sed -i.bak 's/ _USE_FLOAT/ _NO_USE_FLOAT/' Doxyfile

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
endif

ifeq ($(strip $(MODE)),profile)
  CFLAGS += $(PROFILE)
  #LFLAGS += $(PROFILE)
endif

#list of all header and cpp pairs. 
INCLS = my_array.h d_coeff.h spectrum.h  plasma.h tests.h reader.h controller.h

#make lists of source and object files, all headers plus main
SOURCE := $(INCLS:.h=.cpp)
SOURCE += main.cpp
OBJS := $(SOURCE:.cpp=.o)

#header files only (no .cpp)
INCLS += support.h

#add directory prefixes
SOURCE := $(addprefix $(SRCDIR)/, $(SOURCE))
OBJS := $(addprefix $(OBJDIR)/, $(OBJS))
INCLS := $(addprefix include/, $(INCLS))

WARN_STR = "**************Run make clean before changing MODE or TYPE. Run echo_deps if \#include's change or new files added***********"
#reminder about -D options

main : echo_warning $(OBJS)
	$(CC) $(LFLAGS) $(INCLUDE) $(OBJS) $(LIB) -o main
	@echo $(WARN_STR)
	@$(SED_STR)
echo_warning:
	@echo $(WARN_STR)

read_test : $(OBJDIR)/read_test.o
	$(CC) $(INCLUDE) $(SRCDIR)/read_test.cpp $(LIBSDF) -o read_test

#testing makefile commands ;)
debug :
	@echo $(SOURCE)
	@echo " "
	@echo $(OBJS)
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
	@for var in $(SOURCE); do $(CC) $(INCLUDE) -MM $$var |fmt -1 >> dependencies.log 2>&1;\
    done
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

$(OBJDIR)/read_test.o:./$(SRCDIR)/read_test.cpp
	$(CC) $(CFLAGS)  $< -o $@

.PHONY : tar tartest clean veryclean docs

tar:
	tar -cvzf Source.tgz $(SOURCE) $(INCLS) ./files/* Makefile input.deck

clean:
	@rm main $(OBJS)

veryclean:
	@rm main dependencies.log*
	@rm -r $(OBJDIR)

docs:
	doxygen Doxyfile &> Doxy.log
	./redox.sh
	cd latex ; pdflatex --file-line-error --synctex=1 -interaction nonstopmode ./refman.tex &> pdflatex.log
	@echo "Docs built. See Doxy.log and pdflatex.log for details"

