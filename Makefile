

CC = mpic++
SDFPATH = ./SDF
#Path to the SDF libraries.


GIT_VERSION := $(shell git describe --dirty --always --tags)
# This cleverness encodes the git commit version into the source so we can write version number with data

SRCDIR = src
OBJDIR = obj
INCLUDE = -I /usr/local/include/ -I $(SDFPATH)/C/include/ -I ./include/
I = i
#first letter of include dir. don't ask...
LIBSDF = -L /usr/local/lib/ $(SDFPATH)/C/lib/libsdfc.a
LIB := $(LIBSDF)

CFLAGS = -g -c $(INCLUDE) -DVERSION=\"$(GIT_VERSION)\"
DEBUG = -W -Wall -pedantic -D_GLIBCXX_DEBUG
#DEBUG+= -Wno-sign-compare
#DEBUG+= -Wno-unused-parameter
#Comment/uncomment these to hide specific errors...

ifeq ($(strip $(TYPE)),double)
  LIB += -lfftw3 -lm
else
  LIB += -lfftw3f -lm
  CFLAGS += -D_USE_FLOAT
endif

ifeq ($(strip $(MODE)),debug)
  CFLAGS += $(DEBUG)
endif

ifeq ($(strip $(MODE)),test)
  CFLAGS += -DRUN_TESTS_AND_EXIT
endif

#list of all header and cpp pairs. 
INCLS = my_array.h d_coeff.h spectrum.h  plasma.h tests.h reader.h

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

main : $(OBJS)
	$(CC) $(INCLUDE) $(OBJS) $(LIB) -o main

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

#$(OBJS): dependencies.log
dependencies.log: echo_deps
-include dependencies.log


echo_deps :
	@touch dependencies.log
	@rm dependencies.log
 #touch so must exist before rm
	@for var in $(SOURCE); do $(CC) $(INCLUDE) -MM $$var &> depout00001.txt; \
  sed -e $$'s,.h /,.h \\\ \\\n /,g;s,.cpp /,.cpp \\\ \\\n /,g;s,.h $(I),.h \\\ \\\n $(I),g;s,.cpp $(I),.cpp \\\ \\\n $(I),g' depout00001.txt &>depout00002.txt; \
  sed -e '/Xcode.app/d;/boost/d;/SDF/d;/fftw3/d;/mpi/d' depout00002.txt >>dependencies.log 2>&1; \
    done
	@sed -i .bak 's/\\ /\\/' dependencies.log
	@sed -i .bak 's,[a-z/_]*\.o,$(OBJDIR)\/&,' dependencies.log
	@sed -i .bak $$'s,[a-z/_]*\.o,\\\n&,' dependencies.log

	@rm depout00001.txt depout00002.txt

 #-M dumps dependencies to file; we extract those we don't care about, such as stdio, boost and SDF libs. Each src file is appended together. To do this, first split each onto a new line, then remove the lines we don't want, append the object dir to the targets and add blank line

#Create the object directory before it is used, no error if not exists (order-only prereqs, will not rebuild objects if directory timestamp changes)
$(OBJS): | $(OBJDIR)
$(OBJDIR):
	@mkdir -p $(OBJDIR)

$(OBJDIR)/%.o:./$(SRCDIR)/%.cpp
	$(CC) $(CFLAGS)  $< -o $@

$(OBJDIR)/read_test.o:./$(SRCDIR)/read_test.cpp
	$(CC) $(CFLAGS)  $< -o $@

.PHONY : tar tartest clean veryclean

tar:
	tar -cvzf Source.tgz $(SOURCE) $(INCLS) ./files/* Makefile input.deck

clean:
	@rm main $(OBJS) dependencies.log.bak

veryclean:
	@rm main depdendencies.log
	@rm -r $(OBJDIR)


