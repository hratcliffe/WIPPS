

CC = mpic++
SDFPATH = ./SDF
#Path to the SDF libraries.


GIT_VERSION := $(shell git describe --dirty --always --tags)
# This cleverness encodes the git commit version into the source so we can write version number with data

SRCDIR = src
OBJDIR = obj
INCLUDE = -I /usr/local/include/ -I $(SDFPATH)/C/include/ -I ./include/
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

#list of all header and cpp pairs. 
INCLS = my_array.h d_coeff.h spectrum.h reader.h plasma.h tests.h

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

#Create the object directory before it is used, no error if not exists (order-only prereqs, will not rebuild objects if directory timestamp changes)
$(OBJS): | $(OBJDIR)
$(OBJDIR):
	@mkdir -p $(OBJDIR)

#Dependencies

$(OBJDIR)/my_array.o : ./$(SRCDIR)/my_array.cpp $(INCLS)
$(OBJDIR)/spectrum.o : ./$(SRCDIR)/spectrum.cpp $(INCLS)
$(OBJDIR)/d_coeff.o : ./$(SRCDIR)/d_coeff.cpp $(INCLS)
$(OBJDIR)/plasma.o : ./$(SRCDIR)/plasma.cpp $(INCLS)


$(OBJDIR)/%.o:./$(SRCDIR)/%.cpp
	$(CC) $(CFLAGS)  $< -o $@

$(OBJDIR)/read_test.o:./$(SRCDIR)/read_test.cpp
	$(CC) $(CFLAGS)  $< -o $@


.PHONY : tar tartest clean veryclean

tar:
	tar -cvzf Source.tgz $(SOURCE) $(INCLS) ./files/* Makefile input.deck

clean:
	@rm main $(OBJS)

veryclean:
	@rm main
	@rm -r $(OBJDIR)


