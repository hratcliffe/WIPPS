

CC = mpic++
INCLUDE = -I /usr/local/include/ -I ./SDF/C/include/
SRCDIR = src
OBJDIR = obj
CFLAGS = -g -c $(INCLUDE)
DEBUG = -W -Wall -pedantic -D_GLIBCXX_DEBUG
#DEBUG+= -Wno-sign-compare
DEBUG+= -Wno-unused-parameter
#Comment/uncomment these to hide specific errors...
LIB = -L /usr/local/lib/ ./SDF/C/lib/libsdfc.a -lfftw3 -lm

ifeq ($(strip $(MODE)),debug)
  CFLAGS += $(DEBUG)
endif
#links against ncurses libraries

#list of all header and cpp pairs. Add new files here.....
INCLS = my_array.h d_coeff.h

#make lists of source and object files, all headers plus main
SOURCE := $(INCLS:.h=.cpp)
SOURCE += main.cpp
OBJS := $(SOURCE:.cpp=.o)

#header files only (no .cpp)
#INCLS += support.h

#add directory prefixes
SOURCE := $(addprefix $(SRCDIR)/, $(SOURCE))
OBJS := $(addprefix $(OBJDIR)/, $(OBJS))
INCLS := $(addprefix include/, $(INCLS))

main : $(OBJS)
	$(CC) $(INCLUDE) $(OBJS) $(LIB) -o main

#testing makefile commands ;)
debug :
	@echo $(SOURCE)
	@echo " "
	@echo $(OBJS)
	@echo " "
	@echo $(INCLS)


#Create the object directory before it is used, no error if not exists (order-only prereqs, will not rebuild objects if directory timestamp changes)
$(OBJS): | $(OBJDIR)
$(OBJDIR):
	@mkdir -p $(OBJDIR)

#Dependencies

obj/main.o : ./src/main.cpp $(INCLS)
obj/my_array.o : ./src/my_array.cpp $(INCLS)


obj/%.o:./src/%.cpp
	$(CC) $(CFLAGS)  $< -o $@

.PHONY : tar tartest clean veryclean

tar:
	tar -cvzf Source.tgz $(SOURCE) $(INCLS) ./files/* Makefile input.deck

clean:
	@rm main $(OBJS)

veryclean:
	@rm main
	@rm -r $(OBJDIR)


