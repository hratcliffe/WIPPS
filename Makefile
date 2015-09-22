

CC = g++
INCLUDE = -I ./includes/
SRCDIR = src
OBJDIR = obj
CFLAGS = -g -c $(INCLUDE)
DEBUG = -W -Wall -pedantic -D_GLIBCXX_DEBUG
#DEBUG+= -Wno-sign-compare
DEBUG+= -Wno-unused-parameter
#Comment/uncomment these to hide specific errors...


ifeq ($(strip $(MODE)),debug)
  CFLAGS += $(DEBUG)
endif
#links against ncurses libraries
NCURS = -lncurses

#list of all header and cpp pairs. Add new files here.....
INCLS = drawer.h page.h page_items.h sub_page.h dividablepage.h page_overlays.h config.h block.h block_page.h io.h deck_support.h block_set.h regionalisation_consts.h deck_file.h clipboard.h

#Subsets of header files to include
PAGEINCLS = page_items.h drawer.h page.h dividablepage.h block_page.h sub_page.h page_overlays.h
SUPPORTINCLS = cmd_aliases.h deck_support.h config.h regionalisation_consts.h io.h
BLOCKINCLS = block.h block_set.h deck_file.h

#make lists of source and object files, all headers plus main
SOURCE := $(INCLS:.h=.cpp)
SOURCE += decker.cpp
OBJS := $(SOURCE:.cpp=.o)

#header files only (no .cpp)
INCLS += cmd_aliases.h

#add directory prefixes
SOURCE := $(addprefix $(SRCDIR)/, $(SOURCE))
OBJS := $(addprefix $(OBJDIR)/, $(OBJS))
INCLS := $(addprefix includes/, $(INCLS))
PAGEINCLS := $(addprefix includes/, $(PAGEINCLS))
SUPPORTINCLS := $(addprefix includes/, $(SUPPORTINCLS))
BLOCKINCLS := $(addprefix includes/, $(BLOCKINCLS))

decker : $(OBJS)
	$(CC) $(OBJS) -o decker $(NCURS)

#testing makefile commands ;)
debug :
	@echo $(SOURCE)
	@echo " "
	@echo $(OBJS)
	@echo " "
	@echo $(INCLS)

#build testprogram for terminal key map
test : ./obj/curses_test.o ./obj/test_drawer.o
	$(CC) ./obj/curses_test.o ./obj/test_drawer.o -o test $(NCURS)

reset : ./obj/curses_reset.o ./obj/test_drawer.o
	$(CC) ./obj/curses_reset.o ./obj/test_drawer.o -o reset_curses $(NCURS)

#Create the object directory before it is used, no error if not exists (order-only prereqs, will not rebuild objects if directory timestamp changes)
$(OBJS): | $(OBJDIR)
$(OBJDIR):
	@mkdir -p $(OBJDIR)

#Dependencies

obj/block.o : ./src/block.cpp $(BLOCKINCLS) $(SUPPORTINCLS) $(PAGEINCLS)
obj/block_page.o : ./src/block_page.cpp  $(SUPPORTINCLS) $(PAGEINCLS) $(BLOCKINCLS)
obj/block_set.o : ./src/block_set.cpp $(BLOCKINCLS) $(SUPPORTINCLS) $(PAGEINCLS)
obj/clipboard.o:./src/clipboard.cpp ./includes/clipboard.h
obj/config.o : ./src/config.cpp ./includes/drawer.h ./includes/page_items.h $(SUPPORTINCLS)
obj/deck_file.o : ./src/deck_file.cpp ./includes/deck_file.h $(SUPPORTINCLS)
obj/deck_support.o : ./src/deck_support.cpp ./includes/deck_support.h ./includes/config.h
obj/dividablepage.o : ./src/dividablepage.cpp $(PAGEINCLS) $(SUPPORTINCLS)
obj/drawer.o : ./src/drawer.cpp ./includes/drawer.h ./includes/page_items.h ./includes/deck_support.h ./includes/cmd_aliases.h
obj/io.o : ./src/io.cpp $(BLOCKINCLS) $(SUPPORTINCLS) $(PAGEINCLS)
obj/page.o : ./src/page.cpp ./includes/deck_support.h ./includes/page_items.h ./includes/drawer.h ./includes/page.h
obj/page_items.o : ./src/page_items.cpp $(PAGEINCLS) $(SUPPORTINCLS)
obj/page_overlays.o : ./src/page_overlays.cpp $(PAGEINCLS) $(SUPPORTINCLS)
obj/regionalisation_consts.o : ./src/regionalisation_consts.cpp ./includes/regionalisation_consts.h
obj/sub_page.o : ./src/sub_page.cpp $(PAGEINCLS) $(SUPPORTINCLS) $(BLOCKINCLS)

obj/test_drawer.o:./src/test_drawer.cpp ./includes/test_drawer.h
	$(CC) $(CFLAGS)  $< -o $@
obj/curses_test.o:./src/curses_test.cpp ./includes/curses_test.h ./includes/test_drawer.h ./src/test_drawer.cpp
obj/curses_reset.o:./src/curses_reset.cpp ./includes/curses_reset.h ./includes/test_drawer.h ./src/test_drawer.cpp

obj/%.o:./src/%.cpp
	$(CC) $(CFLAGS)  $< -o $@

.PHONY : tar tartest clean veryclean

tar:
	tar -cvzf Source.tgz $(SOURCE) $(INCLS) ./files/* Makefile input.deck

tartest:
	tar -cvzf SourceTest.tgz ./src/curses_test.cpp ./src/test_drawer.cpp ./includes/curses_test.h ./includes/test_drawer.h ./includes/deck_support.h Makefile Keypresses.txt

clean:
	@rm decker $(OBJS)

veryclean:
	@rm decker
	@rm -r $(OBJDIR)


