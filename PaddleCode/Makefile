CPP = g++
FLAGS = -g -Wall -I./lib/ -D_LARGEFILE64_SOURCE -D_FILE_OFFSET_BITS=64 -O3
LIBDIR = ./lib
#LNK_OPTIONS = \
		-framework CoreFoundation\
		-framework IOKit\
		-L../../../../System/Library/Frameworks

all: \
	menuII \
	temp \
    Printing \

menuII: menuII.cpp
	$(CPP) $(FLAGS) $(LNK_OPTIONS) $(LIBDIR)/gbserialFlo.o $(LIBDIR)/SD84.o $(LIBDIR)/algo.o $(LIBDIR)/activegrid.o $(LIBDIR)/loaf.o -o menuII menuII.cpp

temp: temp.cpp
	$(CPP) $(FLAGS) $(LNK_OPTIONS) $(LIBDIR)/gbserialFlo.o $(LIBDIR)/SD84.o $(LIBDIR)/algo.o $(LIBDIR)/activegrid.o $(LIBDIR)/loaf.o -o temp temp.cpp

Printing: Printing.cpp
	$(CPP) $(FLAGS) $(LNK_OPTIONS) $(LIBDIR)/gbserialFlo.o $(LIBDIR)/SD84.o $(LIBDIR)/algo.o $(LIBDIR)/activegrid.o $(LIBDIR)/loaf.o -o Printing Printing.cpp