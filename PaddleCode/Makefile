CPP = g++
FLAGS = -g -Wall -I./lib/ -D_LARGEFILE64_SOURCE -D_FILE_OFFSET_BITS=64 -O3
LIBDIR = ./lib
#LNK_OPTIONS = \
		-framework CoreFoundation\
		-framework IOKit\
		-L../../../../System/Library/Frameworks

all: \
	menuII \
	ltfast\
	Printing

menuII: menuII.cpp
	@echo "==============================================================================================="
	@echo "Compiling menuII.cpp and linking with .o files to produce an executable file called menuII"
	@echo "==============================================================================================="
	$(CPP) $(FLAGS) $(LNK_OPTIONS) $(LIBDIR)/gbserialFlo.o $(LIBDIR)/SD84.o $(LIBDIR)/algo.o $(LIBDIR)/activegrid.o $(LIBDIR)/loaf.o $(LIBDIR)/testloaf.o $(LIBDIR)/pickCorrelations.o -o menuII menuII.cpp

ltfast: menuII.cpp
	@echo "==============================================================================================="
	@echo "Compiling menuII.cpp and linking with .o files to produce an executable file called ltfast"
	@echo "==============================================================================================="
	$(CPP) $(FLAGS) $(LNK_OPTIONS) $(LIBDIR)/gbserialFlo.o $(LIBDIR)/SD84.o $(LIBDIR)/algoltfast.o $(LIBDIR)/activegrid.o $(LIBDIR)/loaf.o $(LIBDIR)/testloaf.o $(LIBDIR)/pickCorrelations.o -o ltfast menuII.cpp

##temp: temp.cpp
##$(CPP) $(FLAGS) $(LNK_OPTIONS) $(LIBDIR)/gbserialFlo.o $(LIBDIR)/SD84.o $(LIBDIR)/algo.o $(LIBDIR)/activegrid.o $(LIBDIR)/loaf.o $(LIBDIR)/pickCorrelations.o -o temp temp.cpp

Printing: Printing.cpp
	@echo "==============================================================================================="
	@echo "Compiling Printing.cpp and linking with .o files to produce an executable file called Printing"
	@echo "==============================================================================================="
	$(CPP) $(FLAGS) $(LNK_OPTIONS) $(LIBDIR)/gbserialFlo.o $(LIBDIR)/SD84.o $(LIBDIR)/algo.o $(LIBDIR)/activegrid.o $(LIBDIR)/loaf.o $(LIBDIR)/pickCorrelations.o Printing.cpp -o Printing

## type make cleaner to delete all temp files, .o files, and executable binary files in your working directory and in lib
cleaner: clean
	rm -f *~ 
	rm -f \#*\# 
	rm -rf *.dSYM
	rm -f $(LIBDIR)/*~
	rm -f $(LIBDIR)/\#*\#
	rm -rf $(LIBDIR)/*.dSYM
	rm -f ../MatlabCode/*~

## type make clean to delete all .o files and executible files in your working directory and in lib
clean:
	rm -f *.o
	rm -f $(LIBDIR)/*.o
	rm -f menuII temp Printing
	rm -f $(LIBDIR)/testloaf $(LIBDIR)/testloaf2 $(LIBDIR)/testPickCorrelations $(LIBDIR)/manyfunction
