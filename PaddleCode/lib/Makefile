CPP = g++
FLAGS = -g -Wall -I./ -DBINARY -O3
PADDLEDIR = ..
#LNK_OPTIONS = \
		-framework CoreFoundation\
		-framework IOKit\
		-L../../../../System/Library/Frameworks

all:  \
	gbserialFlo \
	SD84 \
	activegrid \
	algo \
	loaf \
	testloaf \
	pickCorrelations

#	gbserialdummy \
#	PosixSignalDispatcher \
#	serialstream \
#	serialstreambuf \
#	serialport \

#PosixSignalDispatcher: PosixSignalDispatcher.cpp PosixSignalDispatcher.h
#	$(CPP) $(FLAGS) -c PosixSignalDispatcher.cpp

#serialstream: SerialStream.cc SerialStream.h
#	$(CPP) $(FLAGS) -c SerialStream.cc

#serialstreambuf: SerialStreambuf.cc SerialStreambuf.h
#	$(CPP) $(FLAGS) -c SerialStreamBuf.cc

#serialport: SerialPort.cpp SerialPort.h
#	$(CPP) $(FLAGS) -c SerialPort.cpp

#gbserialdummy: gbserialdummy.cpp gbserialdummy.h
#	$(CPP) $(FLAGS) -c gbserialdummy.cpp

#gbserial: gbserial.cpp gbserial.h
#	$(CPP) $(FLAGS) -c gbserial.cpp

gbserialFlo: gbserialFlo.cpp gbserialFlo.h
	$(CPP) $(FLAGS) -c gbserialFlo.cpp

SD84: SD84.cpp SD84.h
	$(CPP) $(FLAGS) -c SD84.cpp

activegrid: activegrid.cpp activegrid.h
	$(CPP) $(FLAGS) -c activegrid.cpp

algo: algo.cpp algo.h loaf.h 
	$(CPP) $(FLAGS) -c algo.cpp

loaf: loaf.cpp loaf.h
	$(CPP) $(FLAGS) -c loaf.cpp

testloaf: testloaf.cpp
	$(CPP) $(FLAGS) -c testloaf.cpp

pickCorrelations: pickCorrelations.cpp pickCorrelations.h
	$(CPP) $(FLAGS) -c pickCorrelations.cpp

## type make cleaner to delete all temp files, .o files, and executable binary files in your working directory and in PaddleCode
cleaner: clean
	rm -f *~ 
	rm -f \#*\# 
	rm -rf *.dSYM
	rm -f $(PADDLEDIR)/*~
	rm -f $(PADDLEDIR)/\#*\#
	rm -rf $(PADDLEDIR)/*.dSYM
	rm -f ../../MatlabCode/*~

## type make clean to delete all .o files and executible files in your working directory and in PaddleCode
clean:
	rm -f *.o
	rm -f $(PADDLEDIR)/*.o 
	rm -f $(PADDLEDIR)/menuII $(PADDLEDIR)/temp $(PADDLEDIR)/Printing
	rm -f testloaf testloaf2 manyfunction testPickCorrelation