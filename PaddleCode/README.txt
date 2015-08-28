Short Explanation of the Active Grid Code

The scheme of the classes in /lib is the following:
	gbserialFlo.h 	(most inner software)
	SD84.h     	 	(software to communicate with an SD84 board (using gbserialFlo.h))
	activegrid.h	 	(software of the activegrid (using SD84.h))
	algo.h 			(algorithms for control of grid motions in spatial dimension)
	algo3d.h		(algorithms for control of grid motions in both spatial and temporal dimensions)
	pickCorrelations.h	(function pointer program; this is where correlation functions are coded and accessed)
	loaf.h			(custom-built and optimized 3D random number generating and storage structure)
	

The main idea is the following:
The real hardware-software communication is done in gbserialFlo.h and SD84.h. Communication between these low-level programs and high-level algorithms and user interface occurs through activegrid.h. Algorithms that generate servo angles and speeds are contained in algo.h and algo3d.h. Finally, the user interfaces (UI) for the basic/spatially correlated and temporally correlated program sets are run by menuII.cpp and menu3d.cpp, respectively.

To use the active grid in 2D (basic and spatially correlated movements) you work with menuII.cpp. Here you can choose different movements of the active grid. These currently run on an older implementation involving vectors to store angles and other intricate devices which were not possible to recycle for the 3D case.

Methods for these 2D movement algorithms are programmed in algo.cpp. To run these functions on the grid, first compile (using source source.txt - see below), then run the executable menuII.

The executable menu3d runs the user interface program menu3d.cpp. Selecting choice 10 from this interface grants the user access to the spatial and temporal correlations coded in algo3d.cpp.

pickCorrelations.cpp contains the functions accessed by function pointers in algo.cpp and algo3d.cpp (see comments in pickCorrelations.cpp).

loaf.cpp is used by algo3d.cpp to store random numbers for correlation computations in a more efficient way (see comments in loaf.cpp).

Also included are a few test clients that allow more straightforward debugging of the loaf and pickCorrelations implementations.


Other remarks:

-safety: For safety reasons (so that the active grid is not too closed so that it doesn't b break) the method "checkangles()" exists in activegrid.h which should be called before every movement of the servos. At the moment this method works in this way:
	 1. it counts the number of servos with an angle > 80 degrees
	 2. if this number is > 80 then the program will refuse to set the given angles, print a message, and wait for another set to come in

-how to start the program:
     source source.txt 
     ./menuII
     (Attention: There are two Makefiles at the moment: one for all classes in /lib and one for ./menuII. If you use source.txt both Makefiles are called and everything is fine. If you have changed something in the classes in /lib you have to run both Makefiles. Also, both executables - menuII and menu3d - are compiled in the makefiles, so there is no need to compile them separately.)

-------------------

First Steps to install the temperature sensors of the active grid:
I bought 10-30 temperature sensors which can be read out via SD84 board. In the file /TemperatureSensorForTheActiveGrid you can see the first calibration which has been taken with the program temp.cpp
This program just returns the voltage of a chosen analogue channel of SB84 board.

Aim: You should install several such temperature sensors on the active grid and use several analogue channels of SD84 (attention: SD84 is programmed in a way that you have to use the first analogue channel for the first temperature sensor and so on -> so you have to change the channels of the servos which are plugged in a analogue channel at the moment, so you have to change the mapping in the code.)

-------------------

If you have any questions to these programs/classes write me an email:
florian.koehler.email@googlemail.com

Florian Koehler,
September 2011


Questions regarding function pointers, loaf, and temporal correlations may be directed to Kevin Griffin and Nathan Wei:
kevinpg@princeton.edu
nwei@princeton.edu
README updated by Nathan Wei, August 2015


technical info: 

digital servos with brushless motor:
BLS152RB by Futaba
in total: 129 servos
67 on board a
62 on board b

board:
SD84 - 84 Channel Servo Driver Module
power supply: 5V

