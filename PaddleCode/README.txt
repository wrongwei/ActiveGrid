Short Explanation of the Active Grid Code

The scheme of the classes in /lib is the following:
    	   gbserialFlo.h 	 (most inner software)
	   SD84.h     	 (software to communicate with an SD84 board (using gbserialFlo.h))
	   activegrid.h	 (software of the activegrid (using SD84.h))

The main idea is the following:
If you want to write software for the active grid you only use the methods from activegrid.h. So you don't need to worry about the real hardware-software communication done in gbserialFlo.h and SD84.h.

To use the activegrid at the moment you work with menuII.cpp. Here you can choose different movements of the active grid. (option "12" is the most complex movement at the movement).

Methods for movement algorithms are programmed in menu.h.

That's it.



Other remarks:

-safety: For safety reasons (so that the active grid is not too closed so that it doesn't b break) the method "checkangles()" exists in activegrid.h which should be called before every movement of the servos. At the moment this method works in this way:
	 1. it counts the number of servos with an angle > 45 degrees
	 2. if this number is > 80 you should look for another servo movement

-how to start the program:
     source source.txt 
     ./menuII
     (!!!attention: There are two Makefiles at the moment: one for all classes in /lib and one for ./menuII. If you use source.txt both Makefiles are called and everything is fine If you have changed something in the classes in /lib you have to run both Makefiles!!!)

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


technical info: 

digital servos with brushless motor:
BLS152RB by Futaba
in total: 129 servos
67 on board a
62 on board b

board:
SD84 - 84 Channel Servo Driver Module
power supply: 5V

