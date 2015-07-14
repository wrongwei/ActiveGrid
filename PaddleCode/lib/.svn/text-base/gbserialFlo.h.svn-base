
/* 
DON'T FORGET TO INSTALL THE DRIVER FIRST!!!

interface for serial communication

This class provides simple tools to communicate through 
a serial port.  
The tools initialize the port, send data, and receive 
data through the port.  

By using this class, whatever changes need to be made 
to the serial port communication code can be made once, 
by replacing the .cpp file corresponding to the header 
file, rather than throughout whatever other code needs 
to access the serial port.  

The code in the .cpp file is for a Mac.  
Much of the meat of the code comes from Apple's serial 
port example C program for communicating with a modem, 
which is on the web: 
http://developer.apple.com/samplecode/SerialPortSample

115,200 baud, 8 bit words, 2 stop bits, and no parity 
are all hard coded in the function OpenSerialPort, 
since this is what the SD84 board requires.  

An instance of this object makes a connection with a 
serial port when the user calls attach().  Once a 
successful attachment with a serial port is made, the 
instance is permanently connected to that serial port.  
Use the "portnametag" to select which serial port you 
connect to - the object only attaches to serial ports 
whose path name contains this string.  

Greg Bewley
September 2009
 */

#include <stdio.h>
#include <stdlib.h>
#include <memory.h>
#include <time.h>
#include <fcntl.h>
#include <sys/ioctl.h>
#include <errno.h>
#include <paths.h>
#include <termios.h>
#include <sysexits.h>
#include <sys/param.h>
#include <sys/select.h>
#include <sys/time.h>
#include <unistd.h>
#include <math.h>

#include <iostream>
using namespace std; 

#ifndef GBSERIAL
#define GBSERIAL

// by default, serial port path name must contain this: 
#define PORTNAMETAG "serial"

typedef unsigned char uchar; 

class GBSerialPort
{

public: 
	  // creator: 
	GBSerialPort(char * newportnametag = PORTNAMETAG); 
	  // destructor: 
	~GBSerialPort() {}; 
	
	// initialize the port; a single instance of a GBSerialPort can 
	// be succesfully attached to a serial port only once: 
	int  attach(char * newportnametag = PORTNAMETAG); 
	
	  // read something from the port: 
	int  receive(int numofbytes, uchar * inputbuffer); 
	
	  // write something to the port: 
	int  send(uchar * message, int numofbytes); 
	
	// return the pathname of the port: 
	string getnameofpath(); 

 private: 
	const char * devicePath;	
		
	int   fileDescriptor; 
	char  bsdPath[MAXPATHLEN]; 
	
	// the serial port we want has to match this tag: 
	char  portnametag[MAXPATHLEN]; 
	
}; 

inline GBSerialPort::GBSerialPort(char * newportnametag)
{
  
  // this means that you always look for something, 
  // never pick a random port: 
	if (strcmp(newportnametag, "") == 0)
	  strcpy(portnametag, PORTNAMETAG); 
	else
	  strcpy(portnametag, newportnametag); 
	
	fileDescriptor = -1; 
    bsdPath[0] = '\0'; 
}

#endif
