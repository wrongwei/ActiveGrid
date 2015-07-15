
/* 
library to talk to the SD84 servo control board

This class provides easy access to the basic 
functions of the SD84 servo control board.  

The SD84 sends out pulse width modulated signals.  
The servos take the width of the pulse as a 
command to assume a certain angular position.  
There is no way to get feedback from the servo 
on its actual position.  
All the functions provided here interact only 
with the SD84 to set and retrieve the commanded 
servo angles.  

The SD84 is a USB device with an FTDI USB to 
serial converter chip on board: 
http://www.robot-electronics.co.uk/htm/sd84tech.htm

FTDI provides a driver that makes the SD84 look 
like a serial device: 
http://www.ftdichip.com/FTDrivers.htm

The class defined in gbserial.h provides a 
platform independent way to talk through the 
serial port.  

servo numbers start with 1 (not zero).  

Greg Bewley
September 2009
 */

#include <iostream>
#include <string>
#include <cmath>

using namespace std; 

#include "gbserialFlo.h"
#include "SD84constants.h"

//#ifndef SD84
//#define SD84

  // servo angle conversion constants: 
  // angle  90 corresponds to counterclockwise extreme.  
  // angle -90 corresponds to clockwise extreme.  
#define SERVOCENTER  1500  // [microseconds] pulse width to give 0 degrees
                           // 1500 good for Futaba BLS 152
#define MSP180DEG    1200  // [microseconds per 180 degrees] pulse width
                           // 1200 good for Futaba BLS 152

#define DEFAULTGAIN  0.01  // default gain for "smooth moves"

typedef unsigned char uchar; 

class SD84
{
public: 
	  // creator: 
  SD84(char * portnametag = (char *) "", 
	     int newservocenter = SERVOCENTER, 
	     int newmsp180deg   = MSP180DEG); 
	  // destructor: do nothing (the serial port shuts itself down...  ).  
	~SD84() {}; 
	
	  //---- functions to set the angular positions of servos: 
	  // set a single servo to an angle in degrees: 
	int setangle(int servo_number, double angle); 
	// set all servos to the same position: 
	int setangles(double angle); 
	  // set several servos to different positions: 
	int setangles(int first_servo, int num_of_servos, double angles[], bool high_duty);
	  // get a single servo's current commanded angle in degrees from the board: 
	double getangle(int servo_number); 
	// get Temperature (does not work)
	uchar * getTemp(int servo_number); 
	// get Analog Value
	float getAnalog(int numberOfServos); 
	
	  //---- functions to set the rotation rates of servos: 
	  // set the speed of a single servo: 
	int setspeed(int servo_number, double rotationrate); 
	  // set the speed of *all* servos in rotations per second: 
	int setspeeds(double rotationrate); 
	  // set several servos to different speeds: 
	int setspeeds(int first_servo, int num_of_servos, double rotationrates[]); 
	  // The board provides no means to get speeds.  
	  // These functions do not communicate with the board.  
	  // The first function returns the last speed set *globally* using setspeeds (double), 
	  // and will not return the correct speed of an individual servo whose speed was last 
	  // set with setspeeds(int, int, double).  
	  // The second function returns the speed of the specified servo, 
	  // which will always correspond to the last speed the servo was set to.  
	double getspeed(); 
	double getspeed(int servo_number); 
	  // maximum controllable speed is limited by the 1-byte size of the instruction: 
	double getmaxspeed(); 
	
	  //---- miscellaneous other functions: 
	  // set the mode of a channel: 
	int  setmode(int channel, ChannelMode mode); 
	  // set the mode of *all* channels: 
	int  setmodes(ChannelMode mode); 
	
	  // get the value of a digital input channel: 
	int  getinput(int channel); 
	
	  // test whether the board responds correctly using this function, 
	  // processor 1, 2, 3, or 4 should return 84: 
	int  getversion(int processornumber); 
	
	  // send nonsense to the board: 
	void sendnonsense(); 
	
	  //---- functions derived from the above ones, or that operate only locally: 
	  // check whether the port initialized ok: 
	int  checkport(); 
	
	  // get the path to the serial port: 
	string getpath(); 
	
	  // set a new position, and move to this new position with a 
	  // velocity proportional to the difference between positions: 
	void setgain(double newgain); 
	int  smoothmove(int servo_number, double angle); 

    
private: 
	
	  // we communicate with the board through this: 
	GBSerialPort thePort; 
	  // whatever class you use for thePort, it should include the functions: 
	  // int initialize(), int send(char *, int), int receive(int, char *), 
	  // and void getpath(char *).  
	
	  // true if there is something wrong with communications: 
	bool portbroken; 
	
	// path name of the port: 
	string portpath; //[MAXPATHLEN]; 
	
	uchar outputbuffer[256]; 
	uchar inputbuffer[100]; 
	
	  // the last speed set globally, and for each servo: 
	double speed, speeds[MAXSERVOS]; 
	
	  // parameter for smoothmove: 
	double gain; 
	
	  // the servo motion parameters: 
	double servocenter, msp180deg; 
	
	  // for conversion between angles, in degrees, 
	  // and pulse widths, in microseconds.  
	int angle2width(double angle); 
	double width2angle(int width); 
	
	  // convert between changes in pulse widths, in microseconds per cycle, 
	  // to rotation rates, in rotations per second.  
	int rr2dwpc(double rr); 
	double dwpc2rr(int dw); 

}; 

//#endif

inline SD84::SD84(char * portnametag, int newservocenter, int newmsp180deg)
{
	const char * header = "SD84 constructor: "; 
	
	// set the servo parameters, if the user gives no inputs, to the default values: 
	servocenter = newservocenter;  // default SERVOCENTER
	msp180deg   = newmsp180deg;    // default MSP180DEG
	
	cout << header << "calling GBSerialPort.attach...  " << endl; 
	
	if (thePort.attach(portnametag) != 1){
	  portbroken = true; 
	  cout << header << "serial port initialization failed!  " << endl; 
	}
	else
	  {
	    portbroken = false; 
	    
	    // the first three bytes sent to the card are always the same, 
	    // to sync up the communication.  DON'T CHANGE THEM!  
	    outputbuffer[0] = 0xAA; 
	    outputbuffer[1] = 0xA0; 
	    outputbuffer[2] = 0x55; 
	    
	    cout << header << "contacting the SD84 board...  " << endl; 
	    
	    // check that we can talk to the board: 
	    if (getversion(1) == -1){
		cout << header << "problem communicating with the board.  " << endl; 
		portbroken = true; 
	      }
	
	    portpath=thePort.getnameofpath();
    
	    cout << endl << header << endl; 
	    cout << "  Welcome to the SD84 84-channel servo-controller-board controller.  " << endl; 
	    cout << "  Serial port at " << portpath << endl; 
	    cout << "  Using a " << servocenter << " microsecond pulse for servo center, " << endl; 
	    cout << "  and +/- " << msp180deg/2 << " microseconds for +/- 90 degrees.  " << endl; 
	    cout << "  Maximum speed is " << getmaxspeed() << " rotations per second.  " << endl; 
	    cout << "  Minimum speed is " << dwpc2rr(1) << " rotations per second...  " << endl; 
	    
	    // make sure all channels are set to "servo mode": 
	    setmodes(SERVO_MODE); 
	    cout << "Mode is set" << endl; 
	    
	    // the board starts up with speed set to 0 (full speed) on all channels, 
	    // but let's make sure it is so: 
	    cout << endl << header; 
	    setspeeds(2*getmaxspeed()); 
	    cout << "Speed is set" << endl; 
	    
	    // set the default gain for smoothmove: 
	    gain = DEFAULTGAIN; 
	    
	    cout << endl; 
	  }
}

// convert angle in degrees to pulse width in microseconds: 
inline int SD84::angle2width(double angle)
{
  return (int) round((angle * msp180deg / 180.0) + servocenter); 
}

  // convert pulse width in microseconds to angle in degrees: 
inline double SD84::width2angle(int width)
{
	return (double) (width - servocenter) * (180.0 / msp180deg); 
}

  // convert rotation speed to change in pulse width per cycle in microsecs: 
inline int SD84::rr2dwpc(double rr)
{
  return DWMAX;
  //return (int) round(rr * 2.0 * msp180deg * CYCLETIME / 1000.0); 
  //return (int) round(rr * CYCLETIME / (msp180deg * 1000.0)); //flo
}

  // convert change in pulse width per cycle to rotation speed in rot/sec: 
inline double SD84::dwpc2rr(int dw)
{
	return (double) (dw * 1000.0) / (2.0 * msp180deg * CYCLETIME); 
}
