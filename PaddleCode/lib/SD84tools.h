
/* SD84 library of tools for using an SD84 board with servos.  


Greg Bewley
October 2009
 */

#include "SD84.h"

void pollposition(SD84 *board, int servnum, double interval, int count)
{
	double angle; 
	for (int i = 0; i < count; i++)
	{
		angle = board->getangle(servnum); 
		cout << i*interval << " " << angle << endl; 
		usleep(1000000*interval); 
	}
}

int SD84::smoothmove(int servo_number, double angle)
{
	  // get the current commanded position: 
	double oldangle = getangle(servo_number); 
	
	double tempspeed = gain*abs(oldangle - angle); 
//	cout << "SD84::smoothmove: speed for this move: " << tempspeed << " [rot/sec].  " << endl; 
	if (setspeed(servo_number, tempspeed) != 1)
		return -1; 
	
	if (setangle(servo_number, angle) != 1)
		return -1; 
	
	  // return the servo speed to the globally set one: 
//	if (setspeed(servo_number, speed) != 1)
//		return -1; 
	
	return 1; 
}
