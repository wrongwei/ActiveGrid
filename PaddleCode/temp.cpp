#include "SD84.h"
#include <iostream> 
         
using namespace std;  
    
int main(){  
  SD84 * board;  
  board = new SD84((char *)"/dev/cu.usbserial-A2001mHT",1500,1200);   //initialize board 2
  cout << " --- START --- " << endl;
  float voltageSignal;  
  int analogChannel=21;     
  float meanVoltSignal=0;
  float meanVoltSignalSquare=0;
  for(int i=0;i<40;i++){   
    sleep(1);
    voltageSignal = board->getAnalog(analogChannel);   
    cout << "measurement: " << i << " (40)   Channel: " << analogChannel << "  voltage: " << (voltageSignal*10)/2 << " guessed temperature: [C] " << (voltageSignal)/2 << endl;
    meanVoltSignal+=voltageSignal*10/2;
    meanVoltSignalSquare+=(voltageSignal*10*voltageSignal*10/4);
  }
  cout << "mean: " << (float)meanVoltSignal/40 << "  var: " << sqrt((meanVoltSignalSquare/40)-((meanVoltSignal*meanVoltSignal)/(40*40))) << endl;
  return 0;   
} 
 


// 4 mA is 0 mbar -> 0.1 Volt -> 25 value
// 20 mA is 15 mbar -> 0.5 Volt -> 120 value (80 mA would be great for 20 mA)
// resistance: 25 ohm

/*
1500
low high
DC 05

2^7=128

2^8=256
2^9=512
2^10=1024
2^11=2048

1024   1280  1408  1472  1488  1496  1500
2^10 + 2^8 + 2^7 + 2^6 + 2^4 + 2^3 + 2^2
 
101 11011100
5   4+8+16+64+128=220
5   DC
high low
  
low first and than high
*/
