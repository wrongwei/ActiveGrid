/*
library to control the active grid.  

This class provides the tools necessary to make 
the active grid do what you want it to do.  

The active grid is composed of N=129 servos, 
which can be set to angles between -90 and 90 
degrees.  
angle = 0 corresponds to grid flap parallel to 
the flow.  

The servos are arranged on a grid.  

if any error is encountered, all paddles are opened.  


Florian Koehler
florian.koehler.email@googlemail.com
April 2011
 
////////////////////////
 11 rows and 13 columns of servos, with lacking servos at each corner, since the section of the tunnel isn't rectangular.
 In the code the number of servos is taken equal to 13*11=143 in order to make the calculations easier.
 The values corresponding to ghost servos are just useless.
 Values are referenced in arrays 12*14
 Row 0 and column 0 correspond to no servo, so as e.g. (1,1), etc.

 
 I modified this library, so that it only contains functions that take a set of angles in input (given by algo) and return order to the servos (via SD84 and serial).

 I tried to make the structure and the code clearer as they were before.
 
 Florent Lachaussée
 florent.lachaussee@ens.fr
 February 2013

 */


#include <iostream>
#include <string>
#include <cmath>

using namespace std; 

#include "SD84.h"

#define numberOfServos 143 // nb of servos on the rectangular grid
#define Fictive -100 //value of angle for non-existent paddle, e.g. (1,1), (1,13)

//#ifndef ACTIVEGRID
//#define ACTIVEGRID

class activegrid
{
 public: 
  // constructor: 
  activegrid(); 
  // destructor: 
  ~activegrid() {}; 
   

 // set the positions of the servos:
    
  int  setanglesII(double newangle[14][12]); //here the angles are actually set

    
// set the speeds of all of the servos:
    int  setspeeds(double newangleperstep[14][12]);

    
  // set all servo angles to zero degrees: 
  int  opengrid();
    
 // get the position of one servo (the one in memory, not the effective one, to which we can't get through):
   int  getangle(int col, int row);
   double getangleFromBoard(int col, int row);

//show a map of the current angles
    void  showallangles();
    
//show a map of the servos
    void  showservomap();

//map of servo-board-numbers and col/row:
    int servo[14][12];
    
//store our own local copy of the last-set servo angles:
    double angle[14][12];

//to know what the amplitude is, when varying it
    bool high_duty;
    
 private: 
  

  double leftboard_orderofangle[67];//a
  double rightboard_orderofangle[62];//b
  double leftboard_anglestep_orderofangle[67];//a
  double rightboard_anglestep_orderofangle[62];//b


  // two boards control the grid: 
  SD84 * leftboard;// controls the left side of the grid, viewed upstream.  
  SD84 * rightboard;// controls the right side of the grid, viewed upstream.  
  
  double movetime;  // set time it takes for servos to get to a new set angle.
	
  int error();      // opens all paddles and returns -1.  
}; 

//#endif

inline activegrid::activegrid(){

  leftboard = new SD84("/dev/cu.usbserial-A2001mHD",1500,1200);   //initialize board 1
  rightboard = new SD84("/dev/cu.usbserial-A2001mHT",1500,1200);   //initialize board 2
    // The digital output is on channel 84 of rightboard
    // names for the boards in the hall : LB : "/dev/cu.usbserial-A2001mHD"
    // RB : "/dev/cu.usbserial-A2001mHT"
    // in the office : LB: "/dev/cu.usbserial-A2001mH5"
    // RB : "/dev/cu.usbserial-A2001mHa"
  //!!!!!!!!!!!!!!!! BEWARE THE NAMES OF THE DEVICES !!!!!!!!!!!!!!!!!!!!
    
    
  //fill angle and servo 
  for(int row=0;row<12;row++){
    for(int col=0;col<14;col++){
      angle[col][row]=0;
    servo[col][row]=0;
    }
  }

    //mapping of servos : non-zero value of servo for existent paddles.
    servo[2][1]=10038;servo[3][1]=20062;servo[4][1]=10047;servo[5][1]=20059;servo[6][1]=10055;servo[7][1]=20054;servo[8][1]=10061;servo[9][1]=20047;servo[10][1]=10065;servo[11][1]=20039;servo[12][1]=10067;
    servo[1][2]=10028;servo[2][2]=20061;servo[3][2]=10037;servo[4][2]=20058;servo[5][2]=10046;servo[6][2]=20053;servo[7][2]=10054;servo[8][2]=20046;servo[9][2]=10060;servo[10][2]=20038;servo[11][2]=10064;servo[12][2]=20030;servo[13][2]=10066;
    servo[1][3]=20060;servo[2][3]=10027;servo[3][3]=20057;servo[4][3]=10036;servo[5][3]=20052;servo[6][3]=10045;servo[7][3]=20045;servo[8][3]=10053;servo[9][3]=20037;servo[10][3]=10059;servo[11][3]=20029;servo[12][3]=10063;servo[13][3]=20021;
    servo[1][4]=10018;servo[2][4]=20056;servo[3][4]=10026;servo[4][4]=20051;servo[5][4]=10035;servo[6][4]=20044;servo[7][4]=10044;servo[8][4]=20036;servo[9][4]=10052;servo[10][4]=20028;servo[11][4]=10058;servo[12][4]=20020;servo[13][4]=10062;
    servo[1][5]=20055;servo[2][5]=10017;servo[3][5]=20050;servo[4][5]=10025;servo[5][5]=20043;servo[6][5]=10034;servo[7][5]=20035;servo[8][5]=10043;servo[9][5]=20027;servo[10][5]=10051;servo[11][5]=20019;servo[12][5]=10057;servo[13][5]=20012;
    servo[1][6]=10010;servo[2][6]=20049;servo[3][6]=10016;servo[4][6]=20042;servo[5][6]=10024;servo[6][6]=20034;servo[7][6]=10033;servo[8][6]=20026;servo[9][6]=10042;servo[10][6]=20018;servo[11][6]=10050;servo[12][6]=20011;servo[13][6]=10056;
    servo[1][7]=20048;servo[2][7]=10009;servo[3][7]=20041;servo[4][7]=10015;servo[5][7]=20033;servo[6][7]=10023;servo[7][7]=20025;servo[8][7]=10032;servo[9][7]=20017;servo[10][7]=10041;servo[11][7]=20010;servo[12][7]=10049;servo[13][7]=20005;
    servo[1][8]=10004;servo[2][8]=20040;servo[3][8]=10008;servo[4][8]=20032;servo[5][8]=10014;servo[6][8]=20024;servo[7][8]=10022;servo[8][8]=20016;servo[9][8]=10031;servo[10][8]=20009;servo[11][8]=10040;servo[12][8]=20004;servo[13][8]=10048;
    servo[2][9]=10003;servo[3][9]=20031;servo[4][9]=10007;servo[5][9]=20023;servo[6][9]=10013;servo[7][9]=20015;servo[8][9]=10021;servo[9][9]=20008;servo[10][9]=10030;servo[11][9]=20003;servo[12][9]=10039;
    servo[3][10]=10002;servo[4][10]=20022;servo[5][10]=10006;servo[6][10]=20014;servo[7][10]=10012;servo[8][10]=20007;servo[9][10]=10020;servo[10][10]=20002;servo[11][10]=10029;
    servo[4][11]=10001;servo[5][11]=20013;servo[6][11]=10005;servo[7][11]=20006;servo[8][11]=10011;servo[9][11]=20001;servo[10][11]=10019;

    
//    servo[1][1]=20063; // non-existent paddle, but used as output to know the amplitude of motion
    //(for correlatedMovement_steps, motion with amplitude varying with time)
    

  //starting position: set all angles to 0
  opengrid();  
};
