
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
2011
   
///////////////////////
   
11 rows and 13 columns of servos, with lacking servos at each corner, since the section of the tunnel isn't rectangular.
In the code the number of servos is taken equal to 13*11=143 in order to make the calculations easier.
The values corresponding to ghost servos are just useless.
Values are referenced in arrays 12*14
Row 0 and column 0 correspond to no servo, so as e.g. (1,1), etc.
   
I modified this library, so that it only contains functions that take a set of angles or speeds in input (given by algo) and return order to the servos (via SD84 and serial).

I tried to make the structure and the code clearer as they were before.
   
Florent Lachaussée
February 2013

*/

#include <activegrid.h>
#include <curses.h>
#include <time.h>


//set angles of the servos:
int  activegrid::setanglesII(double newangle[14][12]){
    int count = 0;
    //update angle array
    for(int row=0;row<12;row++){
        for(int col=0;col<14;col++){
            angle[col][row]=newangle[col][row];
            if (fabs(angle[col][row]) > 45) count++; // count how many angles exceed 45 degrees
        }
    }
    
    // safety check: no more than 80 paddles should have angles in excess of 45 degrees
    if (count > 80) {
        error(); // throw error method (below)
        return -23; // flee in terror rather than set potentially dangerous angles
    }

    for(int row=11;row>=1;row--){
        for(int col=1;col<=13;col++){
            if(servo[col][row]!=0){
                if(servo[col][row]>=10000 && servo[col][row]<=20000){
                    leftboard_orderofangle[(servo[col][row]-10001)]=angle[col][row];

                }
                else if(servo[col][row]>=20000 && servo[col][row]<=20200){
                    rightboard_orderofangle[(servo[col][row]-20001)]=angle[col][row];
                }
            }
        }
    }
    leftboard->setangles(1,67,leftboard_orderofangle,high_duty);
    rightboard->setangles(1,62,rightboard_orderofangle,high_duty);

 
  return 1;
}


// set the speeds of all of the servos:
int  activegrid::setspeeds(double newangleperstep[14][12]){
    
    for(int row=1;row<=11;row++){
        for(int col=1;col<=13;col++){
            if(servo[col][row]!=0){
                if(servo[col][row]>=10000 && servo[col][row]<=20000){
                    leftboard_anglestep_orderofangle[(servo[col][row]-10001)]=abs(newangleperstep[col][row]);
                                   }
                else if(servo[col][row]>=20000 && servo[col][row]<=20200){
                    rightboard_anglestep_orderofangle[(servo[col][row]-20001)]=abs(newangleperstep[col][row]);
                }
            }
        }
    }
    
    
    leftboard->setspeeds(1,67,leftboard_anglestep_orderofangle);
    rightboard->setspeeds(1,62,rightboard_anglestep_orderofangle);

        return 1;}



// getangle of one servo in memory
int  activegrid::getangle(int col, int row){
  return angle[col][row];
} 

double activegrid::getangleFromBoard(int col, int row){
  if(servo[col][row]!=0){
    if(servo[col][row]>=10000 && servo[col][row]<=20000){
      return leftboard->getangle((servo[col][row]-10000));
    }
    if(servo[col][row]>=20000 && servo[col][row]<=20200){
      return rightboard->getangle((servo[col][row]-20000));
    }
  }
  else
    cout << "In (" << col << "," << row << ") no servo exist" << endl;
  return -999;
}


// set all servos angles to zero degrees: 
int activegrid::opengrid(){
double newangle[14][12];
for(int row=0;row<12;row++){
  for(int col=0;col<14;col++){
    if(servo[col][row]!=0) newangle[col][row]=0;
    else newangle[col][row]=Fictive;
  }
 }
 setanglesII(newangle);
 return 1;
}


//show a map of the current angles
void activegrid::showallangles(){
  cout << "activegrid:showallangles(): Current angle configuration is shown: " << endl;
  cout << "row: " << endl;
  for(int row=11;row>=1;row--){
    cout << row << "\t||\t";
    for(int col=1;col<=13;col++){
      cout << angle[col][row] << "\t";
    }
    cout << endl;
  }
  for(int i=1;i<=120;i++) cout << "-";
  cout << endl;
  cout << "col: \t||\t";
  for(int col=1;col<=13;col++){
    cout << col << "  \t";
  }
  cout << endl;
}

//show a map of the servos
void activegrid::showservomap(){
  string leftname, rightname;
  leftname=leftboard->getpath();
  rightname=rightboard->getpath();
  cout << "activegrid:showservomap(): Map of the servos: \ntwo boards: \n 10001-10067: \t" << leftname << "\n 20001-20062:\t" << rightname << endl;
  for(int row=11;row>=1;row--){
    for(int col=1;col<=13;col++){
      cout << servo[col][row] << "\t";
    }
    cout << endl;
  }
}

//error message
int activegrid::error(){
    const char * header = "activegrid::error: ";
    cout << header << "Cannot set the required angles for safety reasons." << endl;
    //opengrid();
    return -1;
}
