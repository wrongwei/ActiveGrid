
/*
 menu program to choose all kind of options via terminal
 just run it with:
 source source.txt
 ./menuII
 
 Florian Koehler
 florian.koehler.email@googlemail.com
 September 2011
 
 /////////
 
 Always the same goal, with a few changes in the names and new options.
 
 Florent Lachaussée
 florent.lachaussee@ens.fr
 March 2013
 
 Driver from FTDI for usb/serial port must be installed.
 check names of devices in /dev and in activegrid.h
 
 */
/*--------------------------------------------------------------------*/

#include <iostream>
#include <sstream>
#include <unistd.h>
#include <time.h>
#include "lib/algo.h"

/*--------------------------------------------------------------------*/

//Headers needed for signal handling
#define _GNU_SOURCE
#include <stdio.h>
#include <stdlib.h>
#include <signal.h>

/*--------------------------------------------------------------------*/

using namespace std;

//functions:
void wait(float seconds);
int main(int argc, char * const argv[]);
double angle_check(double angle);
double check_frequency(double frequency, double amplitude);


void wait(float seconds) {
    clock_t endwait;
    endwait = clock () + seconds * CLOCKS_PER_SEC ;
    while (clock() < endwait) {}
}

char * myLoaf = NULL;

/*--------------------------------------------------------------------*/
// Install a signal handler to make sure that when the user types
// CTRL-C that memory is freed before the program is terminated
// If you want to terminate the program immediately, type CTRL-backslash
static void mySignalHandler(int iSignal){
  // Free any used memory
  cout << "\nFreeing memory in use\n";
  if (myLoaf != NULL){
    free(myLoaf);
  }
  
  // Kill the program
  cout << "Killing program by sending the signal SIGKILL, which equals 9. ";
  fflush(NULL);
  raise(SIGKILL);
}

/*--------------------------------------------------------------------*/

int main (int argc , char * const argv[]) {
    // create an algo object
    algo alg;
    int choice=-1;
    int i=1;
    
    // Make sure CTRL-C signals are not blocked.
    void (*pfRet)(int);
    sigset_t sSet;
    int iRet;
    iRet = sigemptyset(&sSet);
    if (iRet == -1) {perror(argv[0]); exit(EXIT_FAILURE); }
    iRet = sigaddset(&sSet, SIGINT);
    if (iRet == -1) {perror(argv[0]); exit(EXIT_FAILURE); }
    iRet = sigprocmask(SIG_UNBLOCK, &sSet, NULL);
    if (iRet == -1) {perror(argv[0]); exit(EXIT_FAILURE); }

    // Install mySignalHandler as the handler for CTRL-C signals.
    pfRet = signal(SIGINT, mySignalHandler);
    if (pfRet == SIG_ERR) {perror(argv[0]); exit(EXIT_FAILURE); }
    
    // Welcome
    cout << endl << "Welcome to " << argv[0] << endl;
    cout << "This is the program to run the software of the activegrid." << endl;
    
    while(i!=0){
      cout << "Menu:\n"
	" 1 - see current angle configuration\n"
	" 2 - see current servo map\n"
	" 3 - open the grid\n\n"
	" 4 - set all servos to the same angle\n"
	" 5 - set one certain servo to one certain angle\n"
	" 6 - set one col to one certain angle\n"
	" 7 - set one row to one certain angle\n"
	" 8 - close row 1-5 and open 6-11\n"
	" 9 - get angle of one servo\n\n"
	" 10 - all periodic\n"
	" 11 - chaotic 1 (piecewise periodic)\n"
	" 12 - chaotic 2 (random)\n"
	" 13 - chaotic correlated\n"
	" 14 - chaotic correlated (2 alternating amplitudes)\n"
	" 15 - chaotic correlated (periodic pattern)\n\n"
	" 16 - chaotic correlated in space and correlated in time\n\n"
	" 17 - test the paddles by opening and closing each row and then each column\n\n"
	" 18 - set boundary paddles to constant angle\n\n"
	" 19 - test loaf object\n\n"
	" 20 - end program\n\n";
      cout << "!!! Numbers from 11 to 16 are the most interesting movements !!!" << endl;
      cin >> choice;
      cout << "your choice is  " << choice << "\n";
    
      if((int)choice==1){
	alg.grid.showallangles();
      }
      
      else if((int)choice==2){
	alg.grid.showservomap();
      }
    
      else if((int)choice==3){
	alg.grid.opengrid();
      }
    
      else if((int)choice==4){
	double angle=0;
	cout << "Angle? (-90 - 90)";
	cin >> angle;
	angle=angle_check(angle);
	cout << endl;
	alg.setsameangletoallservos(angle);
      }
        
      else if((int)choice==5){
	int row=0;
	int col=0;
	double angle=0;
	cout << "Col? (1-13)";
	cin >> col;
	cout << endl;
	cout << "Row? (1-11)";
	cin >> row;
	cout << endl;
	cout << "Angle? (-90 - 90)";
	cin >> angle;
	angle=angle_check(angle);
	cout << endl;
	alg.setangleofoneservo(col,row,angle);
      }
        
      else if((int)choice==6){
	int col=0;
	double angle=0;
	cout << "Col? (1-13)";
	cin >> col;
	cout << endl;
	cout << "Angle? (-90 - 90)";
	cin >> angle;
	angle=angle_check(angle);
	cout << endl;
	alg.setanglesofonecolumn(col,angle);
      }
        
        else if((int)choice==7){
            int row=0;
            double angle=0;
            cout << "Row? (1-11)";
            cin >> row;
            cout << endl;
            cout << "Angle? (-90 - 90)";
            cin >> angle;
            angle=angle_check(angle);
            cout << endl;
            alg.setanglesofonerow(row,angle);
        }
        
        else if((int)choice==8){
            alg.openhalfclosehalf();
        }
        
        //Get angle of one servo
        else if((int)choice==9){
            int row=0;
            int col=0;
            cout << "Col? (1-13)";
            cin >> col;
            cout << endl;
            cout << "Row? (1-11)";
            cin >> row;
            cout << endl;
            if(alg.grid.getangle(col,row)!= Fictive )
                cout << "Angle of (" << col << "," << row << "): " << alg.grid.getangle(col,row) << endl;
            else
                cout << "(" << col << "," << row << ") does not exist!" << endl;
        }
        
        else if((int)choice==10){
            double angle=0;
            double frequency=0;
            cout << "maximal angle? (-90 - 90) ";
            cin >> angle;
            angle = angle_check(angle);
            cout << endl;
            cout << "Frequency [Hz]? ";
            cin >> frequency;
            //frequency=check_frequency(frequency,angle);
            if (frequency!=-1) {
                cout << endl;
                wait(1.5);
                
                alg.allperiodic(angle,frequency);
            }
        }
        
        else if((int)choice==11){
            int combine;
            int constant;
            int option = 1;
            cout << "Combine? (1,2,3,4,5,6,7)";
            cin >> combine;
            cout << "Should the area be kept constant? (1,0) ";
            cin >> constant;
            //constant=0;
            
            alg.chaoticMovement(combine,constant,option);
        }
        
        else if((int)choice==12){
            int combine;
            int constant;
            int option = 2;
            float range;
            
            cout << "Range of oscillations? (0->90 degrees) ";
            cin >> range;
            while ( range < 0 || range > 90){
                cout << "Range should be between 0 and 90 degrees!! Try an acceptable value!" <<endl;
                cout << "\nRange of oscillations? (0->90 degrees) ";
                cin >> range;
            }
            
            alg.max_speed = 2 * range / 5;
            alg.max_angle = range;
            alg.min_angle = - range;
            
            cout << "Combine? (1,2,3,4,5,6,7)";
            cin >> combine;
            
            cout << "Should the area be kept constant? (1,0) ";
            cin >> constant;
            //constant=0;
            
            alg.chaoticMovement(combine,constant,option);
        }
        
        else if((int)choice==13){
            float sigma;
            int constant;
            int mode;
            float target_rms;
            int mrow;
            int mcol;
            float alpha;
            double height;
            
            cout << "Choose the correlation function: \n 1 - Gaussian \n 2 - 1/r^2 \n 3 - 1/|r|^3 \n 4 - 1/r^4 \n 5 - top hat \n 6 - top hat with one main paddle \n 7 - top hat with random main paddle \n 8 - top hat with long tail \n 9 - triangular function \n";
            cin >> mode;
            while (mode > 9 || mode < 1) {
                cout << "Invalid choice! Try again! \n";
                cout << "\n Choose the correlation function: \n 1 - Gaussian \n 2 - 1/r^2 \n 3 - 1/|r|^3 \n 4 - 1/r^4 \n 5 - top hat \n 6 - top hat with one main paddle \n 7 - top hat with random main paddle \n 8 - top hat with long tail \n 9 - triangular function \n";
                cin >> mode;
            }
            
            //NEED TO FIX THIS LATER TO MAKE RANGE OF CORRELATION MORE GENERAL FOR ALL FUNCTIONS
            
            if ((int)mode == 1){ // for gaussian
                cout << "Sigma? ";
                cin >> sigma;
            } else if ((int)mode == 6){ // for top hat with one main paddle
                cout << "Choose main paddle: Row? ";
                cin >> mrow;
                while (mrow < 0 || mrow > 12) {
                    cout << "choose an integer between 1 and 11 (included)";
                    cin >> mrow;
                };
                mrow = mrow - 1; //indexing
                cout << "Column? ";
                cin >> mcol;
                while (mcol < 0 || mcol > 14) {
                    cout << "choose an integer between 1 and 13 (included)";
                    cin >> mcol;
                };
                mcol = mcol - 1; //indexing
                
                cout << "Correlation distance? (paddles) ";
                cin >> sigma;
                
            } else if ((int)mode == 5 || (int)mode == 7 || (int)mode == 9) { // for top hat, top hat with random main paddle, and triangle
                cout << "Correlation distance? (paddles) ";
                cin >> sigma;
            } else if ((int)mode == 8){ // for top hat with long tail
                cout << "Perfect correlation distance? (paddles)";
                cin >> alpha;
                cout << "Total correlation distance? (paddles)";
                cin >> sigma;
                while (alpha > sigma) {
                    cout << "Total correlation distance must be greater than perfect correlation distance!";
                    cout << "Perfect correlation distance? (paddles)";
                    cin >> alpha;
                    cout << "Total correlation distance? (paddles)";
                    cin >> sigma;
                }
                cout << "Height of tail? (between 0 and 1)";
                cin >> height;
                while (height < 0 || height > 1) {
                    cout << "Height of tail must be between 0 and 1!";
                    cout << "Height of tail? (between 0 and 1)";
                    cin >> height;
                }
            } else sigma=0;
            
            
            
            cout << "rms of angles? (0->36 degrees) ";
            cin >> target_rms;
            while ( target_rms < 0 || target_rms > 36) {
                cout << "For SAFETY reasons, rms should be between 0 and 36 degrees!! Try an acceptable value!" <<endl;
                cout << "\n rms of angles? (0->36 degrees) ";
                cin >> target_rms;
            }
            
            alg.max_speed = target_rms;
            alg.max_angle = 5 * target_rms /2;
            alg.min_angle = - 5 *target_rms /2;
            
            
            //cout << "Should the area be kept constant? (1,0) ";
            //cin >> constant;
            constant=0;
            
            cout << "\nPreliminary computations in progress, the grid will move soon.\n";
            
            alg.correlatedMovement(constant, sigma, alpha, height, (int) mode, mrow, mcol, target_rms);
            
        }
        
        else if((int)choice==14){
            float sigma1;
            float sigma2;
            int constant;
            int mode;
            float target_rms1;
            float target_rms2;
            double period;
            double duty_cycle;
            
            cout << "Choose the correlation function: \n 1 - Gaussian \n 2 - 1/r^2 \n 3 - 1/|r|^3 \n 4 - 1/r^4 \n";
            cin >> mode;
            while (mode > 4 || mode < 1){
                cout << "Invalid choice! Try again! \n";
                cout << "\n Choose the correlation function: \n 1 - Gaussian \n 2 - 1/r^2 \n 3 - 1/|r|^3 \n 4 - 1/r^4 \n";
                cin >> mode;
            }
            
            if ((int)mode==1){
                cout << "Sigma1? ";
                cin >> sigma1;
                cout << "Sigma2? ";
                cin >> sigma2;}
            else {sigma1 =0; sigma2 = 0;}
            
            cout << "First (biggest) rms of angles? (0->36 degrees) ";
            cin >> target_rms1;
            while ( target_rms1 < 0 || target_rms1 > 36){
                cout << "For SAFETY reasons, rms should be between 0 and 36 degrees!! Try an acceptable value!" <<endl;
                cout << "\n First (biggest) rms of angles? (0->36 degrees) ";
                cin >> target_rms1;
            }
            
            cout << "Second rms of angles? (0->36 degrees) ";
            cin >> target_rms2;
            while ( target_rms2 < 0 || target_rms2 > 36){
                cout << "For SAFETY reasons, rms should be between 0 and 36 degrees!! Try an acceptable value!" <<endl;
                cout << "\n Second rms of angles? (0->36 degrees) ";
                cin >> target_rms2;
            }
            
            alg.max_speed = max(target_rms1,target_rms2);
            alg.max_angle = 5 * max(target_rms1,target_rms2) /2;
            alg.min_angle = - 5 *max(target_rms1,target_rms2) /2;
            
            cout << "Period of the oscillation of amplitude? (> 0.1s)";
            cin >> period;
            while (period <0.1){
                cout << "Period should be more than 0.1 s! Try an acceptable value!" <<endl;
                cout << "Period of the oscillation of amplitude? (> 0.1s)";
                cin >> period;
            }
            
            cout << "Duty cycle for first state? (between 0 and 1)";
            cin >> duty_cycle;
            while (duty_cycle < 0 || duty_cycle > 1){
                cout << "Duy cycle should be between 0 and 1! Try an acceptable value!" <<endl;
                cout << "Duty cycle for first state? (between 0 and 1)";
                cin >> duty_cycle;
            }
            
            
            cout << "Should the area be kept constant? (1,0) ";
            cin >> constant;
            //constant=0;
            
            cout << "\nPreliminary computations in progress, the grid will move soon.\n";
            
            alg.correlatedMovement_steps(constant,sigma1, sigma2, (int) mode, target_rms1, target_rms2, period, duty_cycle);
        }
        
        
        else if((int)choice==15){
            float sigma;
            int constant;
            int mode;
            float target_rms;
            int numberofsteps;
            
            cout << "Choose the correlation function: \n 1 - Gaussian \n 2 - 1/r^2 \n 3 - 1/|r|^3 \n 4 - 1/r^4 \n";
            cin >> mode;
            while (mode > 4 || mode < 1){
                cout << "Invalid choice! Try again! \n";
                cout << "\n Choose the correlation function: \n 1 - Gaussian \n 2 - 1/r^2 \n 3 - 1/|r|^3 \n 4 - 1/r^4 \n";
                cin >> mode;
            }
            
            if ((int)mode==1){
                cout << "Sigma? ";
                cin >> sigma;
            }
            else sigma =0;
            
            cout << "rms of angles? (0->36 degrees) ";
            cin >> target_rms;
            while ( target_rms < 0 || target_rms > 36){
                cout << "For SAFETY reasons, rms should be between 0 and 36 degrees!! Try an acceptable value!" <<endl;
                cout << "\n rms of angles? (0->36 degrees) ";
                cin >> target_rms;
            }
            
            alg.max_speed = target_rms;
            alg.max_angle = 5 * target_rms /2;
            alg.min_angle = - 5 *target_rms /2;
            
            cout << "Number of steps in the periodic pattern (10 per second)? ";
            cin >> numberofsteps;
            
            cout << "Should the area be kept constant? (1,0) ";
            cin >> constant;
            //constant=0;
            
            cout << "\nPreliminary computations in progress, the grid will move soon.\n";
            
            alg.correlatedMovement_periodic(constant,sigma, (int) mode, target_rms, numberofsteps);
        }
        
        // correlated in space and time routine
        
        // ERROR!!!
        // RIGHT NOW THIS FUNCTION DOES EXACTLY WHAT CHOICE 15 DOES.
        // NEED TO CHANGE IT TO ACTUALLY DO THE CORRELATION IN TIME
        else if((int)choice==16){
            int constantArea = 0;
            float spatial_sigma = 0;
            float temporal_sigma = 0;
	    float alpha = 0;
	    double height;
            int typeOfSpatialCorr = 0;
            int typeOfTemporalCorr = 0;
            float target_rms = 0;
            int numberOfSlices = 0;
            
            cout << "Choose the spatial correlation function: \n"
		" 1 - Gaussian \n 2 - 1/r^2 \n 3 - 1/|r|^3 \n 4 - 1/r^4 \n"
		" 5 - Top hat \n"
		//" 6 - True top hat with a fixed main paddle \n"
		//" 7 - True top hat with a random main paddle \n"
		" 8 - Top hat with a long tail \n"
		" 9 - Triangle\n";
            cin >> typeOfSpatialCorr;
            while (typeOfSpatialCorr > 9 || typeOfSpatialCorr < 1 ||
		   typeOfSpatialCorr == 6 || typeOfSpatialCorr == 7){
                cout << "Invalid choice! Try again! \n";
                cout << "\n Choose the spatial correlation function: \n"
		    " 1 - Gaussian \n 2 - 1/r^2 \n 3 - 1/|r|^3 \n 4 - 1/r^4 \n"
		    " 5 - Top hat \n"
		    //" 6 - True top hat with a fixed main paddle \n"
		    //" 7 - True top hat with a random main paddle \n"
		    " 8 - Top hat with a long tail \n"
		    " 9 - Triangle\n";
		cin >> typeOfSpatialCorr;
            }
            
            if (typeOfSpatialCorr == 1){
                cout << "Spatial Sigma? ";
                cin >> spatial_sigma;
            }
	    if (typeOfSpatialCorr == 5 || typeOfSpatialCorr == 6 || typeOfSpatialCorr == 7 ||
		typeOfSpatialCorr == 8 || typeOfSpatialCorr == 9 ){
	      cout << "Spatial Sigma? "; //Rename?
                cin >> spatial_sigma;
            }
            if (typeOfSpatialCorr == 8){
		cout << "Alpha?\nNote: If you want to use a long tail top hat in both the spatial and temporal dimensions, alpha has to be the same in both dimensions\n"; // Describe what alpha is?
		cin >> alpha;
		cout << "Height?\nNote: If you want to use a long tail top hat in both the spatial and temporal dimensions, height has to be the same in both dimensions\n"; //Describe what height is?
		cin >> height;
	    }
            
            cout << "Choose the temporal correlation function: \n"
		" 1 - Gaussian \n 2 - 1/r^2 \n 3 - 1/|r|^3 \n 4 - 1/r^4 \n"
		" 5 - Top hat \n"
		//" 6 - True top hat with a fixed main paddle \n"
		//" 7 - True top hat with a random main paddle \n"
		" 8 - Top hat with a long tail \n"
		" 9 - Triangle\n";
            cin >> typeOfTemporalCorr;
            while (typeOfTemporalCorr > 9 || typeOfTemporalCorr < 1 ||
		   typeOfTemporalCorr == 6 || typeOfTemporalCorr == 7){
                cout << "Invalid choice! Try again! \n";
                cout << "\n Choose the temporal correlation function: \n"
		    " 1 - Gaussian \n 2 - 1/r^2 \n 3 - 1/|r|^3 \n 4 - 1/r^4 \n"
		    " 5 - Top hat \n"
		    //" 6 - True top hat with a fixed main paddle \n"
		    //" 7 - True top hat with a random main paddle \n"
		    " 8 - Top hat with a long tail \n"
		    " 9 - Triangle\n";
                cin >> typeOfTemporalCorr;
            }
            
            if ((int)typeOfTemporalCorr == 1){
                cout << "Temporal Sigma? ";
                cin >> temporal_sigma;
            }
	    if (typeOfTemporalCorr == 5 || typeOfTemporalCorr == 6 || typeOfTemporalCorr == 7 ||
		typeOfTemporalCorr == 8 || typeOfTemporalCorr == 9 ){
	      cout << "Temporal Sigma? "; //Rename?
	      cin >> temporal_sigma;
            }

	    if (typeOfTemporalCorr == 8){
		cout << "Alpha?\nNote: If you want to use a long tail top hat in both the spatial and temporal dimensions, alpha has to be the same in both dimensions\n"; // Describe what alpha is?
		cin >> alpha;
		cout << "Height?\nNote: If you want to use a long tail top hat in both the spatial and temporal dimensions, height has to be the same in both dimensions\n"; //Describe what height is?
		cin >> height;
	    }
            
            cout << "rms of angles? (0->36 degrees) ";
            cin >> target_rms;
            while (target_rms < 0 || target_rms > 36){
                cout << "For SAFETY reasons, rms should be between 0 and 36 degrees!! Try an acceptable value!" <<endl;
                cout << "\n rms of angles? (0->36 degrees) ";
                cin >> target_rms;
            }
            
            alg.max_speed = target_rms;
            alg.max_angle = 5 * target_rms /2;
            alg.min_angle = - 5 *target_rms /2;
            
            cout << "Should the area be kept constant? (1,0) ";
            cin >> constantArea;
            while (constantArea < 0 || constantArea > 1){
                cout << "Choose zero or 1\n";
                cout << "Should the area be kept constant? (1,0)" << endl;
                cin >> constantArea;
            }
            
                        
	    cout << "What is the range of correlation in the temporal"
	      " dimension?\nIn other words, how many time-steps should the temporal"
	      " kernel encompass?\n(Note: one time-step is equal to five grid positions)" << endl;
	    cin >> numberOfSlices;
	    while (numberOfSlices <= 0){
	      cout << "Choose an integer greater than zero!" << endl;
	      cin >> numberOfSlices;
	    }
            
	    cout << "\nPreliminary computations in progress,"
	      " the grid will move soon.\n" << endl;

	    alg.correlatedMovement_correlatedInTime(constantArea, temporal_sigma, temporal_sigma, alpha, height, typeOfTemporalCorr, typeOfTemporalCorr, target_rms, numberOfSlices);
        } 
        
        /* paddle test routine. Opens and closes each row, one at a time. Then opens and closes each column
         one at a time */
        else if((int)choice==17){
            int NUMBER_OF_ROWS = 11;
            int NUMBER_OF_COLUMNS = 13;
            double openAngle = 0;
            double closedAngle = 90;
	    double closedAngle2 = -90;
            int i; //column counter
            int j; // row counter
            
            alg.grid.opengrid();
            wait(10.0);
            
            for (i = 1; i <= NUMBER_OF_COLUMNS; i++)
            {
                alg.setanglesofonecolumn(i, closedAngle);
                wait(2.0);
                alg.setanglesofonecolumn(i, openAngle);
            }
            for (j = 1; j <= NUMBER_OF_ROWS; j++)
            {
                alg.setanglesofonerow(j, closedAngle2);
                wait(2.0);
                alg.setanglesofonerow(j, openAngle);
            }
        }
        
        else if((int)choice==18){
            int topRow = 11;
            double topAngle = 0;
            int bottomRow = 1;
            double bottomAngle = 0;
            double gridAngle = 0;
            
            cout << "Angle of top row? (-90 - 90)";
            cin >> topAngle;
            cout << endl;
            cout << "Angle of bottom row? (-90 - 90)";
            cin >> bottomAngle;
            cout << endl;
            cout << "Angle of rest of grid? (-90 - 90)";
            cin >> gridAngle;
            topAngle = angle_check(topAngle);
            bottomAngle = angle_check(bottomAngle);
            gridAngle = angle_check(gridAngle);
            alg.setsameangletoallservos(gridAngle);
            alg.setanglesofonerow(topRow,topAngle);
            alg.setanglesofonerow(bottomRow,bottomAngle);   
        }

	else if((int)choice==19){
	  testloaf();
	  i = 0;
	}
        
        else if((int)choice>=20 || choice<=0) i=0;
    }
    return 0;
    
}

double angle_check (double angle) {
    if (fabs(angle)>90) {
        cout << "Maximum angle is +-90 degrees!! Angle set to +-90 degrees.";
        if (angle>0) return 90;
        else return -90;
    }
    else return angle;
}

double check_frequency(double frequency, double amplitude) {
    double angleperstep;
    if (frequency<=0) {
        cout << "Frequency should be greater than 0 Hz." << endl;
        return -1;
    }
    angleperstep=2*fabs(amplitude)*frequency/5;
    if (angleperstep>42.8) {
        cout << "Chosen frequency mismatches amplitude, maximum servo speed is exceeded. Frequency is adjusted to: ";
        frequency=5*42.8/(2*amplitude);
        cout << frequency << "Hz.";
        return frequency;
    }
    else return frequency;
}
