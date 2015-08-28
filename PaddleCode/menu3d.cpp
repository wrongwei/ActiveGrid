
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
 
 Dependencies: algo.h, algo3d.h (which rely on pickCorrelations.h, loaf.h,
   activegrid.h, SD84.h, and SD84constants.h, and all the associated .cpp files)
 
 */
/*--------------------------------------------------------------------*/

#include <iostream>
#include <sstream>
#include <unistd.h>
#include <time.h>
#include "lib/algo3d.h"

/*--------------------------------------------------------------------*/

using namespace std;

//functions:
void wait(float seconds);
int main(int argc, char * const argv[]);
double angle_check(double angle);
double check_frequency(double frequency, double amplitude);
activegrid grid;

void wait(float seconds) {
    clock_t endwait;
    endwait = clock () + seconds * CLOCKS_PER_SEC ;
    while (clock() < endwait) {}
}

int main (int argc , char * const argv[]) {
    // create an algo object
    //algo alg;
    int choice=-1;
    int i=1;
        
    // Welcome
    cout << endl << "Welcome to " << argv[0] << endl;
    cout << "This is the program to run the software of the activegrid." << endl;
    //activegrid grid;    
    //    grid.opengrid();
    //grid.high_duty = false;
    //    cout << grid.high_duty << endl;
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
        " 10 - correlated in space and correlated in time\n\n"
        " 11 - grid test protocol\n\n"
        " 12 - set boundary paddles to constant angle\n\n"
        " 13 - test loaf object\n\n"
        " 14 - end program\n\n";
        cin >> choice;
        cout << "your choice is  " << choice << "\n";

        if((int)choice==1){
            grid.showallangles();
        }
        
        else if((int)choice==2){
	    grid.showservomap();
        }
        
        else if((int)choice==3){
	    grid.opengrid();
	}
        /* need to make a method in algo3d that takes an algo object as a parameter
        else if((int)choice==4){
            double angle=0;
            cout << "Angle? (-90 - 90)";
            cin >> angle;
            angle=angle_check(angle);
            cout << endl;
            alg.setsameangletoallservos(angle);
        }
        */

        /* need to make a method in algo3d that takes an algo object as a parameter
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
        */

        /* need to make a method in algo3d that takes an algo object as a parameter
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
	*/
        
	/* need to make a method in algo3d that takes an algo object as a parameter
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
        */

	/* need to make a method in algo3d that takes an algo object as a parameter
	   else if((int)choice==8){
            alg.openhalfclosehalf();
        }
        */

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
            if(grid.getangle(col,row)!= Fictive )
                cout << "Angle of (" << col << "," << row << "): " << grid.getangle(col,row) << endl;
            else
                cout << "(" << col << "," << row << ") does not exist!" << endl;
        }
        
        // correlated in space and time routine
        else if((int)choice==10){
            int constantArea = 0;
            float spatial_sigma = 0;
            float temporal_sigma = 0;
            float spatial_alpha = 0;
            float temporal_alpha = 0;
            float spatial_height = 0;
            float temporal_height = 0;
            int typeOfSpatialCorr = 0;
            int typeOfTemporalCorr = 0;
            float target_rms = 0;
            int numberOfSlices = 0;
            
            // CorrParameters temporalCorrParam;
            // CorrParameters spatialCorrParam;
            
            cout << "Choose the spatial correlation function: \n"
            " 0 - Random \n 1 - Gaussian \n 2 - 1/r^2 \n 3 - 1/|r|^3 \n 4 - 1/r^4 \n"
            " 5 - Top hat \n"
            //" 6 - True top hat with a fixed main paddle \n"
            //" 7 - True top hat with a random main paddle \n"
            " 8 - Top hat with a long tail \n"
            " 9 - Triangle\n"
            " 10 - Unsharp filter\n";
            cin >> typeOfSpatialCorr;
            while (typeOfSpatialCorr > 10 || typeOfSpatialCorr < 0 ||
                   typeOfSpatialCorr == 6 || typeOfSpatialCorr == 7){
                cout << "Invalid choice! Try again! \n";
                cout << "\n Choose the spatial correlation function: \n"
                " 0 - Random \n 1 - Gaussian \n 2 - 1/r^2 \n 3 - 1/|r|^3 \n 4 - 1/r^4 \n"
                " 5 - Top hat \n"
                //" 6 - True top hat with a fixed main paddle \n"
                //" 7 - True top hat with a random main paddle \n"
                " 8 - Top hat with a long tail \n"
                " 9 - Triangle\n"
                " 10 - Unsharp filter\n";
                cin >> typeOfSpatialCorr;
            }
            
            if (typeOfSpatialCorr == 1 || typeOfSpatialCorr > 4){
                cout << "Spatial Sigma? (0->7) "; // maximum length of 7 because this is the biggest sigma loaf can handle at the moment
                cin >> spatial_sigma;
                while (spatial_sigma < 0 || spatial_sigma > 7) {
                    cout << "Spatial sigma must be between 0 and 7. Try an acceptable value!" <<endl;
                    cout << "\nSpatial Sigma? (0->7) ";
                    cin >> spatial_sigma;
                }
            }
            if (typeOfSpatialCorr == 8){
                cout << "Spatial Alpha? This is the width of part with correlation of 1 (usually spatial_alpha=0)\n";
                cin >> spatial_alpha;
                cout << "Spatial Height? This is how tall the tail is. 1 would be a correlation of 1 and 0 would be no correlation in the tail.\n";
                cin >> spatial_height;
            }
            else if (typeOfSpatialCorr == 10) {
                cout << "Spatial Alpha? (This sets the radius of the positive portion of the unsharp kernel)\n";
                cin >> spatial_alpha;
                cout << "Spatial Depth? (This sets the magnitude of the negative portion of the unsharp kernel, should be between 0 and 1)\n";
                cin >> spatial_height;
            }
            
            cout << "Choose the temporal correlation function: \n"
            " 0 - Random \n 1 - Gaussian \n 2 - 1/r^2 \n 3 - 1/|r|^3 \n 4 - 1/r^4 \n"
            " 5 - Top hat \n"
            //" 6 - True top hat with a fixed main paddle \n"
            //" 7 - True top hat with a random main paddle \n"
            " 8 - Top hat with a long tail \n"
            " 9 - Triangle\n"
            " 10 - Unsharp filter\n";
            cin >> typeOfTemporalCorr;
            while (typeOfTemporalCorr > 10 || typeOfTemporalCorr < 0 ||
                   typeOfTemporalCorr == 6 || typeOfTemporalCorr == 7){
                cout << "Invalid choice! Try again! \n";
                cout << "\n Choose the temporal correlation function: \n"
                " 0 - Random \n 1 - Gaussian \n 2 - 1/r^2 \n 3 - 1/|r|^3 \n 4 - 1/r^4 \n"
                " 5 - Top hat \n"
                //" 6 - True top hat with a fixed main paddle \n"
                //" 7 - True top hat with a random main paddle \n"
                " 8 - Top hat with a long tail \n"
                " 9 - Triangle\n"
                " 10 - Unsharp filter\n";
                cin >> typeOfTemporalCorr;
            }
            
            if ((int)typeOfTemporalCorr == 1 || (int)typeOfTemporalCorr == 10){
                cout << "Temporal Sigma? (has units of 0.1 seconds times spacing) ";
                cin >> temporal_sigma;
            }
            if (typeOfTemporalCorr == 5 || typeOfTemporalCorr == 6 || typeOfTemporalCorr == 7 ||
                typeOfTemporalCorr == 8 || typeOfTemporalCorr == 9 ){
                cout << "Temporal Sigma? (has units of 0.1 seconds times spacing) ";
                cin >> temporal_sigma;
            }
            
            if (typeOfTemporalCorr == 8){
                cout << "Temporal Alpha? This is the width of the part with correlation of 1 (usually temporal_alpha=0)\n";
                cin >> temporal_alpha;
                cout << "Temporal Height? This is how tall the tail is. 1 would be a correlation of 1 and 0 would be no correlation in the tail.\n";
                cin >> temporal_height;
            }
            
            else if (typeOfTemporalCorr == 10) {
                cout << "Temporal Alpha? (This sets the radius of the positive portion of the unsharp kernel)\n";
                cin >> temporal_alpha;
                cout << "Temporal Depth? (This sets the magnitude of the negative portion of the unsharp kernel, should be between 0 and 1)\n";
                cin >> temporal_height;
            }
            
            cout << "rms of angles? (0->50 degrees) "; // 30 * sqrt(2) = 42.4, which is about the maximum servo speed (42.8 degrees per 0.1 second interval)
            cin >> target_rms;
            while (target_rms < 0 || target_rms > 50){
                cout << "For SAFETY reasons, rms should be between 0 and 50 degrees!! Try an acceptable value!" <<endl;
                cout << "\n rms of angles? (0->50 degrees) ";
                cin >> target_rms;
            }
            /* NOT IMPLEMENTED - CAUSES PERIODIC FREEZING ISSUES IN THE GRID
            cout << "Should the area be kept constant? (1,0) ";
            cin >> constantArea;
            while (constantArea < 0 || constantArea > 1){
                cout << "Choose zero or 1\n";
                cout << "Should the area be kept constant? (1,0)" << endl;
                cin >> constantArea;
            }*/
            
            // calculate width of temporal kernel (number of time-slices to analyze at a time)
            if ((int)typeOfTemporalCorr == 1) { // calculate width of Gaussian kernel
                numberOfSlices = ceil(6 * temporal_sigma) + 1; // odd number with 3 standard deviations on either side of the middle slice
            }
            else if ((int)typeOfTemporalCorr > 1 && (int)typeOfTemporalCorr <= 4) { // user sets width of 1/r^n kernel manually
                cout << "How many time-slices should the temporal kernel include? (number must be odd) ";
                cin >> numberOfSlices;
                while (numberOfSlices <= 0 || numberOfSlices % 2 == 0) { // input must be odd and > 0 so we can define a middle slice and have nice symmetry
                    cout << "Input must be an odd integer greater than zero!" << endl;
                    cin >> numberOfSlices;
                }
            }
            else { // top hats, triangle, etc. - return values of 0 outside of one standard deviation
                numberOfSlices = (2 * temporal_sigma) + 1; // so kernel just needs to have 1 standard deviation on either side of the middle slice
            }
            if (numberOfSlices % 2 == 0) numberOfSlices++; // make sure numberOfSlices is odd
            cout << "\nThe range of the temporal correlation is " << numberOfSlices << " time-steps.\nNote: one time-step is equal to 'SPACING' grid positions, where the value of SPACING is set manually in algo.cpp. \nEach grid position lasts 0.1 seconds." << endl;
            
            cout << "\nInitializing preliminary computations...\n" << endl;
            
	    correlatedMovement_correlatedInTime(constantArea, spatial_sigma, temporal_sigma, spatial_alpha, temporal_alpha, spatial_height, temporal_height, typeOfSpatialCorr, typeOfTemporalCorr, target_rms, numberOfSlices);
        }
        
        /* paddle test routine. Opens and closes each row, one at a time. Then opens and closes each column
         one at a time */
        
	/* need to write algo3d setanaglesofonecolumn first
	else if((int)choice==11){
            int NUMBER_OF_ROWS = 11;
            int NUMBER_OF_COLUMNS = 13;
            double openAngle = 0;
            double closedAngle = 90;
            double closedAngle2 = -90;
            int i; //column counter
            int j; // row counter
            
            grid.opengrid();
            cout << "\nThe test code will begin in 10 seconds. Quick, run over and watch the grid!" << endl;
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
	*/
	
	/* need to write setanglesofonerow in algo3d
        else if((int)choice==12){
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
	*/
        
        else if((int)choice==13){
            testloaf();
            i = 0;
        }
        
        else if((int)choice>=14 || choice<=0) i=0;
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
