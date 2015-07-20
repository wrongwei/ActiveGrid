/*------------------------------------------------------------------------*/
/* Library to create sets of angles and angular speeds, in order to give  */
/* to the grid a specific motion.                                         */
/* There are a lot of functions, but also a lot of possible motions.      */
/*                                                                        */
/* Florent Lachaussée
 * florent.lachaussee@ens.fr
 * March 2013                                                             */
/* algo.h                                                                 */
/*------------------------------------------------------------------------*/

// Protecting against multiple inclusion ...
#ifndef ALGO_INCLUDED
#define ALGO_INCLUDED

/*------------------------------------------------------------------------*/

#include <iostream>
#include <vector>
#include <cstdlib>
#include <time.h>
#include <cmath>
#include <fstream>
#include "activegrid.h"
#include "loaf.h"

#define range_of_corr 7

class algo{
public:
    algo(); //constructor
    ~algo(){}; //destructor
    
    // create an active grid object
    activegrid grid;
    
    // functions for motion (basic ones)
    int  setsameangletoallservos(double angle);
    int  setsamespeedtoallservos(double angleperstep);
    double checkorientationofservo(int col, int row, double angle);
    int  setangleofoneservo(int col, int row, double angle);
    int  setanglesofonerow(int tunedrow, double angle);
    int  setanglesofonecolumn(int tunedcol, double angle);
    int  openhalfclosehalf();
    
    // principle functions to pass the angle orders
    int setanglestoallservos(float * positions, float * anglesteps, int combine, int constant, float rms);
    int setanglestoallservosII(float * positions, float * anglesteps, int constant, float rms);// light version
    int setanglestoallservosIII(float angles[13][11], float steps[13][11], int constant, float rms);// 3d version
    
    // functions computing "complicated" motions
    int allperiodic(double angle, double frequency);
    int chaoticMovement(int combine, int constant, int option);
    int correlatedMovement (int constant, float sigma, float alpha, double height, int mode, int mrow, int mcol, float target_rms);
    int correlatedMovement_steps (int constant, float sigma1, float sigma2, int mode, float target_rms1, float target_rms2, double period, double duty_cycle);
    int correlatedMovement_periodic (int constant, float sigma, int mode, float target_rms, int numberofsteps);
    int correlatedMovement_correlatedInTime(int constantArea, float spatial_sigma, float temporal_sigma, float alpha, double height, int typeOfSpatialCorr, int typeOfTemporalCorr, float target_rms, int numberOfSlices);
    
    //keeps the projected area constant
    void area(double actpos[14][12], float rms);
    double angle_constant_area;
    
    
    // output file for computed angles
    ofstream anglefile;
    ofstream angleperiod;
    
    // computation of convolution used for corrleated motion is easier with arrays,
    // so we reintroduce temporarily array notations. Conversion tools follow.
    int columns[numberOfServos];
    int rows[numberOfServos];
    
    
    //useful quantities for computations
    
    float min_speed,max_speed;
    float min_amplitude,max_amplitude;
    double min_angle, max_angle;
    float speed;
    float frequency;
    float amplitude;
    int updatetimeinmus;//in mu sec
    double * old_angle;
    double * new_angle;
    float norm, norm1, norm2;
    
    //methods (used for chaotic and correlated movement)
    void updateOneWing(int WingNumber);
    void updateOneWing2(int WingNumber);
    void run(float actpos[], float actstep[], int option);
    void runcorr(float actpos[], float actstep[], float sigma, float alpha,
		 double height, int mode, int mrow, int mcol, float correction,
		 float norm, float oldpos[], float oldstep[], float err[]);
    void runcorr_3D(float newslice[][11], loaf& myLoaf, int halfLoaf, int upperTimeBound, float spaceSigma, float timeSigma, float alpha,
                    double height, int spaceMode, int timeMode, int mrow, int mcol, float correction);
    float compute_rms(int option);
    float compute_rmscorr(float sigma, int mode, float alpha, double height,
			  int mrow, int mcol);
    float compute_rmscorr_3D(float spaceSigma, float timeSigma, int spaceMode, int timeMode, float alpha, double height,
                             int mrow, int mcol, int halfLoaf, int upperTimeBound);
    
    // helper methods for correlation procedures
    void initialize_pos_step(float actpos[], float actstep[], float oldpos[], float oldstep[], int i);
    void gaussian2d(float actpos[], float oldpos[], float actstep[], float sigma, float correction,
                    float norm, float interpos[][11], float err[], int i);
    void inverse_r_to_n_2d(float actpos[], float oldpos[], float actstep[], float correction,
                           float interpos[][11], int n, int i);
    void tophat2d(float actpos[], float oldpos[], float actstep[], float correction,
                  float interpos[][11], float sigma, int i);
    void truetophat2d(float actpos[], float oldpos[], float actstep[], float correction,
                      float interpos[][11], float sigma, int mrow, int mcol, int i);
    void tophatlongtail2d(float actpos[], float oldpos[], float actstep[], float correction,
                          float interpos[][11], float sigma, float alpha, double height, int i);
    void triangle2d(float actpos[], float oldpos[], float actstep[], float correction,
                    float interpos[][11], float sigma, int i);
    void safety_check(float actpos[], float actstep[], float err[], bool isGaussian, int i);
    
    // to store speeds and positions after random computation and before convolution
    vector<float> * positions_random;
    vector<float> * steps_random;
    int * actualpositioninvector;
    
    // to store speeds and angles to do the convolution and to send orders
    float positions[numberOfServos];
    float anglesteps[numberOfServos];
    float old_positions[numberOfServos];
    float old_steps[numberOfServos];
    float err[numberOfServos];
    
    // This function performs the modulo operation but will always have NON NEGATIVE remainder
    int modulo(int numerator, int denominator);
    
    
private:
    
};

inline algo::algo(){
    srand(time(NULL));//initialize random seed
    cout << "algo::algo()" << endl;
    
    cout << "numberOfServos: " << numberOfServos << endl;
    positions_random = new vector<float>[numberOfServos]; // pointer to (ie array of) vectors
    steps_random = new vector<float>[numberOfServos];
    actualpositioninvector = new int[numberOfServos];
    
    old_angle = new double [numberOfServos];
    new_angle= new double [numberOfServos];
    
    //get speed,amplitude:
    //maximal possible angle-per-step-motion (for updaterate 10Hz): 42.8 degrees/step, limited by the servo speed
    min_speed=10.;
    max_speed=40.0;
    // for piecewise periodic motion
    min_amplitude=30;
    max_amplitude=90; 
    // max and min angle for random motion
    // note that the loaf data type has its own values for min and max angle so
    // changing the values below will not affect the values in loaf.cpp
    min_angle=-50;
    max_angle=50;
    
    updatetimeinmus = 100000;//in mu sec !!!!!!!! SHOULD BE 100 000 !!!!!!!!!!!!!
    
    // initializations
    for(int i=0;i<numberOfServos;i++){
        old_angle[i]=0;
        new_angle[i]=0;
        
        positions[i]=0;
        anglesteps[i]=0;
        old_positions[i]=0;
        old_steps[i]=0;
        err[i]=0;
        
        positions_random[i].clear();
        steps_random[i].clear();
        positions_random[i].push_back(-100);
        steps_random[i].push_back(-100);
        actualpositioninvector[i]=0;
        
        columns[i]=i%13;// !! indexed by integers from 0 to 12 !!
        rows[i]= floor(i/13);
    }
    grid.high_duty = false;
    
    
    cout << "algo::algo() END" << endl;
}

#endif // ALGO_INCLUDED
