/***************************** ALGO3D ******************************
 
 Methods and functions for temporal correlation control code
 These functions were adapted from the old algo.cpp / algo.h, which was
 written before 2015 and debugged/refactored/optimized July-August 2015.
 
 These methods are entirely independent of algo.h and menuII, since their
 functions are different enough to merit a separate, cleaner executable.
 
 These methods are accessed by menu3d.cpp, which handles routing depending on
 user inputs. They depend on pickCorrelations.cpp/pickCorrelations.h,
 which contain function pointers to the different correlation kernels.
 
 Note that this file contains multiple versions of runcorr_3D. This is because
 our computers were not fast enough to run the function pointer implementation
 for very large temporal kernels (sigma = 50). This workaround should ONLY be
 used for specific kernels where speed is an issue; otherwise, the standard
 function pointer runcorr_3D should be used for readability and simplicity.
 This also applies for the unsharp kernel, which requires special absolute value
 correlation handling when combined with itself in order to act as expected.
 
 Dependencies: pickCorrelations.h/.cpp, activegrid.h/.cpp, loaf.h/.cpp
    (only first order dependencies are listed)
 
 Written by Kevin Griffin and Nathan Wei, July-August 2015
 
 ********************************************************************/

#include <algo3d.h>
#include <assert.h>
#include <time.h>
#include <curses.h>

/*--------------------------------------------------------------------*/

//Headers needed for signal handling
#define _GNU_SOURCE
#include <stdio.h>
#include <stdlib.h>
#include <signal.h>

/*--------------------------------------------------------------------*/
// Install a signal handler to make the program print out some useful
// statistics before terminating.
// If you want to terminate the program without calling the signal handler,
// type CTRL-backslash to exit immediately

// private global variables used by the signal handler to keep track of out of bounds paddles
char * myLoaf = NULL; //necessary?????
static int outOfBoundsCount = 0; // number of paddles that try to move more than 40 degrees per time step
static int over90orminus90count = 0; // numer of paddles that try to move to angles > 90 or angles < -90
static int numberOfAnglesSet = 0; // total number of angles set. Allows you to calcuate percent of paddles out of bounds

// install the signal handler
static void mySignalHandler(int iSignal){
    cout << "Percent of paddles that wanted to move more than 40 degrees in .1 seconds but were constrained: " << ((float)outOfBoundsCount / numberOfAnglesSet) * 100 << "%" << endl;
    cout << "Percent of angles that wanted to be > 90 or < -90 but were clipped: " << ((float)over90orminus90count / numberOfAnglesSet) * 100 << "%" << endl;
    cout << "\nRemember to save the angle file!" << endl;
    // Kill the program
    cout << "Exiting the program.\n";
    fflush(NULL);
    exit(0);
}

/*--------------------------------------------------------------------*/

// global variables needed for algo3d method operation
// (in place of algo.h variable storage implementation in the 2D version)
extern activegrid grid;
int updatetimeinmus = 100000; // time to calculate correlation and send angles to grid = 0.1 seconds
float max_speed = 40.0; // maximum speed of servos = 42.8 degrees / 0.1 seconds
ofstream anglefile; // file to write angles to (for MATLAB processing, if desired)
float norm; // normalization for correlation functions

// set differents angles to all servos (angles array provided by runcorr_3D/correlatedMovement_correlatedInTime)
static int setanglestoallservosIII(float angles[13][11], float steps[13][11], int constant, float rms){
    double newangle[14][12];
    double newangleperstep[14][12];
    for(int row = 1; row < 12; row++) {
        for(int col = 1; col < 14; col++) {
            if (grid.servo[col][row]!=0){
                // converting 13x11 array into 14x12 array (necessary for servo indexing in activegrid.cpp)
                newangle[col][row]=angles[(col-1)][(row-1)];
                newangleperstep[col][row]=steps[(col-1)][(row-1)];
            }
            else { // for non-existing servos
                if (row!=2) {
                    newangle[col][row]=Fictive;
                    newangleperstep[col][row]=0;
                }
            }
        }
    }
    
    //if(constant==1) {area(newangle, rms);} // NOT IMPLEMENTED because it caused periodic freezing of the grid
    
    // set speeds and angles to grid (methods in activegrid.cpp)
    grid.setspeeds(newangleperstep);
    grid.setanglesII(newangle);
    
    // write all angles to file
    for (int j = 0; j < 11; j++) {
        for (int i = 0; i < 13; i++)
            anglefile << "    " << angles[i][j];
    }
    anglefile << endl;
    
    return 1;
}

// determines the positions of the servos by doing a 3D correlation. Stores these
// positions in a 2D array, which lives in correlatedMovement_correlatedInTime.
// Unlike its predecessor, this method does not deal with step-setting. That is done entirely by the client that calls it.
static void runcorr_3D(float newslice[][11], loaf* myLoaf, int halfLoaf, float spaceSigma, float timeSigma, float spaceAlpha, float timeAlpha, float spaceHeight,
                float timeHeight, int spaceMode, int timeMode, int mrow, int mcol, float correction) {
    
    // convolution to create correlation between paddles
    // periodic boundary conditions are used. For explanation of auto/cross-correlation see wikipedia
    
    float crumb = 0;
    float bound;
    if (spaceMode <= 4) bound = range_of_corr;
    else bound = spaceSigma;
    
    // Declare function pointers for the spatial and temporal correlation functions
    float (*pfSpatialCorr)(int j, int k, float spaceSigma, float spaceAlpha, float height);
    float (*pfTemporalCorr)(int t, float timeSigma, float timeAlpha, float height);
    pfSpatialCorr = pickSpatialCorr(spaceMode);
    pfTemporalCorr = pickTemporalCorr(timeMode);
    
    // Loop through servos and calculate/create correlations, using helper methods
    for (int row = 0; row < 11; row++) {
        for (int col = 0; col < 13; col++){
            newslice[col][row] = 0; // start each angle at zero, then add in results of correlation
            for (int j = -bound; j <= bound; j++) { // range of neighbours used to compute convolution
                for (int k = -bound; k <= bound; k++) { // j and k refer to the shift
                    for (int t = -halfLoaf; t <= halfLoaf; t++) { // t taken from the center of the loaf
                        crumb = myLoaf->Loaf_access(j + col, k + row, t + halfLoaf);
                        /* This code is for absolute value correlations (i.e. 47 degrees is seen as perfectly correlated with -47 degrees)
                         if (crumb < 0)
                         crumb = -crumb;
                         */
                        //cout << (pfTemporalCorr(t, timeSigma, height) * pfSpatialCorr(j, k, spaceSigma, height)) << endl; // debugging
                        // multiply original angle by correction factor, spatial correlation function, and temporal correlation function,
                        // which eventually gives the convolution product of kernel matrix with random angles matrix
                        newslice[col][row] += (correction * crumb * pfSpatialCorr(j, k, spaceSigma, spaceAlpha, spaceHeight) * pfTemporalCorr(t, timeSigma, timeAlpha, timeHeight));
                    }
                }
            }
            /* This code is for absolute value correlations
             crumb = myLoaf->Loaf_access(col, row, halfLoaf);
             if (crumb < 0)
             newslice[col][row] = -newslice[col][row];
             */
            // normalization by coefficient calculated in correlatedMovement_correlatedInTime
            newslice[col][row] = newslice[col][row] / norm;
        }
    }
    myLoaf->Loaf_slice(); // remove oldest slice and add new slice - now we're ready for another iteration!
}

// Faster version of runcorr_3D, for long tail in both directions
// Necessary to avoid timing issues when running tests on lt5.2lt50
static void ltfast(float newslice[][11], loaf* myLoaf, int halfLoaf, float spaceSigma, float timeSigma, float spaceAlpha, float timeAlpha, float spaceHeight, float timeHeight, float correction) {
    
    // convolution to create correlation between paddles
    // periodic boundary conditions are used
    float crumb;
    float dist;
    float abs_t;
    float spatialCorrFactor;
    float temporalCorrFactor;
    
    // Loop through servos and calculate/create correlations, using inlined longtail/longtail function
    // runs faster than accessing function pointers, but is not sustainable programmatically
    for (int row = 0; row < 11; row++) {
        for (int col = 0; col < 13; col++){
            newslice[col][row] = 0; // start each angle at zero, then add in results of correlation
            for (int j = -spaceSigma; j <= spaceSigma; j++) { // range of neighbours used to compute convolution
                for (int k = -spaceSigma; k <= spaceSigma; k++) { // j and k refer to the shift
                    for (int t = -halfLoaf; t <= halfLoaf; t++) { // t taken from the center of the loaf
                        crumb = myLoaf->Loaf_access(j + col, k + row, t + halfLoaf);
                        dist = sqrt((j*j)+(k*k));
                        //spatial correlation inlined for top hat long tail
                        if (dist <= spaceAlpha)
                            spatialCorrFactor = 1;
                        else if (dist <= spaceSigma)
                            spatialCorrFactor = spaceHeight;
                        else spatialCorrFactor = 0;
                        //temporal correlation inlined for top hat long tail
                        abs_t = abs(t);
                        if (abs_t <= timeAlpha)
                            temporalCorrFactor = 1;
                        else if (abs_t <= timeSigma)
                            temporalCorrFactor = timeHeight;
                        else
                            temporalCorrFactor = 0;
                        // calculate correlation (convolution product of kernel matrix with random angles matrix)
                        newslice[col][row] += (correction * crumb * spatialCorrFactor * temporalCorrFactor);
                    }
                }
            }
            // normalization by coefficient calculated in correlatedMovement_correlatedInTime
            newslice[col][row] = newslice[col][row] / norm;
        }
    }
    myLoaf->Loaf_slice(); // remove oldest slice and add new slice
}

// Modified version of runcorr_3D, for unsharp in both directions
// Necessary to avoid issue with 2 unsharps having negative signs that multiply to become positive
// Also, in order for unsharp to work (i.e. not look like top hat long tail), it needs the absolute value definition of correlation
static void unsharp(float newslice[][11], loaf* myLoaf, int halfLoaf, float spaceSigma, float timeSigma, float spaceAlpha, float timeAlpha, float spaceDepth, float timeDepth, float correction) {
    
    // convolution to create correlation between paddles
    // periodic boundary conditions are used
    float crumb;
    float dist;
    float abs_t;
    float spatialCorrFactor;
    float temporalCorrFactor;
    
    // Loop through servos and calculate/create correlations, using hard-coded double unsharp with angle sign correction
    for (int row = 0; row < 11; row++) {
        for (int col = 0; col < 13; col++){
            newslice[col][row] = 0; // start each angle at zero, then add in results of correlation
            for (int j = -spaceSigma; j <= spaceSigma; j++) { // range of neighbours used to compute convolution
                for (int k = -spaceSigma; k <= spaceSigma; k++) { // j and k refer to the shift
                    for (int t = -halfLoaf; t <= halfLoaf; t++) { // t taken from the center of the loaf
                        crumb = myLoaf->Loaf_access(j + col, k + row, t + halfLoaf);
                        // Implement absolute value correlations (i.e. 47 degrees is seen as perfectly correlated with -47 degrees)
                        if (crumb < 0)
                            crumb = -crumb; // make everything positive before correlating
                        dist = sqrt((j*j)+(k*k));
                        //spatial correlation inlined for top hat long tail
                        if (dist <= spaceAlpha)
                            spatialCorrFactor = 1;
                        else if (dist <= spaceSigma)
                            spatialCorrFactor = -spaceDepth;
                        else spatialCorrFactor = 0;
                        //temporal correlation inlined for top hat long tail
                        abs_t = abs(t);
                        if (abs_t <= timeAlpha)
                            temporalCorrFactor = 1;
                        else if (abs_t <= timeSigma) {
                            temporalCorrFactor = -timeDepth;
                            if (dist > spaceAlpha) // if outside both alphas but inside both sigmas, then flip sign
                                crumb = -crumb; // since the 2 negative kernels would otherwise multiply to give + sign, not what we want
                        }
                        else
                            temporalCorrFactor = 0;
                        // calculate correlation (convolution product of kernel matrix with random angles matrix)
                        newslice[col][row] += (correction * crumb * spatialCorrFactor * temporalCorrFactor);
                    }
                }
            }
            // Absolute value correlations: after correlating on positive angles, return negative signs where necessary
            crumb = myLoaf->Loaf_access(col, row, halfLoaf);
            if (crumb < 0)
                newslice[col][row] = -newslice[col][row]; // if angle was originally negative, make it so again
            
            // normalization by coefficient calculated in correlatedMovement_correlatedInTime
            newslice[col][row] = newslice[col][row] / norm;
        }
    }
    myLoaf->Loaf_slice(); // remove oldest slice and add new slice
}

/* takes a random 3D sequence and computes its std dev (=rms). It helps calculate the correction
 coefficent that is needed to give to the output the desired rms value of angles. */
static float compute_rmscorr_3D(float spaceSigma, float timeSigma, int spaceMode, int timeMode, float spaceAlpha, float timeAlpha, float spaceHeight, float timeHeight, int mrow, int mcol, int halfLoaf) {
    cout << "compute_rmscorr_3D is running tests now" << endl << "Countdown:" << endl;
    // set up test parameters
    float mean = 0;
    float rms = 0;
    int trials = 4000;
    float slice[13][11] = {{0}}; // this array begins life stuffed with zeros
    float slicestorage[13][11][trials];
    loaf testLoaf = loaf(halfLoaf*2 + 1); // bake test loaf of width = numberOfSlices (recomposed from halfLoaf)
    
    // takes a random correlated sequence of angles, without correction, and executes 4000 sample runs of runcorr_3D
    for (int t = 0; t < trials; t++) {
        runcorr_3D(slice, &testLoaf, halfLoaf, spaceSigma, timeSigma, spaceAlpha, timeAlpha, spaceHeight, timeHeight, spaceMode, timeMode, mrow, mcol, 1);
        if (t % 100 == 0) cout << (4000 - t) / 100 << endl; // countdown to finish (UI)
        for (int row = 0; row < 11; row++) {
            for (int col = 0; col < 13; col++){
                mean += slice[col][row] / (numberOfServos * trials); // calculate mean as we go
                slicestorage[col][row][t] = slice[col][row]; // store angle values for future use in rms calculation
            }
        }
    }
    // calculate variance from previously-found mean and angle measurements
    for (int t = 0; t < trials; t++) {
        for (int row = 0; row < 11; row++) {
            for (int col = 0; col < 13; col++){
                rms += pow(slicestorage[col][row][t] - mean, (int) 2) / (numberOfServos * trials);
            }
        }
    }
    rms = sqrt(rms); // rms is the sqrt of variance
    
    return rms;
}

/* Movement of the paddles that is correlated in space and in time
 * This method currently operates on a time-schedule of 1 new correlation computation to the grid every 0.1 seconds.
 * In principle, this can be easily modified to separate the time-scale of the grid from the time-scale of the correlations
 *  (e.g. by running one correlation every 0.5 seconds but still sending angles to the grid every 0.1 seconds). The current schedule
 *  gives the highest rate of correlation computations we can handle, given the nature of the runcorr_3D algorithm and our current computers.
 * We define the term "constraining" as the occasion when a paddle is given a new angle to move to that it cannot reach without exceeding the
 *  limits of its servo (e.g. an amplitude of 60 degrees when the maximum angular distance it can travel in 0.1 seconds is 42.8 degrees). Such
 *  noncompliant paddles will be "constrained" to move at the greatest allowed servo speed (e.g. 4.28 deg/s instead of 6 deg/s) and thus will
 *  not reach their target destinations. This was a necessary comprimise in order to maximize the grid's correlation computation rate.
 */
int correlatedMovement_correlatedInTime(int constantArea, float spatial_sigma, float temporal_sigma, float spatial_alpha, float temporal_alpha, float spatial_height, float temporal_height, int typeOfSpatialCorr, int typeOfTemporalCorr, float target_rms, int numberOfSlices){

    /* To install mySignalHandler (which prints statistics when you type CTRL-C)
    you must make sure that CTRL-C signals are not blocked.
    If you are having problems with the signal handler (i.e. CTRL-C signals are blocked)
    you can just comment out the following block and the block after this one where the signal handler is actually installed.*/
    void (*pfRet)(int);
    sigset_t sSet;
    int iRet;
    iRet = sigemptyset(&sSet);
    if (iRet == -1) {perror("ERROR installing signal handler in correlatedMovement_correlatedInTime"); exit(EXIT_FAILURE); }
    iRet = sigaddset(&sSet, SIGINT);
    if (iRet == -1) {perror("ERROR installing signal handler in correlatedMovement_correlatedInTime"); exit(EXIT_FAILURE); }
    iRet = sigprocmask(SIG_UNBLOCK, &sSet, NULL);
    if (iRet == -1) {perror("ERROR installing signal handler in correlatedMovement_correlatedInTime"); exit(EXIT_FAILURE); }
    
    // Install mySignalHandler as the handler for CTRL-C signals.
    pfRet = signal(SIGINT, mySignalHandler);
    if (pfRet == SIG_ERR) {perror("error installing signal handler in correlatedMovement_correlatedInTime"); exit(EXIT_FAILURE); }

    // special method selection (for slow or unconventional kernels that need to be inlined)
    bool ltfast_on = false;
    bool unsharp_on = false;
    if (typeOfSpatialCorr == 8 && typeOfTemporalCorr == 8) {
        cout << "\n***********************************\n";
        cout << "*      Running ltfast method      *";
        cout << "\n***********************************" << endl;
        ltfast_on = true;
    }
    if (typeOfSpatialCorr == 10 && typeOfTemporalCorr == 10) {
        cout << "\n************************************\n";
        cout << "*      Running double unsharp      *";
        cout << "\n************************************" << endl;
        unsharp_on = true;
    }
    // NOTE: if you're interested, you could re-implement true top hats by inlining them here as well
    
    // UI to avoid accidentally overwriting angle files
    ifstream ifile("angleservo_cM_cIT.txt");
    int overwriteFile;
    if (ifile) {
        cout << "WARNING: An angle file for this program already exists. 1: continue and overwrite file, 0: kill program" << endl;
        cin >> overwriteFile;
        if (!overwriteFile)
            exit(0);
    }
    
    anglefile.open("angleservo_cM_cIT.txt", ios::out | ios::trunc); // file to plot angles in function of time
    
    // create (bake?) Loaf object using constructor
    loaf freshLoaf = loaf(numberOfSlices);
    
    float oldslice[13][11] = {{0}}; // stores the last configuration of paddles that was sent to the grid
    float step_size[13][11] = {{0}}; // stores the step size needed to get to the next configuration
    float newslice[13][11] = {{0}}; // stores the step size needed to get to the next configuration
    
    float rms;
    float correction = 1;
    int i = 1; // grid number counter
    int SPACING = 1; // number of interpolations needed to keep servo speed under its max value, in worst case
    
    float amplitude; // for steps/speeds calculations and safety checks, below
    float diff; // difference between two doubles (used with epsilon in place of == operator, which doesn't perform well on doubles)
    float EPSILON = 0.1; // error tolerance for double comparisons (just be within a tenth of a degree)
    int halfLoaf = numberOfSlices / (int) 2; // determine where to set middle slice of loaf (assumes numberOfSlices is odd)
    
    // Normalization calculations (complicated because different correlation functions have different methods for normalization)
    norm = 0; // master norm parameter that's accessible/used from runcorr_3D
    float bound;
    if (typeOfSpatialCorr <= 4 || typeOfSpatialCorr == 10 || typeOfSpatialCorr == 0) bound = range_of_corr;
    else bound = spatial_sigma;
    /* Declare function pointers for the spatial correlation functions
     First, we declare a pointer *pfSpatial/TemporalCorr to a function with 2 arguments j and k (or 1 argument t), plus parameters necessary for function operation (currently sigma and height, though more may need to be added if different kernels or a customizable unsharp are to be implemented).
     Then we call the function pickSpatialCorr (from pickCorrelations.cpp). This function returns a pointer to the function of choice (determined by 'mode'),
     which can then be used by simply calling pfSpatial/TemporalCorr(j, k). */
    
    /* NOTE: We attempted to import sigma, alpha, and height using the pickSpatialCorr/pickTemporalCorr functions, saving them as private global variables in pickCorrelations.cpp.
     This allowed us to simply call pfSpatialCorr(j, k) and pfTemporalCorr(t), removing excess arguments from the runcorr and computermscorr functions, and changing the implementation so that (except in the periodic case) the function pointers would only be created once and then would be passed to the necessary methods as arguments.
     We found this method to be much cleaner, but unfortunately we introduced a bug somewhere along the line. Rather than spend days trying to track it down, we reverted to our previous version and continued taking data. It might be worth re-implementing this feature, if you have time and want to make the code more streamlined. For reference, our previous attempt should still be in a series of commits on GitHub. ~ Nathan Wei and Kevin Griffin
     */
    
    float (*pfSpatialCorr)(int j, int k, float spatial_sigma, float spatial_alpha, float spatial_height);
    float (*pfTemporalCorr)(int t, float temporal_sigma, float temporal_alpha, float temporal_height);
    pfSpatialCorr = pickSpatialCorr(typeOfSpatialCorr);
    pfTemporalCorr = pickTemporalCorr(typeOfTemporalCorr);
    
    // Correlation function work for finding normalization
    // Note: this is different from the previous implementation, which had mysteriously different logic for each function
    for (int j = -bound; j <= bound; j++) { // range of neighbors used to compute normalization/convolution
        for (int k = -bound; k <= bound; k++) { // j and k refer to the shift
            for (int t = -halfLoaf; t <= halfLoaf; t++) {
                norm += (pfSpatialCorr(j, k, spatial_sigma, spatial_alpha, spatial_height) * pfTemporalCorr(t, temporal_sigma, temporal_alpha, temporal_height));
            }
        }
    }
    
    // makes a random correlated sequence of angles, with the same parameters but without correction
    // computes its mean and rms value of angles. This is done so that the rms correction factor can be
    // determined before the angles have been produced
    rms = compute_rmscorr_3D(spatial_sigma, temporal_sigma, typeOfSpatialCorr, typeOfTemporalCorr, spatial_alpha, temporal_alpha, spatial_height, temporal_height, 1, 1, halfLoaf);
    correction = target_rms / rms; // correction factor
    cout << "Done! Correction factor is " << correction << endl << "Setting up timing..." << endl;
    cout << "Done! Starting grid motions" << endl;
    
    //timing:
    // timing uses the standard timeval structure. a timeval struct holds seconds and remaining microseconds. This time is the number of seconds and remaining microseconds since 01.01.1970. Note: once microseconds reaches 10000000, seconds increments and microseconds is set to zero
    timeval startTime; // declare a structure for holding the time that the last slice of angles was sent to the grid
    timeval currentTime; // declare a structure for holding the current time
    long usecElapsed; // a varaible for holding the difference between currentTime and startTime
    gettimeofday(&startTime,0); // initialize startTime with the current time
    
    // ------------ (unclear whether this is actually necessary)
    gettimeofday(&currentTime,0);
    usecElapsed = (currentTime.tv_sec - startTime.tv_sec)*1000000 + ((signed long)currentTime.tv_usec - (signed long)startTime.tv_usec);
    while (usecElapsed <= updatetimeinmus){
        gettimeofday(&currentTime,0);
        usecElapsed = (currentTime.tv_sec - startTime.tv_sec)*1000000 + ((signed long)currentTime.tv_usec - (signed long)startTime.tv_usec);
    }
    while (usecElapsed > updatetimeinmus){
        gettimeofday(&currentTime,0);
        usecElapsed = (signed long)currentTime.tv_usec - (signed long)startTime.tv_usec;
    }
    // ------------
    
    bool firstTime = true; // prevent timing error message from showing up during first iteration
    // main loop: give angle orders (runs until ctrl-C is issued by user, which is caught by signal handler in menuII)
    while(0==0){
        
        //freshLoaf.Loaf_printFullArray(); // debugging
        cout << "\nGrid #" << i << " "; // print grid number
        i += 1;
        // get new slice of angles (using either inlined algorithm or standard function pointer algorithm)
        if (ltfast_on)
            ltfast(newslice, &freshLoaf, halfLoaf, spatial_sigma, temporal_sigma, spatial_alpha, temporal_alpha,
                   spatial_height, temporal_height, correction);
        else if (unsharp_on)
            unsharp(newslice, &freshLoaf, halfLoaf, spatial_sigma, temporal_sigma, spatial_alpha, temporal_alpha,
                    spatial_height, temporal_height, correction);
        else
            runcorr_3D(newslice, &freshLoaf, halfLoaf, spatial_sigma, temporal_sigma, spatial_alpha, temporal_alpha,
                       spatial_height, temporal_height, typeOfSpatialCorr, typeOfTemporalCorr, 0, 0, correction);
        
        // store necessary servo speeds after carrying out safety checks
        for (int row = 0; row < 11; row++) {
            for (int col = 0; col < 13; col++){
                numberOfAnglesSet++; // total number of paddles moved, since the beginning of time (global variable)
                // angle safety processing: do not exceed angle of 90 degrees
                if (newslice[col][row]>90){
                    newslice[col][row]=90;
                    over90orminus90count++; // keep track of number of angles that need to be kept from exceeding +/- 90 degrees
                    //cout << "+"; // print visual representation within histogram if you want
                }
                else if (newslice[col][row]<-90){
                    newslice[col][row]=-90;
                    over90orminus90count++;
                    //cout << "-";
                }
                
                amplitude = newslice[col][row] - oldslice[col][row]; // calculate the amplitude between the old and the new angles
                if (fabs(amplitude)/(max_speed) > SPACING) { // constrain paddles that have to move too far for the servos to handle
                    cout << "*"; // constraining histogram - 1 * = 1 paddle constrained by max speed
                    outOfBoundsCount++; // keep track of number of constrained paddles
                    if (amplitude > 0) step_size[col][row] = max_speed;
                    else if (amplitude < 0) step_size[col][row] = -max_speed;
                }
                // move at constant speed (equal-sized steps over SPACING interpolations) to the target angle
                // this implemenation assumes that servos don't have a minimum speed or angle
                else step_size[col][row] = amplitude/(SPACING);
            }
        }
        
        /* Create SPACING timeslices to separate old and new configurations, and feed each one to the grid in succession.
         * This decouples the grid updating timescale from the correlation computation timescale, and also means we call the
         *   computationally-expensive runcorr_3D method once every SPACING grid configurations.
         * SPACING is currently set to 1, in order to get the highest frequency of correlations possible. With faster computers or a more efficient
         *   runcorr_3D algorithm, increasing SPACING to improve the resolution of the grid motions (while still running 10 correlations per second)
         *   might be possible. */
        for (int t = 0; t < SPACING; t++) {
            //cout << " " << (t+1); // print interpolation number
            
            // compute new intermediate grid position, with steps necessary to attain it
            // adding the calculated servo speed every 0.1 seconds SPACING times guarantees the angle will reach its target by the end of the iteration
            for (int row = 0; row < 11; row++) {
                for(int col = 0; col < 13; col++){
                    if (step_size[col][row] != 0) { // don't bother checking servos that have already arrived
                        diff = fabs(newslice[col][row] - oldslice[col][row]); // determination of approximate equality for doubles
                        if (diff > EPSILON) oldslice[col][row] += step_size[col][row]; // not equal -> add another step
                        else step_size[col][row] = 0; // paddle has arrived; tell servo not to move any more
                    }
                }
            }
            
            //set position of each servo (within 0.1 seconds)
            gettimeofday(&currentTime,0); // set currentTime to hold the current time
            usecElapsed = (currentTime.tv_sec - startTime.tv_sec)*1000000 + ((signed long)currentTime.tv_usec - (signed long)startTime.tv_usec);// useconds elapsed since startTime
            
            if (usecElapsed > updatetimeinmus && !firstTime){ // no need to wait because runcorr took more than .1 sec
                cout << "Time Elapsed is greater than .1 sec.  Time Elapsed = " << usecElapsed;
            }
            else if (usecElapsed < 0){
                assert(0); // assert because something bizarre happened (maybe the timer overflowed somehow) - safety for a case that shouldn't ever occur
            }
            else {
                cout << " " << usecElapsed << "\u03BCs"; // print computation time in microseconds (at end of histogram)
                while (usecElapsed < updatetimeinmus){ // we need to wait
                    gettimeofday(&currentTime,0);
                    usecElapsed = (currentTime.tv_sec - startTime.tv_sec)*1000000 + ((signed long)currentTime.tv_usec - (signed long)startTime.tv_usec);
                }
            }
            // record new start time to compare against for next iteration
            gettimeofday(&startTime,0);
            
            if (firstTime) firstTime = false; // set to false for remainder of run
            
            //----------------
            setanglestoallservosIII(oldslice, step_size, constantArea, target_rms); // for motion - send generated 2D arrays to grid communication method
        }
    }
    anglefile.close(); // never reaches this point
    return 0; // never reaches this point
}
