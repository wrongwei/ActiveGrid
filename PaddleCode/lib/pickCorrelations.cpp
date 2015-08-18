/*------------------------------------------------------------------------*/
/* Author: Kevin Griffin                                                  */
/* kevinpg@princeton.edu                                                  */
/* 22/07/15                                                               */
/* pickCorrelations.cpp                                                   */
/*------------------------------------------------------------------------*/
#include "pickCorrelations.h"
#include <assert.h>
#include <cmath>
/*------------------------------------------------------------------------*/
/* For details of how to use this module, please see testPickCorrelations.cpp
 * The basic idea is that you have an int that specifies which type of correlation
 * you want to perform in the spatial dimension. You call pickSpatialCorr and
 * accepts that int as an arg and retruns you a ptr to the appropriate correlation
 * for that integer. For example, let's say you pass pickSpatialCorr the integer 1.
 * pickSpatialCorr will return a ptr to the function gaussianSpatialCorrelation.
 * pcikTemporalCorr works in the same way but returns a ptr to the appropriate
 * temporal correlation function. So, which the argument 1, it returns a ptr to
 * gaussian TemporalCorr.
 */
/*------------------------------------------------------------------------*/
/* set up some global variables that are initialized in pickTemporalCorr
 and pickSpatialCorr and then used in the correlation functions */
static float spatialSigmaLimit = 0;
static float temporalSigmaLimit = 0;
static float spatialAlphaLimit = 0;
static float temporalAlphaLimit = 0;
static float spatialHeight = 0;
static float temporalHeight = 0;
static float twoTimesSpatialSigmaSquared = 0;
static float twoTimesTemporalSigmaSquared = 0;
static float unsharpCenterPaddleHeightSpatial = 0;
static float unsharpCenterPaddleHeightTemporal = 0;

/*------------------------------------------------------------------------*/

// random. In other words, no correlation
// We have not implemented randomTemporalCorr in menuII because this would lead to many angles
// being asked to move more than 40 degrees per sec
float randomSpatialCorr(int j, int k){
    if (j == 0 && k == 0)
        return 1;
    return 0;
}

float randomTemporalCorr(int t){
    if (t == 0)
        return 1;
    return 0;
}

// Gaussian convolution function
float gaussianSpatialCorr(int j, int k){
    return expf(-(((j*j) + (k*k)) / (twoTimesSpatialSigmaSquared)));
}

float gaussianTemporalCorr(int t){
    return expf(-((t*t) / (twoTimesTemporalSigmaSquared)));
}

// 1/r^2 convolution function
float inverseSquareSpatialCorr(int j, int k){
    return pow( (sqrt(j*j + k*k) + 1) , -2);
}

float inverseSquareTemporalCorr(int t){
    return pow( (abs(t) + 1) , -2);
}

// 1/r^3 convolution function
float inverseCubeSpatialCorr(int j, int k){
    return pow( (sqrt(j*j + k*k) + 1) , -3);
}

float inverseCubeTemporalCorr(int t){
    return pow( (abs(t) + 1) , -3);
}

// 1/r^4 convolution function
float inverseQuarticSpatialCorr(int j, int k){
    return pow( (sqrt(j*j + k*k) + 1) , -4);
}

float inverseQuarticTemporalCorr(int t){
    return pow( (abs(t) + 1) , -4);
}

// top hat convolution function
float topHatSpatialCorr(int j, int k){
    double dist = sqrt((j*j) + (k*k));
    if (dist <= spatialSigmaLimit) return 1.0;
    else return 0;
}

float topHatTemporalCorr(int t){
    if (abs(t) <= temporalSigmaLimit) return 1.0;
    else return 0;
}

// true top hat with one main paddle, no wrapping around
float trueTopHatSpatialCorr(int j, int k){
    //this function is not ready yet, so it just throws an error
    assert(0);
    return 0;
}

float trueTopHatTemporalCorr(int t){
    //this function is not ready yet, so it just throws an error
    assert(0);
    return 0;
}

// true top hat with one randomly chosen paddle
float trueTopHatRandomSpatialCorr(int j, int k){
    //this function is not ready yet, so it just throws an error
    assert(0);
    return 0;
}

float trueTopHatRandomTemporalCorr(int t){
    //this function is not ready yet, so it just throws an error
    assert(0);
    return 0;
}

float topHatLongTailSpatialCorr(int j, int k){
    double dist = sqrt((j*j)+(k*k));
    if (dist <= spatialAlphaLimit)
        return 1.0;
    if (dist <= spatialSigmaLimit)
        return spatialHeight;
    return 0;
}

float topHatLongTailTemporalCorr(int t){
    float absoluteVal_t = abs(t);
    if (absoluteVal_t <= temporalAlphaLimit)
        return 1.0;
    if (absoluteVal_t <= temporalSigmaLimit)
        return temporalHeight;
    return 0;
}

float triangleSpatialCorr(int j, int k){
    double dist = sqrt((j*j)+(k*k));
    if (dist <= spatialSigmaLimit)
        return ((-1) / spatialSigmaLimit * dist + 1);
    return 0;
}

float triangleTemporalCorr(int t){
    float absoluteVal_t = abs(t);
    if (absoluteVal_t <= temporalSigmaLimit)
        return (-1 / temporalSigmaLimit * absoluteVal_t + 1);
    return 0;
}

// unsharp correlations. Unsharp correlations are equal to the image in the
// center and are equal to a negative gaussian everywhere else.
// The effect for this correlation is that the grid is made less correlated,
// because panels will be made more different to the neighbors; if you think of it
// as a form of image processing it makes the image sharper, which is the opposite
// of the gaussian (which blurs the image)
// Height is the scaling sharpness. Higher height means a sharper image (or sharper contrast)
float unsharpSpatialCorr(int j, int k){
    float dist = sqrt(j*j + k*k);
    if (dist <= spatialAlphaLimit)
        return unsharpCenterPaddleHeightSpatial;
    //the negative sign for the tail is added by special logic in runcorr
    if (dist <= spatialSigmaLimit)
        return spatialHeight;
    return 0;
}

float unsharpTemporalCorr(int t){
    if (t <= temporalAlphaLimit){
        return unsharpCenterPaddleHeightTemporal;
    }
    //the negative sign for the tail is added by special logic in runcorr
    if (t <= temporalSigmaLimit)
        return temporalHeight;
    return 0;
}

/*------------------------------------------------------------------------*/

float (*pickSpatialCorr(int typeOfSpatialCorr, float spatial_sigma, float spatial_alpha, float spatial_height)) (int j, int k) {
    // validate parameters
    assert (typeOfSpatialCorr >= 0 && typeOfSpatialCorr <= 10);
    
    spatialSigmaLimit = spatial_sigma;
    spatialAlphaLimit = spatial_alpha;
    spatialHeight = spatial_height;
    twoTimesSpatialSigmaSquared = 2 * spatial_sigma * spatial_sigma;
    unsharpCenterPaddleHeightSpatial = 1;
    
    if(typeOfSpatialCorr == 0)
        return &randomSpatialCorr;
    else if(typeOfSpatialCorr == 1)
        return &gaussianSpatialCorr;
    else if(typeOfSpatialCorr == 2)
        return &inverseSquareSpatialCorr;
    else if(typeOfSpatialCorr == 3)
        return &inverseCubeSpatialCorr;
    else if(typeOfSpatialCorr == 4)
        return &inverseQuarticSpatialCorr;
    else if(typeOfSpatialCorr == 5)
        return &topHatSpatialCorr;
    else if(typeOfSpatialCorr == 6)
        return &trueTopHatSpatialCorr;
    else if(typeOfSpatialCorr == 7)
        return &trueTopHatSpatialCorr;
    else if(typeOfSpatialCorr == 8)
        return &topHatLongTailSpatialCorr;
    else if(typeOfSpatialCorr == 9)
        return &triangleSpatialCorr;
    else if(typeOfSpatialCorr == 10) {
        // calculate areas of positive and negative sections of unsharp kernel
        float inside_area = 0;
        float outside_area = 0;
        float dist = 0;
        for (int row = -spatial_sigma; row <= spatial_sigma; row++) {
            for (int col = -spatial_sigma; col <= spatial_sigma; col++) {
                dist = sqrt((col*col) + (row*row));
                if (dist <= spatial_alpha) inside_area++;
                else if (dist <= spatial_sigma) outside_area++;
            }
        }
        // calculate scaling of positive region (so that integral over kernel is zero)
        //unsharpCenterPaddleHeightSpatial = outside_area * spatial_height / inside_area;
        return &unsharpSpatialCorr;
    }
            
    // should never reach this point
    assert(1);
    return NULL;
}
                    
/*------------------------------------------------------------------------*/

/* The function pickTemporalCorr accepts an integer named typeOfTemporalCorr, with associated parameters needed
 * for computations, and returns a pointer to a function (which accepts 2 ints as parameters and returns a float)
 * In other words, this function returns a correlation function.
 */

float (*pickTemporalCorr(int typeOfTemporalCorr, float temporal_sigma, float temporal_alpha, float temporal_height)) (int t){
    // validate parameters
    assert (typeOfTemporalCorr >= 0 && typeOfTemporalCorr <= 10);
    
    temporalSigmaLimit = temporal_sigma;
    temporalAlphaLimit = temporal_alpha;
    temporalHeight = temporal_height;
    twoTimesTemporalSigmaSquared = 2 * temporal_sigma * temporal_sigma;
    unsharpCenterPaddleHeightTemporal = 1;
    
    if(typeOfTemporalCorr == 0)
        return &randomTemporalCorr;
    else if(typeOfTemporalCorr == 1)
        return &gaussianTemporalCorr;
    else if(typeOfTemporalCorr == 2)
        return &inverseSquareTemporalCorr;
    else if(typeOfTemporalCorr == 3)
        return &inverseCubeTemporalCorr;
    else if(typeOfTemporalCorr == 4)
        return &inverseQuarticTemporalCorr;
    else if(typeOfTemporalCorr == 5)
        return &topHatTemporalCorr;
    else if(typeOfTemporalCorr == 6)
        return &trueTopHatTemporalCorr;
    else if(typeOfTemporalCorr == 7)
        return &trueTopHatRandomTemporalCorr;
    else if(typeOfTemporalCorr == 8)
        return &topHatLongTailTemporalCorr;
    else if(typeOfTemporalCorr == 9)
        return &triangleTemporalCorr;
    else if(typeOfTemporalCorr == 10) {
        // calculate scaling of positive region (so that integral over kernel is zero)
        //unsharpCenterPaddleHeightTemporal = temporal_height * (temporal_sigma - temporal_alpha) / (temporal_alpha + 0.5);
        return &unsharpTemporalCorr;
    }
    
    // should never reach this point
    assert(1);
    return NULL;
}
