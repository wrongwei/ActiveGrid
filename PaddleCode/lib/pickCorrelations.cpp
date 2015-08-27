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
/* For details of how to use this module, please see testPickCorrelations.cpp.
 * The basic idea is that you have an int that specifies which type of correlation
 * you want to perform in the spatial dimension. You call pickSpatialCorr, which
 * accepts that int as an arg and retruns you a ptr to the appropriate correlation
 * for that integer. For example, let's say you pass pickSpatialCorr the integer 1.
 * pickSpatialCorr will return a ptr to the function gaussianSpatialCorrelation.
 * pickTemporalCorr works in the same way but returns a ptr to the appropriate
 * temporal correlation function. So, given the argument 5, it returns a ptr to
 * topHatTemporalCorr.
 * NOTE: All correlation functions require the parameters sigma, alpha, and height
 * in order for the function pointer implementation to work (even if they're not used)
 */
/*------------------------------------------------------------------------*/


// Random (i.e. no correlation)
// We have not implemented randomTemporalCorr in menuII because this would lead to many angles
// being asked to move more than 40 degrees per 0.1 seconds (i.e. lots of constraining and problems)
float randomSpatialCorr(int j, int k, float spatial_sigma, float spatial_alpha, float height){
    if (j == 0 && k == 0)
        return 1;
    return 0;
}
// WEIRD AND DANGEROUS: USE AT YOUR OWN RISK
float randomTemporalCorr(int t, float temporal_sigma, float temporal_alpha, float height){
    if (t == 0)
        return 1;
    return 0;
}

// Gaussian convolution function
float gaussianSpatialCorr(int j, int k, float spatial_sigma, float spatial_alpha, float height){
    return expf(-(((j*j) + (k*k)) / (2 * spatial_sigma*spatial_sigma)));
}

float gaussianTemporalCorr(int t, float temporal_sigma, float temporal_alpha, float height){
    return expf(-((t*t) / (2 * temporal_sigma*temporal_sigma)));
}

// 1/r^2 convolution function
float inverseSquareSpatialCorr(int j, int k, float spatial_sigma, float spatial_alpha, float height){
    return pow( (sqrt(j*j + k*k) + 1) , -2);
}

float inverseSquareTemporalCorr(int t, float temporal_sigma, float temporal_alpha, float height){
    return pow( (abs(t) + 1) , -2);
}

// 1/r^3 convolution function
float inverseCubeSpatialCorr(int j, int k, float spatial_sigma, float spatial_alpha, float height){
    return pow( (sqrt(j*j + k*k) + 1) , -3);
}

float inverseCubeTemporalCorr(int t, float temporal_sigma, float temporal_alpha, float height){
    return pow( (abs(t) + 1) , -3);
}

// 1/r^4 convolution function
float inverseQuarticSpatialCorr(int j, int k, float spatial_sigma, float spatial_alpha, float height){
    return pow( (sqrt(j*j + k*k) + 1) , -4);
}

float inverseQuarticTemporalCorr(int t, float temporal_sigma, float temporal_alpha, float height){
    return pow( (abs(t) + 1) , -4);
}

// top hat convolution function
float topHatSpatialCorr(int j, int k, float spatial_sigma, float spatial_alpha, float height){
    double dist = sqrt((j*j) + (k*k));
    if (dist <= spatial_sigma) return 1.0;
    else return 0;
}

float topHatTemporalCorr(int t, float temporal_sigma, float temporal_alpha, float height){
    if (abs(t) <= temporal_sigma) return 1.0;
    else return 0;
}

// true top hat with one main paddle, no wrapping around (NOT IMPLEMENTED IN 3D)
float trueTopHatSpatialCorr(int j, int k, float spatial_sigma, float spatial_alpha, float height){
    return 0;
}

float trueTopHatTemporalCorr(int t, float temporal_sigma, float temporal_alpha, float height){
    return 0;
}

// true top hat with one randomly chosen paddle (NOT IMPLEMENTED IN 3D)
float trueTopHatRandomSpatialCorr(int j, int k, float spatial_sigma, float spatial_alpha, float height){
    return 0;
}

float trueTopHatRandomTemporalCorr(int t, float temporal_sigma, float temporal_alpha, float height){
    return 0;
}

// top hat with center height of 1 (width alpha) and tails between alpha and sigma of height 'height'
float topHatLongTailSpatialCorr(int j, int k, float spatial_sigma, float spatial_alpha, float height){
    float dist = sqrt((j*j)+(k*k));
    if (dist <= spatial_alpha) return 1.0;
    if (dist <= spatial_sigma) return height;
    return 0;
}

float topHatLongTailTemporalCorr(int t, float temporal_sigma, float temporal_alpha, float height){
    float abs_t = abs(t);
    if (abs_t <= temporal_alpha) return 1.0;
    if (abs_t <= temporal_sigma) return height;
    return 0;
}

// triangular correlation function, with max height of 1 and base of width sigma
float triangleSpatialCorr(int j, int k, float spatial_sigma, float spatial_alpha, float height){
    double dist = sqrt((j*j)+(k*k));
    if (dist <= spatial_sigma)
        return ((-1) / spatial_sigma * dist + 1);
    return 0;
}

float triangleTemporalCorr(int t, float temporal_sigma, float temporal_alpha, float height){
    if (abs(t) <= temporal_sigma)
        return ((-1) / temporal_sigma * abs(t) + 1);
    return 0;
}

// unsharp correlations. Unsharp correlations are equal to the image in the
// center and have a negative long tail outside alpha (but within sigma).
// The effect for this correlation is that the grid is made less correlated,
// because panels will be made more different to the neighbors; if you think of it
// as a form of image processing it makes the image sharper, which is the opposite
// of the gaussian (which blurs the image)
// Depth is the scaling sharpness (0->1). Larger depth means a sharper image (or sharper contrast)
float unsharpSpatialCorr(int j, int k, float spatial_sigma, float spatial_alpha, float depth){
    float dist = sqrt(j*j + k*k);
    if (dist <= spatial_alpha) return 1;
    if (dist <= spatial_sigma) return -depth;
    return 0;
    //return height * -expf(-(((j*j) + (k*k)) / (2 * spatial_sigma*spatial_sigma))); // previous Gaussian implementation
}

float unsharpTemporalCorr(int t, float temporal_sigma, float temporal_alpha, float depth){
    if (t <= temporal_alpha) return 1;
    if (t <= temporal_sigma) return -depth;
    return 0;
    //return height * -expf(-((t*t) / (2 * temporal_sigma*temporal_sigma))); // previous Gaussian implementation
}

/*------------------------------------------------------------------------*/

/* The function pickSpatialCorr accepts an integer named typeOfSpatialCorr
 * and returns a pointer to a function (which accepts 2 ints and 3 floats as parameters
 * and returns a float)
 * In other words, this function returns a correlation function, which can then be used
 *  as an ordinary function (by calling it with the pointer name and feeding it arguments)
 */
float (*pickSpatialCorr(int typeOfSpatialCorr)) (int j, int k, float spatial_sigma, float spatial_alpha, float height){
    // validate parameters
    assert (typeOfSpatialCorr >= 0 && typeOfSpatialCorr <= 10);
    
    if(typeOfSpatialCorr == 0)
        return randomSpatialCorr;
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
    else if(typeOfSpatialCorr == 10)
        return &unsharpSpatialCorr;
    
    // should never reach this point
    assert(1);
    return NULL;
}

/*------------------------------------------------------------------------*/

/* The function pickTemporalCorr accepts an integer named typeOfTemporalCorr
 * and returns a pointer to a function (which accepts an int and 3 floats as parameters
 * and returns a float)
 * In other words, this function returns a correlation function (see above).
 */
float (*pickTemporalCorr(int typeOfTemporalCorr)) (int t, float temporal_sigma, float temporal_alpha, float height){
    // validate parameters
    assert (typeOfTemporalCorr >= 0 && typeOfTemporalCorr <= 10);
    
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
    else if(typeOfTemporalCorr == 10)
        return &unsharpTemporalCorr;
    
    // should never reach this point
    assert(1);
    return NULL;
}
