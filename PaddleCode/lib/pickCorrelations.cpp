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

// Gaussian convolution function
float gaussianSpatialCorr(int j, int k, float *ptr_to_norm, float *ptr_to_norm1, float spatial_sigma, float height){
    //cout << "1s" << endl; // debugging
    return (float)(exp(-((j*j) + (k*k)) / (2 * spatial_sigma*spatial_sigma)) / *ptr_to_norm);
}

float gaussianTemporalCorr(int t, float *ptr_to_norm, float *ptr_to_norm1, float temporal_sigma, float height){
    //cout << "1t" << endl; // debugging
    return (float)(exp(-(t*t)/(2*temporal_sigma*temporal_sigma)) / *ptr_to_norm);
}

// 1/r^2 convolution function
float inverseSquareSpatialCorr(int j, int k, float *ptr_to_norm, float *ptr_to_norm1, float spatial_sigma, float height){
    //cout << "2s" << endl; // debugging
    return pow( (sqrt(j*j + k*k) + 1) , -2);
}

float inverseSquareTemporalCorr(int t, float *ptr_to_norm, float *ptr_to_norm1, float temporal_sigma, float height){
    //cout << "2t" << endl; // debugging
    return pow( (abs(t) + 1) , -2);
}

// 1/r^3 convolution function
float inverseCubeSpatialCorr(int j, int k, float *ptr_to_norm, float *ptr_to_norm1, float spatial_sigma, float height){
    //cout << "3s" << endl; // debugging
    return pow( (sqrt(j*j + k*k) + 1) , -3);
}

float inverseCubeTemporalCorr(int t, float *ptr_to_norm, float *ptr_to_norm1, float temporal_sigma, float height){
    //cout << "3t" << endl; // debugging
    return pow( (abs(t) + 1) , -3);
}

// 1/r^4 convolution function
float inverseQuarticSpatialCorr(int j, int k, float *ptr_to_norm, float *ptr_to_norm1, float spatial_sigma, float height){
    //cout << "4s" << endl; // debugging
    return pow( (sqrt(j*j + k*k) + 1) , -4);
}

float inverseQuarticTemporalCorr(int t, float *ptr_to_norm, float *ptr_to_norm1, float temporal_sigma, float height){
    //cout << "4t" << endl; // debugging
    return pow( (abs(t) + 1) , -4);
}

// top hat convolution function
float topHatSpatialCorr(int j, int k, float *ptr_to_norm, float *ptr_to_norm1, float spatial_sigma, float height){
    //cout << "5s" << endl; // debugging
    double dist = sqrt((j*j) + (k*k));
    if (dist <= spatial_sigma){
        (*ptr_to_norm1)++;
        return 1.0;
    }
    else 
        return 0;
}

float topHatTemporalCorr(int t, float *ptr_to_norm, float *ptr_to_norm1, float temporal_sigma, float height){
    //cout << "5t" << endl; // debugging
    if (abs(t) <= temporal_sigma){
        (*ptr_to_norm1)++;
        return 1.0;
    }
    else 
        return 0;
}

// true top hat with one main paddle, no wrapping around
float trueTopHatSpatialCorr(int j, int k, float *ptr_to_norm, float *ptr_to_norm1, float spatial_sigma, float height){
    //cout << "6s" << endl; // debugging
    return 0;
}

float trueTopHatTemporalCorr(int t, float *ptr_to_norm, float *ptr_to_norm1, float temporal_sigma, float height){
    //cout << "6t" << endl; // debugging
    return 0;
}

// true top hat with one randomly chosen paddle
float trueTopHatRandomSpatialCorr(int j, int k, float *ptr_to_norm, float *ptr_to_norm1, float spatial_sigma, float height){
    //cout << "7s" << endl; // debugging
    return 0;
}

float trueTopHatRandomTemporalCorr(int t, float *ptr_to_norm, float *ptr_to_norm1, float temporal_sigma, float height){
    //cout << "7t" << endl; // debugging
    return 0;
}

float topHatLongTailSpatialCorr(int j, int k, float *ptr_to_norm, float *ptr_to_norm1, float alpha, float height){
    //cout << "8s" << endl; // debugging
    double dist = sqrt((j*j)+(k*k));
    if (dist <= alpha) {
        (*ptr_to_norm)++;
        return 1.0;
    }
    else {
        (*ptr_to_norm1)++;
        return height;
    }
}

float topHatLongTailTemporalCorr(int t, float *ptr_to_norm, float *ptr_to_norm1, float alpha, float height){
    //cout << "8t" << endl; // debugging
    if (abs(t) <= alpha) {
        (*ptr_to_norm)++;
        return 1.0;
    }
    else {
        (*ptr_to_norm1)++;
        return height;
    }
}

float triangleSpatialCorr(int j, int k, float *ptr_to_norm, float *ptr_to_norm1, float spatial_sigma, float height){
    //cout << "9s" << endl; // debugging
    double dist = sqrt((j*j)+(k*k));
    if (dist <= spatial_sigma) {
        (*ptr_to_norm1)++;
        return (-1) / spatial_sigma * dist + 1;
    }
    return 0;
}

float triangleTemporalCorr(int t, float *ptr_to_norm, float *ptr_to_norm1, float temporal_sigma, float height){
    //cout << "9t" << endl; // debugging
    if (abs(t) <= temporal_sigma) {
        (*ptr_to_norm1)++;
        return (-1) / temporal_sigma * abs(t) + 1;
  }
    return 0;
}

/*------------------------------------------------------------------------*/

float (*pickSpatialCorr(int typeOfSpatialCorr)) (int j, int k, float *ptr_to_norm, float *ptr_to_norm1, float spatial_sigma, float height){
  // validate parameters
  assert (typeOfSpatialCorr > 0 && typeOfSpatialCorr <= 9);
  
  if(typeOfSpatialCorr == 1)
    return gaussianSpatialCorr;
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
  else
    return &triangleSpatialCorr;
  
  // should never reach this point
  assert(1);
  return NULL;
}

/*------------------------------------------------------------------------*/

/* The function pickTemporalCorr accepts an integer named typeOfTemporalCorr
 * and returns a pointer to a function (which accepts an int, int, float*,float as parameters
 * and returns a float)
 * In other words, this function returns a correlation function.
 */
float (*pickTemporalCorr(int typeOfTemporalCorr)) (int t, float *ptr_to_norm, float *ptr_to_norm1, float temporal_sigma, float height){
  // validate parameters
  assert (typeOfTemporalCorr > 0 && typeOfTemporalCorr <= 9);
  
  if(typeOfTemporalCorr == 1)
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
  else
    return &triangleTemporalCorr;
  
  // should never reach this point
  assert(1);
  return NULL;
}
