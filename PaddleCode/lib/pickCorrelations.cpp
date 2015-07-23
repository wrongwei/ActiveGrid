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

float gaussianSpatialCorr(int j, int k, float *ptr_to_norm, float *ptr_to_norm1, float spatial_sigma, float height){
  return (float)(exp(-((j*j) + (k*k)) / (2 * spatial_sigma*spatial_sigma)) / *ptr_to_norm);
}

float gaussianTemporalCorr(int t, float *ptr_to_norm, float *ptr_to_norm1, float temporal_sigma, float height){
  return (float)(exp(-(t*t)/(2*temporal_sigma*temporal_sigma)) / *ptr_to_norm);
}

float inverseSquareSpatialCorr(int j, int k, float *ptr_to_norm, float *ptr_to_norm1, float spatial_sigma, float height){
  return pow( (sqrt(j*j + k*k) + 1) , -2);
}

float inverseSquareTemporalCorr(int t, float *ptr_to_norm, float *ptr_to_norm1, float temporal_sigma, float height){
  return pow( (abs(t) + 1) , -2);
}

float inverseCubeSpatialCorr(int j, int k, float *ptr_to_norm, float *ptr_to_norm1, float spatial_sigma, float height){
  return pow( (sqrt(j*j + k*k) + 1) , -3);
}

float inverseCubeTemporalCorr(int t, float *ptr_to_norm, float *ptr_to_norm1, float temporal_sigma, float height){
  return pow( (abs(t) + 1) , -3);
}

float inverseQuarticSpatialCorr(int j, int k, float *ptr_to_norm, float *ptr_to_norm1, float spatial_sigma, float height){
  return pow( (sqrt(j*j + k*k) + 1) , -4);
}

float inverseQuarticTemporalCorr(int t, float *ptr_to_norm, float *ptr_to_norm1, float temporal_sigma, float height){
  return pow( (abs(t) + 1) , -4);
}

float topHatSpatialCorr(int j, int k, float *ptr_to_norm, float *ptr_to_norm1, float spatial_sigma, float height){
    double dist = sqrt((j*j) + (k*k));
    if (dist <= spatial_sigma){
        *ptr_to_norm1++;
        return 1.0;
    }
    else return 0;
}

float topHatTemporalCorr(int t, float *ptr_to_norm, float *ptr_to_norm1, float *ptr_to_norm1float temporal_sigma, float height){
    if (fabs(t) <= temporal_sigma){
        *ptr_to_norm1++;
        return 1.0;
    }
    else return 0;
}

// true top hat with one main paddle, no wrapping around
float trueTopHatSpatialCorr(int j, int k, float *ptr_to_norm, float *ptr_to_norm1, float spatial_sigma, float height){
  return 0;
}

float trueTopHatTemporalCorr(int t, float *ptr_to_norm, float *ptr_to_norm1, float temporal_sigma, float height){
  return 0;
}

// true top hat with one randomly chosen paddle
float trueTopHatRandomSpatialCorr(int j, int k, float *ptr_to_norm, float *ptr_to_norm1, float spatial_sigma, float height){
  return 0;
}

float trueTopHatRandomTemporalCorr(int t, float *ptr_to_norm, float *ptr_to_norm1, float temporal_sigma, float height){
  return 0;
}

float topHatLongTailSpatialCorr(int j, int k, float *ptr_to_norm, float *ptr_to_norm1, float alpha, float height){
    double dist = sqrt((j*j))+(k*k));
    if (dist <= alpha) {
        *ptr_to_norm++;
        return 1.0;
    } else {
        *ptr_to_norm1++;
        return height;
    }
}

float topHatLongTailTemporalCorr(int t, float *ptr_to_norm, float *ptr_to_norm1, float temporal_sigma, float height){
  return 0;
}

float triangleSpatialCorr(int j, int k, float *ptr_to_norm, float *ptr_to_norm1, float spatial_sigma, float height){
    double dist = sqrt((j*j))+(k*k));
    if (dist <= sigma) {
        *ptr_to_norm1++;
        return (-1) / spatial_sigma * dist + 1;
    }
}

float triangleTemporalCorr(int t, float *ptr_to_norm, float *ptr_to_norm1, float temporal_sigma, float height){
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
