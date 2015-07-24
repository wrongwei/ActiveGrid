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
 * guassian TemporalCorr. 
 */
/*------------------------------------------------------------------------*/

// gaussian correlation function in the spatial dimensions
float gaussianSpatialCorr(int j, int k, float *ptr_to_norm, float spatial_sigma){
  return (float)(exp(-((j*j) + (k*k)) / (2 * spatial_sigma*spatial_sigma)));
}

// gaussian correlation function in the temporal dimension
float gaussianTemporalCorr(int t, float *ptr_to_norm, float temporal_sigma){
  return (float)(exp(-(t*t)/(2*temporal_sigma*temporal_sigma)));
}

// correlation function in the spatial dimensions
float inverseSquareSpatialCorr(int j, int k, float *ptr_to_norm, float spatial_sigma){
  return pow( (sqrt(j*j + k*k) + 1) , -2);
}

// correlation function in the temporal dimension
float inverseSquareTemporalCorr(int t, float *ptr_to_norm, float temporal_sigma){
  return pow( (abs(t) + 1) , -2);
}

// correlation function in the spatial dimensions
float inverseCubeSpatialCorr(int j, int k, float *ptr_to_norm, float spatial_sigma){
  return pow( (sqrt(j*j + k*k) + 1) , -3);
}

// correlation function in the temporal dimension
float inverseCubeTemporalCorr(int t, float *ptr_to_norm, float temporal_sigma){
  return pow( (abs(t) + 1) , -3);
}

// correlation function in the spatial dimensions
float inverseQuarticSpatialCorr(int j, int k, float *ptr_to_norm, float spatial_sigma){
  return pow( (sqrt(j*j + k*k) + 1) , -4);
}

// correlation function in the temporal dimension
float inverseQuarticTemporalCorr(int t, float *ptr_to_norm, float temporal_sigma){
  return pow( (abs(t) + 1) , -4);
}

// correlation function in the spatial dimensions
float topHatSpatialCorr(int j, int k, float *ptr_to_norm, float spatial_sigma){
  return 0;
}

// correlation function in the temporal dimension
float topHatTemporalCorr(int t, float *ptr_to_norm, float temporal_sigma){
  return 0;
}

// correlation function in the spatial dimensions
float trueTopHatSpatialCorr(int j, int k, float *ptr_to_norm, float spatial_sigma){
  return 0;
}

// correlation function in the temporal dimension
float trueTopHatTemporalCorr(int t, float *ptr_to_norm, float temporal_sigma){
  return 0;
}

// correlation function in the spatial dimensions
float trueTopHatRandomSpatialCorr(int j, int k, float *ptr_to_norm, float spatial_sigma){
  return 0;
}

// correlation function in the temporal dimension
float trueTopHatRandomTemporalCorr(int t, float *ptr_to_norm, float temporal_sigma){
  return 0;
}

// correlation function in the spatial dimensions
float topHatLongTailSpatialCorr(int j, int k, float *ptr_to_norm, float spatial_sigma){
  return 0;
}

// correlation function in the temporal dimension
float topHatLongTailTemporalCorr(int t, float *ptr_to_norm, float temporal_sigma){
  return 0;
}

// correlation function in the spatial dimensions
float triangleSpatialCorr(int j, int k, float *ptr_to_norm, float spatial_sigma){
  return 0;
}

// correlation function in the temporal dimension
float triangleTemporalCorr(int t, float *ptr_to_norm, float temporal_sigma){
  return 0;
}

/*------------------------------------------------------------------------*/
/* The function pickSpatialCorr accepts an integer named typeOfSpatialCorr
 * and returns a pointer to a function (which accepts an int, int, float*,float as parameters
 * and returns a float)
 * In other words, this function returns a correlation function.
 */
float (*pickSpatialCorr(int typeOfSpatialCorr)) (int j, int k, float *ptr_to_norm, float spatial_sigma){
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
float (*pickTemporalCorr(int typeOfTemporalCorr)) (int t, float *ptr_to_norm, float temporal_sigma){
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
