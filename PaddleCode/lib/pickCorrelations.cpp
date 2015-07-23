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

float gaussianSpatialCorr(int j, int k, float *ptr_to_norm, float spatial_sigma){
  return (float)exp(-(pow((float)j,(int) 2)+ pow((float)k,(int)2))/(2* pow(spatial_sigma,2)));
}

float gaussianTemporalCorr(int t, float *ptr_to_norm, float temporal_sigma){
  return (float)exp(-(pow((float)t,(int) 2))/(2* pow(temporal_sigma,2)));
}

float inverseSquareSpatialCorr(int j, int k, float *ptr_to_norm, float spatial_sigma){
  return pow( (sqrt(j*j + k*k) + 1) , -2);
}

float inverseSquareTemporalCorr(int t, float *ptr_to_norm, float temporal_sigma){
  return pow( (abs(t) + 1) , -2);
}

float inverseCubeSpatialCorr(int j, int k, float *ptr_to_norm, float spatial_sigma){
  return pow( (sqrt(j*j + k*k) + 1) , -3);
}

float inverseCubeTemporalCorr(int t, float *ptr_to_norm, float temporal_sigma){
  return pow( (abs(t) + 1) , -3);
}

float inverseQuarticSpatialCorr(int j, int k, float *ptr_to_norm, float spatial_sigma){
  return pow( (sqrt(j*j + k*k) + 1) , -4);
}

float inverseQuarticTemporalCorr(int t, float *ptr_to_norm, float temporal_sigma){
  return pow( (abs(t) + 1) , -4);
}

float topHatSpatialCorr(int j, int k, float *ptr_to_norm, float spatial_sigma){
  return 0;
}

float topHatTemporalCorr(int t, float *ptr_to_norm, float temporal_sigma){
  return 0;
}

float trueTopHatSpatialCorr(int j, int k, float *ptr_to_norm, float spatial_sigma){
  return 0;
}

float trueTopHatTemporalCorr(int t, float *ptr_to_norm, float temporal_sigma){
  return 0;
}

float trueTopHatRandomSpatialCorr(int j, int k, float *ptr_to_norm, float spatial_sigma){
  return 0;
}

float trueTopHatRandomTemporalCorr(int t, float *ptr_to_norm, float temporal_sigma){
  return 0;
}

float topHatLongTailSpatialCorr(int j, int k, float *ptr_to_norm, float spatial_sigma){
  return 0;
}

float topHatLongTailTemporalCorr(int t, float *ptr_to_norm, float temporal_sigma){
  return 0;
}

float triangleSpatialCorr(int j, int k, float *ptr_to_norm, float spatial_sigma){
  return 0;
}

float triangleTemporalCorr(int t, float *ptr_to_norm, float temporal_sigma){
  return 0;
}

/*------------------------------------------------------------------------*/

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
