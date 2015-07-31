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

// random. In other words, no correlation
// We have not implemented randomTemporalCorr in menuII because this would lead to many angles
// being asked to move more than 40 degrees per sec
float randomSpatialCorr(int j, int k, float spatial_sigma, float height){
    if (j == 0 && k == 0)
	return 1;
    return 0;
}

float randomTemporalCorr(int t, float temporal_sigma, float height){
    if (t == 0)
	return 1;
    return 0;
}

// Gaussian convolution function
float gaussianSpatialCorr(int j, int k, float spatial_sigma, float height){
    return expf(-(((j*j) + (k*k)) / (2 * spatial_sigma*spatial_sigma)));
}

float gaussianTemporalCorr(int t, float temporal_sigma, float height){
    return expf(-((t*t) / (2 * temporal_sigma*temporal_sigma)));
}

// 1/r^2 convolution function
float inverseSquareSpatialCorr(int j, int k, float spatial_sigma, float height){
    return pow( (sqrt(j*j + k*k) + 1) , -2);
}

float inverseSquareTemporalCorr(int t, float temporal_sigma, float height){
    return pow( (abs(t) + 1) , -2);
}

// 1/r^3 convolution function
float inverseCubeSpatialCorr(int j, int k, float spatial_sigma, float height){
    return pow( (sqrt(j*j + k*k) + 1) , -3);
}

float inverseCubeTemporalCorr(int t, float temporal_sigma, float height){
    return pow( (abs(t) + 1) , -3);
}

// 1/r^4 convolution function
float inverseQuarticSpatialCorr(int j, int k, float spatial_sigma, float height){
    return pow( (sqrt(j*j + k*k) + 1) , -4);
}

float inverseQuarticTemporalCorr(int t, float temporal_sigma, float height){
    return pow( (abs(t) + 1) , -4);
}

// top hat convolution function
float topHatSpatialCorr(int j, int k, float spatial_sigma, float height){
    double dist = sqrt((j*j) + (k*k));
    if (dist <= spatial_sigma) return 1.0;
    else return 0;
}

float topHatTemporalCorr(int t, float temporal_sigma, float height){
    if (abs(t) <= temporal_sigma) return 1.0;
    else return 0;
}

// true top hat with one main paddle, no wrapping around
float trueTopHatSpatialCorr(int j, int k, float spatial_sigma, float height){
    return 0;
}

float trueTopHatTemporalCorr(int t, float temporal_sigma, float height){
    return 0;
}

// true top hat with one randomly chosen paddle
float trueTopHatRandomSpatialCorr(int j, int k, float spatial_sigma, float height){
    return 0;
}

float trueTopHatRandomTemporalCorr(int t, float temporal_sigma, float height){
    return 0;
}

float topHatLongTailSpatialCorr(int j, int k, float alpha, float height){
    double dist = sqrt((j*j)+(k*k));
    if (dist <= alpha) return 1.0;
    else return height;
}

float topHatLongTailTemporalCorr(int t, float alpha, float height){
    if (abs(t) <= alpha) return 1.0;
    else return height;
}

float triangleSpatialCorr(int j, int k, float spatial_sigma, float height){
    double dist = sqrt((j*j)+(k*k));
    if (dist <= spatial_sigma) return ((-1) / spatial_sigma * dist + 1);
    else return 0;
}

float triangleTemporalCorr(int t, float temporal_sigma, float height){
    if (abs(t) <= temporal_sigma) return ((-1) / temporal_sigma * abs(t) + 1);
    else return 0;
}

// unsharp correlations. Unsharp correlations are equal to the image in the
// center and are equal to a negative gaussian everywhere else. 
// The effect for this correlation is that the grid is made less correlated,
// because panels will be made more different to the neighbors; if you think of it
// as a form of image processing it makes the image sharper, which is the opposite
// of the gaussian (which blurs the image)
// Height is the scaling sharpness. Higher height means a sharper image (or sharper contrast)
float unsharpSpatialCorr(int j, int k, float spatial_sigma, float height){
    if (j == 0 && k == 0){
	return 225; //(7*2+1)*(7*2+1) width of spatial ker squared
    }
    return height * -expf(-(((j*j) + (k*k)) / (2 * spatial_sigma*spatial_sigma)));
}

float unsharpTemporalCorr(int t, float temporal_sigma, float height){
    if (t == 0){
	int widthOfTempKer = ceil(6*temporal_sigma)+1;
	if (widthOfTempKer % 2)
	    widthOfTempKer++;
	return widthOfTempKer;
    }
    return height * -expf(-((t*t) / (2 * temporal_sigma*temporal_sigma)));
}

/*------------------------------------------------------------------------*/

float (*pickSpatialCorr(int typeOfSpatialCorr)) (int j, int k, float spatial_sigma, float height){
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
 * and returns a pointer to a function (which accepts an int, int, float*,float as parameters
 * and returns a float)
 * In other words, this function returns a correlation function.
 */
float (*pickTemporalCorr(int typeOfTemporalCorr)) (int t, float temporal_sigma, float height){
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
