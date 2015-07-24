/*------------------------------------------------------------------------*/
/* Author: Kevin Griffin                                                  */
/* kevinpg@princeton.edu                                                  */
/* 22/07/15                                                               */
/* pickCorrelations.h                                                     */
/*------------------------------------------------------------------------*/

#ifndef PICK_CORRELATIONS_INCLUDED
#define PICK_CORRELATIONS_INCLUDED
#include <iostream>
using namespace std;
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

// The spatial correlation fuctions. pickSpatialCorr returns
// pointers to these functions.
float gaussianSpatialCorr(int j, int k, float *ptr_to_norm, float *ptr_to_norm1, float spatial_sigma, float height);
float inverseSquareSpatialCorr(int j, int k, float *ptr_to_norm, float *ptr_to_norm1, float spatial_sigma, float height);
float inverseCubeSpatialCorr(int j, int k, float *ptr_to_norm, float *ptr_to_norm1, float spatial_sigma, float height);
float inverseQuarticSpatialCorr(int j, int k, float *ptr_to_norm, float *ptr_to_norm1, float spatial_sigma, float height);
float topHatSpatialCorr(int j, int k, float *ptr_to_norm, float *ptr_to_norm1, float spatial_sigma, float height);
float trueTopHatSpatialCorr(int j, int k, float *ptr_to_norm, float *ptr_to_norm1, float spatial_sigma, float height);
float trueTopHatRandomSpatialCorr(int j, int k, float *ptr_to_norm, float *ptr_to_norm1, float spatial_sigma, float height);
float topHatLongTailSpatialCorr(int j, int k, float *ptr_to_norm, float *ptr_to_norm1, float spatial_sigma, float height);
float triangleSpatialCorr(int j, int k, float *ptr_to_norm, float *ptr_to_norm1, float spatial_sigma, float height);

// The temporal correlation fuctions. pickTemporalCorr returns
// pointers to these functions.
float gaussianTemporalCorr(int t, float *ptr_to_norm, float *ptr_to_norm1, float temporal_sigma, float height);
float inverseSquareTemporalCorr(int t, float *ptr_to_norm, float *ptr_to_norm1, float temporal_sigma, float height);
float inverseCubeTemporalCorr(int t, float *ptr_to_norm, float *ptr_to_norm1, float temporal_sigma, float height);
float inverseQuarticTemporalCorr(int t, float *ptr_to_norm, float *ptr_to_norm1, float temporal_sigma, float height);
float topHatTemporalCorr(int t, float *ptr_to_norm, float *ptr_to_norm1, float temporal_sigma, float height);
float trueTopHatTemporalCorr(int t, float *ptr_to_norm, float *ptr_to_norm1, float temporal_sigma, float height);
float trueTopHatRandomTemporalCorr(int t, float *ptr_to_norm, float *ptr_to_norm1, float temporal_sigma, float height);
float topHatLongTailTemporalCorr(int t, float *ptr_to_norm, float *ptr_to_norm1, float temporal_sigma, float height);
float triangleTemporalCorr(int t, float *ptr_to_norm, float *ptr_to_norm1, float temporal_sigma, float height);

/*------------------------------------------------------------------------*/
/* The function pickSpatialCorr accepts an integer named typeOfSpatialCorr
 * and returns a pointer to a function (which accepts an int, int, float*,float as parameters
 * and returns a float)
 * In other words, this function returns a correlation function.
 */
float (*pickSpatialCorr(int typeOfSpatialCorr)) (int j, int k, float *ptr_to_norm, float *ptr_to_norm1, float spatial_sigma, float height);

/* The function pickTemporalCorr accepts an integer named typeOfTemporalCorr
 * and returns a pointer to a function (which accepts an int, int, float*,float as parameters
 * and returns a float)
 * In other words, this function returns a correlation function.
 */
float (*pickTemporalCorr(int typeOfTemporalCorr)) (int t, float *ptr_to_norm, float *ptr_to_norm1, float temporal_sigma, float height);
/*------------------------------------------------------------------------*/
#endif // PICK_CORRELATIONS_INCLUDED
