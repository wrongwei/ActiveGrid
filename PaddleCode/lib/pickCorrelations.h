/*------------------------------------------------------------------------*/
/* Author: Kevin Griffin                                                  */
/* kevinpg@princeton.edu                                                  */
/* 22/07/15                                                               */
/* pickCorrelations.h                                                     */
/*------------------------------------------------------------------------*/

#ifndef PICK_CORRELATIONS_INCLUDED
#define PICK_CORRELATIONS_INCLUDED

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

// The spatial correlation fuctions. pickSpatialCorr returns
// pointers to these functions.
float gaussianSpatialCorr(int j, int k, float *ptr_to_norm, float spatial_sigma);
float inverseSquareSpatialCorr(int j, int k, float *ptr_to_norm, float spatial_sigma);
float inverseCubeSpatialCorr(int j, int k, float *ptr_to_norm, float spatial_sigma);
float inverseQuarticSpatialCorr(int j, int k, float *ptr_to_norm, float spatial_sigma);
float topHatSpatialCorr(int j, int k, float *ptr_to_norm, float spatial_sigma);
float trueTopHatSpatialCorr(int j, int k, float *ptr_to_norm, float spatial_sigma);
float trueTopHatRandomSpatialCorr(int j, int k, float *ptr_to_norm, float spatial_sigma);
float topHatLongTailSpatialCorr(int j, int k, float *ptr_to_norm, float spatial_sigma);
float traingleSpatialCorr(int j, int k, float *ptr_to_norm, float spatial_sigma);

// The temporal correlation fuctions. pickTemporalCorr returns
// pointers to these functions.
float gaussianTemporalCorr(int t, float *ptr_to_norm, float temporal_sigma);
float inverseSquareTemporalCorr(int t, float *ptr_to_norm, float temporal_sigma);
float inverseCubeTemporalCorr(int t, float *ptr_to_norm, float temporal_sigma);
float inverseQuarticTemporalCorr(int t, float *ptr_to_norm, float temporal_sigma);
float topHatTemporalCorr(int t, float *ptr_to_norm, float temporal_sigma);
float trueTopHatTemporalCorr(int t, float *ptr_to_norm, float temporal_sigma);
float trueTopHatRandomTemporalCorr(int t, float *ptr_to_norm, float temporal_sigma);
float topHatLongTailTemporalCorr(int t, float *ptr_to_norm, float temporal_sigma);
float triangleTemporalCorr(int t, float *ptr_to_norm, float temporal_sigma);

/*------------------------------------------------------------------------*/
/* The function pickSpatialCorr accepts an integer named typeOfSpatialCorr
 * and returns a pointer to a function (which accepts an int, int, float*,float as parameters
 * and returns a float)
 * In other words, this function returns a correlation function.
 */
float (*pickSpatialCorr(int typeOfSpatialCorr)) (int j, int k, float *ptr_to_norm, float spatial_sigma);

/* The function pickTemporalCorr accepts an integer named typeOfTemporalCorr
 * and returns a pointer to a function (which accepts an int, int, float*,float as parameters
 * and returns a float)
 * In other words, this function returns a correlation function.
 */
float (*pickTemporalCorr(int typeOfTemporalCorr)) (int t, float *ptr_to_norm, float temporal_sigma);

/*------------------------------------------------------------------------*/
#endif // PICK_CORRELATIONS_INCLUDED
