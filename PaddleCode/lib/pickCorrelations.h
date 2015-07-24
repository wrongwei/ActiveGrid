/*------------------------------------------------------------------------*/
/* Author: Kevin Griffin                                                  */
/* kevinpg@princeton.edu                                                  */
/* 22/07/15                                                               */
/* pickCorrelations.h                                                     */
/*------------------------------------------------------------------------*/

#ifndef PICK_CORRELATIONS_INCLUDED
#define PICK_CORRELATIONS_INCLUDED

/*------------------------------------------------------------------------*/
// The correlation functions
float gaussianSpatialCorr(int j, int k, float *ptr_to_norm, float *ptr_to_norm1, float spatial_sigma, float height);
float gaussianTemporalCorr(int t, float *ptr_to_norm, float *ptr_to_norm1, float temporal_sigma, float height);
float inverseSquareSpatialCorr(int j, int k, float *ptr_to_norm, float *ptr_to_norm1, float spatial_sigma, float height);
float inverseSquareTemporalCorr(int t, float *ptr_to_norm, float *ptr_to_norm1, float temporal_sigma, float height);
float inverseCubeSpatialCorr(int j, int k, float *ptr_to_norm, float *ptr_to_norm1, float spatial_sigma, float height);
float inverseCubeTemporalCorr(int t, float *ptr_to_norm, float *ptr_to_norm1, float temporal_sigma, float height);
float inverseQuarticSpatialCorr(int j, int k, float *ptr_to_norm, float *ptr_to_norm1, float spatial_sigma, float height);
float inverseQuarticTemporalCorr(int t, float *ptr_to_norm, float *ptr_to_norm1, float temporal_sigma, float height);
float topHatSpatialCorr(int j, int k, float *ptr_to_norm, float *ptr_to_norm1, float spatial_sigma, float height);
float topHatTemporalCorr(int t, float *ptr_to_norm, float *ptr_to_norm1, float temporal_sigma, float height);
float trueTopHatSpatialCorr(int j, int k, float *ptr_to_norm, float *ptr_to_norm1, float spatial_sigma, float height);
float trueTopHatTemporalCorr(int t, float *ptr_to_norm, float *ptr_to_norm1, float temporal_sigma, float height);
float trueTopHatRandomSpatialCorr(int j, int k, float *ptr_to_norm, float *ptr_to_norm1, float spatial_sigma, float height);
float trueTopHatRandomTemporalCorr(int t, float *ptr_to_norm, float *ptr_to_norm1, float temporal_sigma, float height);
float topHatLongTailSpatialCorr(int j, int k, float *ptr_to_norm, float *ptr_to_norm1, float spatial_sigma, float height);
float topHatLongTailTemporalCorr(int t, float *ptr_to_norm, float *ptr_to_norm1, float temporal_sigma, float height);
float traingleSpatialCorr(int j, int k, float *ptr_to_norm, float *ptr_to_norm1, float spatial_sigma, float height);
float triangleTemporalCorr(int t, float *ptr_to_norm, float *ptr_to_norm1, float temporal_sigma, float height);

/*------------------------------------------------------------------------*/
// comment
float (*pickSpatialCorr(int typeOfSpatialCorr)) (int j, int k, float *ptr_to_norm, float *ptr_to_norm1, float spatial_sigma, float height);
// comment
float (*pickTemporalCorr(int typeOfTemporalCorr)) (int t, float *ptr_to_norm, float *ptr_to_norm1, float temporal_sigma, float height);

/*------------------------------------------------------------------------*/
#endif // PICK_CORRELATIONS_INCLUDED
