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
float gaussianSpatialCorr(int j, int k, double dist, float spatial_sigma);
float gaussianTemporalCorr(int t, double dist, float temporal_sigma);
float inverseSquareSpatialCorr(int j, int k, double dist, float spatial_sigma);
float inverseSquareTemporalCorr(int t, double dist, float temporal_sigma);
float inverseCubeSpatialCorr(int j, int k, double dist, float spatial_sigma);
float inverseCubeTemporalCorr(int t, double dist, float temporal_sigma);
float inverseQuarticSpatialCorr(int j, int k, double dist, float spatial_sigma);
float inverseQuarticTemporalCorr(int t, double dist, float temporal_sigma);
float topHatSpatialCorr(int j, int k, double dist, float spatial_sigma);
float topHatTemporalCorr(int t, double dist, float temporal_sigma);
float trueTopHatSpatialCorr(int j, int k, double dist, float spatial_sigma);
float trueTopHatTemporalCorr(int t, double dist, float temporal_sigma);
float trueTopHatRandomSpatialCorr(int j, int k, double dist, float spatial_sigma);
float trueTopHatRandomTemporalCorr(int t, double dist, float temporal_sigma);
float topHatLongTailSpatialCorr(int j, int k, double dist, float spatial_sigma);
float topHatLongTailTemporalCorr(int t, double dist, float temporal_sigma);
float traingleSpatialCorr(int j, int k, double dist, float spatial_sigma);
float triangleTemporalCorr(int t, double dist, float temporal_sigma);

/*------------------------------------------------------------------------*/
// comment
float (*pickSpatialCorr(int typeOfSpatialCorr)) (int j, int k, double dist, float spatial_sigma);
// comment
float (*pickTemporalCorr(int typeOfTemporalCorr)) (int t, double dist, float temporal_sigma);

/*------------------------------------------------------------------------*/
#endif // PICK_CORRELATIONS_INCLUDED
