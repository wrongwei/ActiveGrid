/*------------------------------------------------------------------------*/
/* Author: Kevin Griffin                                                  */
/* kevinpg@princeton.edu                                                  */
/* 22/07/15                                                               */
/* testPickCorrelations.cpp                                               */
/*------------------------------------------------------------------------*/

// This is a test method for pickCorrelations.cpp

/*------------------------------------------------------------------------*/

// Compile with g++ testPickCorrelations.cpp pickCorrelations.cpp -o testPickCorrelations

/*------------------------------------------------------------------------*/
#include <iostream>
#include "pickCorrelations.h"
using namespace std;
/*------------------------------------------------------------------------*/
int main(){
  int j = 2;
  int k = 3;
  int t = 10;
  int spatial_sigma = 5;
  int temporal_sigma = 5;
  float norm = 1.2;
  float *ptr_to_norm = &norm;

  // Decalre function pointers for the spacial and temporal correlation functions
  float (*pfSpatialCorr)(int j, int k, float *ptr_to_norm, float spatial_sigma);
  float (*pfTemporalCorr)(int t, float *ptr_to_norm, float temporal_sigma);
  
  cout << endl << "Testing Temporal Correlations..." << endl;

  for (int typeOfTemporalCorr = 1; typeOfTemporalCorr <= 9; typeOfTemporalCorr++){

    // initialize these function pointers
    pfTemporalCorr = pickTemporalCorr(typeOfTemporalCorr);
  
    //Call the correlation functions using the previously initialized function pointers
    cout << "Type of temporal correlation = " << typeOfTemporalCorr;
    cout << "   Value of the correlation = " << pfTemporalCorr(t, ptr_to_norm, temporal_sigma) << endl;
  }

  cout << endl << "Testing Spatial Correlations..." << endl;

  for (int typeOfSpatialCorr = 1; typeOfSpatialCorr <= 9; typeOfSpatialCorr++){

    // initialize these function pointers
    pfSpatialCorr = pickSpatialCorr(typeOfSpatialCorr);
  
    //Call the correlation functions using the previously initialized function pointers
    cout << "Type of spatial correlation = " << typeOfSpatialCorr;
    cout << "   Value of the correlation = " << pfSpatialCorr(j, k, ptr_to_norm, spatial_sigma) << endl;
  }
  
  cout << endl;
  
  return 0;
}
