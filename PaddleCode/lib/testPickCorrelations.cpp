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
/* For the details of what the functions in pickCorrelations.cpp do, see the comments in that file.
 * This file will help you test those functions and teach you how to use them 
 */
/*------------------------------------------------------------------------*/
int main(){
  // Declaring some variables that will be used as arguments of the correlation functions
  int j = 2;
  int k = 3;
  int t = 10;
  int spatial_sigma = 5;
  int temporal_sigma = 5;
  float norm = 1.2;
  float norm1 = 1.2;
  int height = 3;
  // Declare a pointer to norm because we want the correlation functions to be
  // able to incrment norm if they need to (I think some top hat functions need
  // this feature
  float *ptr_to_norm = &norm; // NOTE: ptr_to_norm and &norm are both float pointers to the same float (named norm)
  float *ptr_to_norm1 = &norm1; // NOTE: ptr_to_norm and &norm are both float pointers to the same float (named norm)

  // Decalre function pointers for the spacial and temporal correlation functions
  // These are variables that you will use to store the name of a function that
  // pickTemporalCorr or pickSpatialCorr returns.
  float (*pfSpatialCorr)(int j, int k, float *ptr_to_norm, float *ptr_to_norm1, float spatial_sigma, float height); // a function pointer named pfSpatialCorr
  float (*pfTemporalCorr)(int t, float *ptr_to_norm, float *ptr_to_norm1, float temporal_sigma, float height); // a function pointer named pfTemporalCorr
  
  // Now let's assign values to these functions pointers so that they actually point somewhere
  pfSpatialCorr = pickSpatialCorr(1);
  // I have picked 1(which corresponds to gaussian), so now pfSpatialCorr = the name of the function gaussianSpatialCorr
  // So I can use pfSpatialCorr in the same way that I can use gaussianSpatialCorr
  // For example,
  cout << "Test if the pointers are working..." << endl;
  cout << "This value: " << gaussianSpatialCorr(j, k, ptr_to_norm, ptr_to_norm1, spatial_sigma, height) << endl;
  cout << "Should be equivalent to this value: " << pfSpatialCorr(j, k, ptr_to_norm, ptr_to_norm1, spatial_sigma, height) << endl;

  // If you understand the above example, then you are ready to use this module!

  // The loops below just allow us to test that all the correlation functions are working properly
  cout << endl << "Testing Temporal Correlations..." << endl;

  for (int typeOfTemporalCorr = 1; typeOfTemporalCorr <= 9; typeOfTemporalCorr++){

    // initialize these function pointers
    pfTemporalCorr = pickTemporalCorr(typeOfTemporalCorr);
  
    //Call the correlation functions using the previously initialized function pointers
    cout << "Type of temporal correlation = " << typeOfTemporalCorr;
    cout << "   Value of the correlation = " << pfTemporalCorr(t, ptr_to_norm, ptr_to_norm1, temporal_sigma, height) << endl;
  }

  cout << endl << "Testing Spatial Correlations..." << endl;

  for (int typeOfSpatialCorr = 1; typeOfSpatialCorr <= 9; typeOfSpatialCorr++){

    // initialize these function pointers
    pfSpatialCorr = pickSpatialCorr(typeOfSpatialCorr);
  
    //Call the correlation functions using the previously initialized function pointers
    cout << "Type of spatial correlation = " << typeOfSpatialCorr;
    cout << "   Value of the correlation = " << pfSpatialCorr(j, k, ptr_to_norm, ptr_to_norm1, spatial_sigma, height) << endl;
  }
  
  cout << endl;
  
  return 0;
}
