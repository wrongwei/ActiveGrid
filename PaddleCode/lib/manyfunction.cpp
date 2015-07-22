/*------------------------------------------------------------------------*/
/* Author: Kevin Griffin                                                  */
/* kevinpg@princeton.edu                                                  */
/* 22/07/15                                                               */
/* manyfunction.cpp                                                       */
/*------------------------------------------------------------------------*/
// Compile with g++ -O3 manyfunction.cpp loaf.cpp -o manyfunction
/*------------------------------------------------------------------------*/
#include <iostream>
#include <cmath>
#include "loaf.h"
/*------------------------------------------------------------------------*/
int main(){
  float newslice[13][11]/* deleted this initialization = {0} */;
  float crumb = 0;
  int norm = 13;
  int timeSigma = 30;
  int spaceSigma = 25;
  int correction = 9;
  enum {range_of_corr = 7};
  enum {halfLoaf = 5000};
  loaf myLoaf(10001);
  for (int col = 0; col < 13; col++) {
    for (int row = 0; row < 11; row++) {
      newslice[col][row] = 0;
      for (int j = -range_of_corr; j <= range_of_corr; j++) {
	for (int k = -range_of_corr; k <= range_of_corr;k++) {
	  for (int t = -halfLoaf; t <= halfLoaf; t++) {
	    crumb = myLoaf.Loaf_access(j, k, (t + halfLoaf));
	    newslice[col][row] += correction * (float) exp( -(pow((float) j, (int) 2) + pow((float) k, (int) 2)) / (2 * pow(spaceSigma, 2))) * (float)exp( -(pow((float) t, 2) / (2 * pow(timeSigma, 2)))) * crumb / norm;
	  }
	}
      }
    }
  }
  //+  myLoaf.Loaf_print();
  return 0;
}
