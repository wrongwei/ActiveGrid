/*------------------------------------------------------------------------*/
/* Author: Kevin Griffin                                                  */
/* kevinpg@princeton.edu                                                  */
/* 16/07/15                                                               */
/* testloaf.cpp                                                           */
/*------------------------------------------------------------------------*/

#include "loaf.h"
#include <iostream>

/*------------------------------------------------------------------------*/

/* A test client for loaf.cpp */

int main(void)
{
  cout << "Test client for loaf.cpp is running..." << endl << endl;
  
  Loaf_T myLoaf = Loaf_bake(5);

  // Test Loaf_print
  Loaf_print(myLoaf);

  // Test Loaf_set and Loaf_print
  Loaf_set(myLoaf,  2, 0, 0,      0);
  Loaf_set(myLoaf,  2, 1, 0,      1);
  Loaf_set(myLoaf,  2, 2, 0,      2);
  Loaf_set(myLoaf,  2, 3, 0,      3);
  Loaf_set(myLoaf,  0, 0, 2, -88.88);
  Loaf_set(myLoaf,  0,10, 2,  88.88);
  Loaf_set(myLoaf, 12, 0, 2,  48.84);
  Loaf_set(myLoaf, 12,10, 2, -48.84);
  Loaf_set(myLoaf, 12,10, 4,      3);
  Loaf_set(myLoaf, 12, 9, 4,      2);
  Loaf_set(myLoaf, 12, 8, 4,      1);
  Loaf_print(myLoaf);

  //These calls should trigger assert to fail (causing program to abort)
  //Loaf_set(myLoaf, -1, 1, 1,      30);
  //Loaf_set(myLoaf, 13, 1, 1,      30);
  //Loaf_set(myLoaf,  1,-1, 4,      30);
  //Loaf_set(myLoaf,  1,11, 4,      30);
  //Loaf_set(myLoaf,  1, 1,-1,      30);
  //Loaf_set(myLoaf,  1, 1, 5,      30);

  // Test Loaf_access
  cout << Loaf_access(myLoaf, 0,0,2) << endl;
  cout << Loaf_access(myLoaf, 0,10,2) << endl;
  cout << Loaf_access(myLoaf, 12,0,2) << endl;
  cout << Loaf_access(myLoaf, 12,10,2) << endl;

  // Test Loaf_slice
  Loaf_slice(myLoaf);
  Loaf_print(myLoaf);
  
  Loaf_eat(myLoaf);

  cout << "Done testing client for loaf.cpp!" << endl;
  return 0;
}
