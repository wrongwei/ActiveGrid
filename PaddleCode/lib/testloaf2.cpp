/*------------------------------------------------------------------------*/
/* Author: Kevin Griffin                                                  */
/* kevinpg@princeton.edu                                                  */
/* 16/07/15                                                               */
/* testloaf2.cpp                                                           */
/*------------------------------------------------------------------------*/

#include "loaf.h"
#include <iostream>

/*------------------------------------------------------------------------*/

// A stand alone test client for loaf.cpp that does not run in menuII
// To compile: g++ testloaf2.cpp -o testloaf2
// To run: ./testloaf2

int main(int argc, char *argv[])
{
  loaf myLoaf = loaf(5);
  cout << "Test client for loaf.cpp is running..." << endl << endl;
  
  // Test Loaf_print
  myLoaf.Loaf_print();
  
  // Test Loaf_set and Loaf_print
  cout << "Testing Loaf_set and Loaf_print..." << endl;
  myLoaf.Loaf_set( 2, 0, 0,      0);
  myLoaf.Loaf_set( 2, 1, 0,      1);
  myLoaf.Loaf_set( 2, 2, 0,      2);
  myLoaf.Loaf_set( 2, 3, 0,      3);
  myLoaf.Loaf_set( 0, 0, 2, -88.88);
  myLoaf.Loaf_set( 0,10, 2,  88.88);
  myLoaf.Loaf_set(12, 0, 2,  48.84);
  myLoaf.Loaf_set(12,10, 2, -48.84);
  myLoaf.Loaf_set(12,10, 4,      3);
  myLoaf.Loaf_set(12, 9, 4,      2);
  myLoaf.Loaf_set(12, 8, 4,      1);
  myLoaf.Loaf_print();
  
  //These calls should trigger assert to fail (causing program to abort)
  //myLoaf.Loaf_set(-1, 1, 1, 30);
  //myLoaf.Loaf_set(13, 1, 1, 30);
  //myLoaf.Loaf_set( 1,-1, 4, 30);
  //myLoaf.Loaf_set( 1,11, 4, 30);
  //myLoaf.Loaf_set( 1, 1,-1, 30);
  //myLoaf.Loaf_set( 1, 1, 5, 30);
  
  // Testing Loaf_access
  cout << "Testing Loaf_access..." << endl;
  cout << "Should print -88.88:   " << myLoaf.Loaf_access(0,0,2) << endl;
  cout << "Should print  88.88:   " << myLoaf.Loaf_access(0,10,2) << endl;
  cout << "Should print  48.88:   " << myLoaf.Loaf_access(12,0,2) << endl;
  cout << "Should print -48.88:   " << myLoaf.Loaf_access(12,10,2) << "\n\n";
  
  // Test Loaf_slice
  cout << "Testing Loaf_slice..." << endl;
  myLoaf.Loaf_slice();
  myLoaf.Loaf_print();
    
  myLoaf.~loaf();
  
  cout << "Done testing loaf.cpp!" << endl << endl;
  return 0;
}
