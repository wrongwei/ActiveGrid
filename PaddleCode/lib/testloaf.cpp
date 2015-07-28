/*------------------------------------------------------------------------*/
/* Author: Kevin Griffin                                                  */
/* kevinpg@princeton.edu                                                  */
/* 16/07/15                                                               */
/* testloaf.cpp                                                           */
/*------------------------------------------------------------------------*/

#include "loaf.h"
#include <iostream>

/*------------------------------------------------------------------------*/

// A test client for loaf.cpp
// To run this test client, run menuII and enter the corresponding number

void testloaf(void)
{
  loaf myLoaf = loaf(5);
  cout << "Test client for loaf.cpp is running..." << endl << endl;
  
  // Test Loaf_print
  myLoaf.Loaf_print();
  
  // Test Loaf_set and Loaf_print
  cout << "Testing Loaf_set and Loaf_print..." << endl;
  myLoaf.Loaf_set( 2, 0, 0,  77.7);
  myLoaf.Loaf_set( 2, 1, 0,  77.7);
  myLoaf.Loaf_set( 2, 2, 0,  77.7);
  myLoaf.Loaf_set( 2, 3, 0,     3);
  myLoaf.Loaf_set( 0, 0, 2, -12.3);
  myLoaf.Loaf_set( 0,10, 2,  12.3);
  myLoaf.Loaf_set(12, 0, 2,  32.1);
  myLoaf.Loaf_set(12,10, 2, -32.1);
  myLoaf.Loaf_set(12,10, 4,     3);
  myLoaf.Loaf_set(12, 9, 4,     2);
  myLoaf.Loaf_set(12, 8, 4,     1);
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
  cout << "Should print -12.3:   " << myLoaf.Loaf_access(0,0,2) << endl;
  cout << "Should print  12.3:   " << myLoaf.Loaf_access(0,10,2) << endl;
  cout << "Should print  32.1:   " << myLoaf.Loaf_access(12,0,2) << endl;
  cout << "Should print -32.1:   " << myLoaf.Loaf_access(12,10,2) << "\n\n";
  
  // Test Loaf_slice
  cout << "Testing Loaf_slice..." << endl;
  myLoaf.Loaf_slice();
  myLoaf.Loaf_print();

  // Testing Loaf_access after a slice
  cout << "Testing Loaf_access after a slice..." << endl;
  cout << "Should print -12.3:   " << myLoaf.Loaf_access(0,0,1) << endl;
  cout << "Should print  12.3:   " << myLoaf.Loaf_access(0,10,1) << endl;
  cout << "Should print  32.1:   " << myLoaf.Loaf_access(12,0,1) << endl;
  cout << "Should print -32.1:   " << myLoaf.Loaf_access(12,10,1) << "\n\n";
    
  // Testing Loaf_set and Loaf_printFullArray
  cout << "Testing Loaf_set and Loaf_printFullArray..." << endl;
  myLoaf.Loaf_set(-7, -7, 1, -12.3);
  myLoaf.Loaf_set(-7, 17, 1,  12.3);
  myLoaf.Loaf_set(19, -7, 1,  32.1);
  myLoaf.Loaf_set(19, 17, 1, -32.1);

  myLoaf.Loaf_printFullArray();
  myLoaf.Loaf_printLowLevelMirror();

  myLoaf.~loaf();
  
  cout << "Done testing loaf.cpp!" << endl << endl;
}
