/*------------------------------------------------------------------------*/
/* Author: Kevin Griffin                                                  */
/* kevinpg@princeton.edu                                                  */
/* 16/07/15                                                               */
/* loaf.h                                                                 */
/*------------------------------------------------------------------------*/

// Protecting against multiple inclusion ...
#ifndef LOAF_INCLUDED
#define LOAF_INCLUDED

/*------------------------------------------------------------------------*/

#include <stddef.h>
#include <iostream>
using namespace std;

/*------------------------------------------------------------------------*/

/* A Loaf_T is a data structure for storing the servo positions of multiple
 * time steps. The servo positions in a particular time step are called a slice.
 * A Loaf_T object contains multiple slices. 
 * This structure is used by runcorr_3D to do correlations in space and time.
 * Multiple Loaf_T 's will have distinct values and exisit in different
 * locations in memory. So remember, everytime you create a Loaf_T object,
 * you must destroy it when you are done using it (see the method Loaf_eat). 
 * The Loaf_T object contains slices. Each slice is an array of floats.
 * so a Loaf_T is pointer to an array of floats. Thus a Loaf_T is a float**
 */

typedef float **Loaf_T;

/*------------------------------------------------------------------------*/

class loaf
{
 public:
  //constructor:
  loaf(int numberOfSlices);
  // destructor:
  ~loaf() {};

/*------------------------------------------------------------------------*/
// functions...
/*------------------------------------------------------------------------*/

/* Returns a new Loaf_T object that contains numberOfSlices randomized slices
 * Each slice is a 2D array of the angles of all the servos.
 * These 2D arrays are 13 x 11 */

//Loaf_T Loaf_bake(int numberOfSlices);

/*------------------------------------------------------------------------*/

/* Returns the angle value of the servo in myLoaf in column i,
 * in row j, in the time slice t
 * i goes from 0 to 12
 * j goes from 0 to 10
 * t = 0 for the oldest slice and t = (numberOfSlices - 1) in the newest slice.
 */

  float Loaf_access(/*Loaf_T myLoaf, */int i, int j, int t);

/*------------------------------------------------------------------------*/

/* Sets the angle value of the servo in myLoaf in column i,
 * in row j, in the time slice t
 * i goes from 1 to 13
 * j goes from 1 to 11
 * t = 0 for the oldest slice and t = (numberOfSlices - 1) in the newest slice.
 * Returns nothing.
 * NOTE: this method is only for testing the LOAF_T data type.
 *       RUNCORR_3D SHOULD NEVER CALL Loaf_set
 */

void Loaf_set(/*Loaf_T myLoaf, */int i, int j, int t, float newAngle);

/*------------------------------------------------------------------------*/

/* Frees the oldest slice in myLoaf and adds a new slice. The new slice contains
 * random angle positions. Returns nothing */

 void Loaf_slice(/*Loaf_T myLoaf*/);

/*------------------------------------------------------------------------*/

/* Prints out the contents of myLoaf. Returns nothing */

void Loaf_print(/*Loaf_T myLoaf*/);

/*------------------------------------------------------------------------*/

/* Loaf_eat frees all memory used by myLoaf and returns nothing */

// void Loaf_eat(/*Loaf_T myLoaf*/);

/*------------------------------------------------------------------------*/
//  variables...
/*------------------------------------------------------------------------*/
 Loaf_T myLoaf;      // a new Loaf_T object

 private:

 // The number of slices is initialized by the constructor
 int NUMBER_OF_SLICES;
 
 // max and min servo angles
 //If float values are desired, replace the enum statements with the commented const float ones
 enum {MAX_ANGLE = 90};
 // const float MAX_ANGLE = 90.0;
 enum {MIN_ANGLE = -90};
 // const float MIN_ANGLE = -90.0;
};

#endif // LOAF_INCLUDED
