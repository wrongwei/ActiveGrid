/*------------------------------------------------------------------------*/
/* Author: Kevin Griffin                                                  */
/* kevinpg@princeton.edu                                                  */
/* 16/07/15                                                               */
/* loaf.cpp                                                               */
/*------------------------------------------------------------------------*/

#include "loaf.h"
#include <stddef.h>
#include <stdlib.h>
#include <iostream>
#include <stdio.h>
#include <assert.h>

/*------------------------------------------------------------------------*/

// number of servos
enum {NUMBER_OF_SERVOS = 143};

/*------------------------------------------------------------------------*/

/* High Level Comment:
 * A Loaf_T is a data structure for storing the servo positions of multiple
 * time steps. The servo positions in a particular time step are called a slice.
 * A Loaf_T object contains multiple slices. 
 * This structure is used by runcorr_3D to do correlations in space and time.
 * Multiple Loaf_T 's will have distinct values and exisit in different
 * locations in memory. So remember, everytime you create a Loaf_T object,
 * you must destroy it when you are done using it (see the method Loaf_eat). 
 * The Loaf_T object contains slices. Each slice is an array of floats.
 * so a Loaf_T is pointer to an array of floats. Thus a Loaf_T is a float**
 */

/* Low Level Comment: 
 * How the Loaf_T works:
 * The Loaf_T is a float**
 * It points to the beginning of an array of float*
 * These point to floats that represent the servo positions at certain time step
 */

/* Assumptions:
 * 1D servo array starts with zero and ends with 142
 * Servos in a 2D array are referred to by order pair (column,row)
 * starting with 0,0 in the bottom right (looking upstream of the wind tunnel)
 * 1,0 is immediately to the right of 0,0
 */

/*------------------------------------------------------------------------*/

/* Returns a new Loaf_T object that contains numberOfSlices randomized slices
 * Each slice is a 2D array of the angles of all the servos.
 * These 2D arrays are 13 x 11 */
loaf::loaf(int numberOfSlices)
	  /*Loaf_T loaf::Loaf_bake(int numberOfSlices)*/
{
  //  Loaf_T myLoaf;      // a new Loaf_T object
  int iSlice;         // counter variable
  int iServo;         // counter variable
  
  // store the number of slices in a global variable so that other
  // methods won't have to keep asking for the number of slices
  NUMBER_OF_SLICES = numberOfSlices;

  // Allocate memory for the array containing the names of slices
  // So this is an array of float*
  myLoaf = (Loaf_T)malloc(NUMBER_OF_SLICES * sizeof(float*));
  // Allocate memory for the slices themselves
  // So these are many float* stored in the array allocated above above.
  for (iSlice =0; iSlice < NUMBER_OF_SLICES; iSlice++)
    {
      myLoaf[iSlice] = (float*)malloc(NUMBER_OF_SERVOS * sizeof(float));
      // Fill the current slice with random angles
      for (iServo = 0; iServo < NUMBER_OF_SERVOS; iServo++)
	{
	  myLoaf[iSlice][iServo] = ((float)rand()/RAND_MAX)*
	    ((float)MAX_ANGLE-(float)MIN_ANGLE) + (float)MIN_ANGLE;
	}
    }
  //  return myLoaf;
}


/*------------------------------------------------------------------------*/

/* Returns the angle value of the servo in myLoaf in column i,
 * in row j, in the time slice t
 * i goes from 0 to 12
 * j goes from 0 to 10
 * t = 0 for the oldest slice and t = (NUMBER_OF_SLICES - 1) in the newest slice.
 */

float loaf::Loaf_access(/*Loaf_T myLoaf, */int i, int j, int t)
{
  // validate parameters
  // validation of parameters is slow, so commenting out the assert statements
  // may speed up the code consideraby
  assert(i <= 12);
  assert(i >= 0);
  assert(j <= 10);
  assert(j >= 0);
  assert(t < NUMBER_OF_SLICES);
  assert(t >= 0);

  return myLoaf[t][j*13 + i];
}

/*------------------------------------------------------------------------*/

/* Sets the angle value of the servo in myLoaf in column i,
 * in row j, in the time slice t
 * i goes from 0 to 12
 * j goes from 0 to 10
 * t = 0 for the oldest slice and t = (NUMBER_OF_SLICES - 1) in the newest slice.
 * Returns nothing.
 * NOTE: this method is only for testing the LOAF_T data type.
 *       RUNCORR_3D SHOULD NEVER CALL Loaf_set
 */

void loaf::Loaf_set(/*Loaf_T myLoaf, */int i, int j, int t, float newAngle)
{
  // validate parameters
  // validation of parameters is slow, so commenting out the assert statements
  // may speed up the code consideraby
  assert(i <= 12);
  assert(i >= 0);
  assert(j <= 10);
  assert(j >= 0);
  assert(t < NUMBER_OF_SLICES);
  assert(t >= 0);

  myLoaf[t][j*13 + i] = newAngle;
}

/*------------------------------------------------------------------------*/

/* Frees the oldest slice in myLoaf and adds a new slice. The new slice contains
 * random angle positions. Returns nothing */

void loaf::Loaf_slice(/*Loaf_T myLoaf*/)
{
  float* newSlice;
  int i;
  int iServo;

  // Allocate memory for a new slice and fill it with random numbers
  newSlice = (float*)malloc(NUMBER_OF_SERVOS * sizeof(float));
  for (iServo = 0; iServo < NUMBER_OF_SERVOS; iServo++)
    {
      newSlice[iServo] = ((float)rand()/RAND_MAX)*
	((float)MAX_ANGLE-(float)MIN_ANGLE) + (float)MIN_ANGLE;
    }
  
  // Free the memory used by the oldest slice
  free(myLoaf[0]);

  // Update the pointers in the loaf so that the previous 2nd oldest is now the
  // 1st oldest slice, the 3rd is now the 2nd oldest, etc.
  for(i = 0; i < (NUMBER_OF_SLICES - 1); i++)
    {
      myLoaf[i]=myLoaf[i + 1];
    }
  
  // Add the new slice to the loaf
  myLoaf[NUMBER_OF_SLICES - 1] = newSlice;
}

/*------------------------------------------------------------------------*/

/* Prints out the contents of myLoaf. Returns nothing */

void loaf::Loaf_print(/*Loaf_T myLoaf*/)
{
  int iSlice;
  int iRow;
  int iColumn;

  for (iSlice = 0; iSlice < NUMBER_OF_SLICES; iSlice ++)
    {
      cout << "Time Slice " << iSlice << endl << endl;
      for (iRow = 10; iRow >= 0; iRow--)
	{
	  printf(" %2i| ", iRow);
	  for (iColumn = 0; iColumn < 13; iColumn++)
	    {
	      printf("%6.2f  ", Loaf_access(/*myLoaf, */iColumn, iRow, iSlice));
	    }
	  cout << endl;
	}
      cout << "    --------------------------------------------------"
	   << "-----------------------------------------------------\n";
      for (iColumn = 0; iColumn < 13; iColumn++)
	printf("%8i", iColumn);
      cout << endl;
    }
}

/*------------------------------------------------------------------------*/

/* Loaf_eat frees all memory used by myLoaf and returns nothing */

//loaf::~loaf()
/*void loaf::Loaf_eat(Loaf_T myLoaf)*/
/*{
  int iSlice;
  for (iSlice = 0; iSlice < NUMBER_OF_SLICES; iSlice++)
    free(myLoaf[iSlice]);
  free(myLoaf);
}*/