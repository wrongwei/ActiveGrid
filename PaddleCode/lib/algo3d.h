/********************* ALGO3D HEADER FILE *********************
 
 Includes methods and functions for controlling grid correlations
   in the temporal dimension (hence "3D")
 Referenced by menu3d.cpp
 Written by Kevin Griffin and Nathan Wei, July-August 2015
 
 **************************************************************/

#ifndef ALGO3D_INCLUDED
#define ALGO3D_INCLUDED

#include <iostream>
#include <cstdlib>
#include <time.h>
#include <cmath>
#include <fstream>
#include "activegrid.h"
#include "loaf.h"
#include "pickCorrelations.h"

#define range_of_corr 7

//declare the testloaf function
void testloaf(void);

// methods in algo3d.cpp

/* Because these are static methods (i.e. private methods for use in algo3d.cpp only), they should not be
   declared in this algo3d.h file. If you want to make them public, uncomment these lines.
static int setanglestoallservosIII(float angles[13][11], float steps[13][11], int constant, float rms);

static void runcorr_3D(float newslice[][11], loaf* myLoaf, int halfLoaf, float spaceSigma, float timeSigma,
		       float spaceAlpha, float timeAlpha, float spaceHeight, float timeHeight, int spaceMode,
		       int timeMode, int mrow, int mcol, float correction);

static void ltfast(float newslice[][11], loaf* myLoaf, int halfLoaf, float spaceSigma, float timeSigma,
		   float spaceAlpha, float timeAlpha, float spaceHeight, float timeHeight, float correction);

static void unsharp(float newslice[][11], loaf* myLoaf, int halfLoaf, float spaceSigma, float timeSigma,
		    float spaceAlpha, float timeAlpha, float spaceHeight, float timeHeight, float correction);

static float compute_rmscorr_3D(float spaceSigma, float timeSigma, int spaceMode, int timeMode, float spaceAlpha,
				float timeAlpha, float spaceHeight, float timeHeight, int mrow, int mcol, int halfLoaf);
 */

int correlatedMovement_correlatedInTime(int constantArea, float spatial_sigma, float temporal_sigma, float spatial_alpha,
					float temporal_alpha, float spatial_height, float temporal_height,
					int typeOfSpatialCorr, int typeOfTemporalCorr, float target_rms, int numberOfSlices);

#endif // ALGO3D_INCLUDED
