MATLAB Code Directory for Active Grid Data Processing
readme.txt created: 27 August 2015 by Nathan Wei
readme.txt last updated: 27 August 2015 by Nathan Wei
Scripts written by Greg Bewley (2009-2011), Horace Zhang and Jessie Liu (2014), Nathan Wei and Kevin Griffin (2015)

Note: Many of these scripts require specific directory setups, such as specific folder name formats and data or calibration files in specific locations. This was done to make data processing simpler and less work-intensive for the user. Please refer to the headers and comments on each script to see how your files and data should be organized for the programs to run correctly.


****************************** Code Architecture Outline ******************************
This is a list of the most useful scripts and their dependencies.
Bullets and indents delineate different levels of the architecture.

makeallstats.m: generates statistics from a set of .dat files and saves them as a single .mat workspace
* loadvelocityff.m: transforms raw voltages to velocities
  - loadLabviewbin.m: loads velocities from Labview-generated .dat files
  - extractprobeinfo.m: receives probe information from a calibration file
  - extractcalibdata.m: receives calibration information and data from a calibration file
  - applytempcorrection.m: corrects voltages for temperatures
  - processcalibrationff.m: finds coefficients for curve fit using calibration file data
    + findKingslaw.m: finds King’s Law coefficients for calibration data
    + calibration file matching format in sample_calibration.m
      ~ velocity calibration file from Dantec StreamWare software
  - applyKingslaw.m: uses King’s Law coefficients to get velocities from voltages

makeallstats_edec.m: calls makeallstats as a function multiple times to make a set of .mat workspaces
* makeallstats.m: same function and dependencies as above

loglogcorf.m: generates important statistics and correlation function plots from .mat workspace created by makeallstats
* hwils.m: calculates the integral length scale of a correlation function

checkIso.m: finds isotropy coefficient (using deprecated calibration method)
* loadvelocityff.m: same function and dependencies as above
* calibparams.m: finds parameters for yaw (angular) calibration

CheckIsoV2.m: finds isotropy coefficient (using improved calibration method, courtesy of Dr. Willem van de Water)
* loadvelocityxwire.m: loadvelocityff.m for x-wire probes
  - loadLabviewbin.m: see above
  - extractprobeinfo.m: see above
  - extractcalibdata.m: see above
  - applytempcorrection.m: see above
  - processcalibrationxwire.m: returns matrix of calibration curve fit parameters (rather than one set)
    + calibration file matching format in sample_calibration_xwire.m
      ~ n velocity calibration files from n StreamWare velocity calibrations at different angles

CorrCheck.m; statistics and plots from angle files (generated by menuII during operation)

EdecfitClient.m: calculates parameters needed to plot energy decay curves from .mat workspaces
* Edecfit.m: plots energy decay function (with or without error bars) for energy decay data sets

cumulativeaverageenergy.m: plot cumulative average over time to see how convergent the mean is by the end of your data set

getrmsangle.m: calculate the RMS angle from a number of angle files in a directory


************************* Alphabetical Summary of Contents *************************
Comprehensive list of the files in this directory, including summaries of functionalities for most

- applyKingslaw.m: uses voltages and King’s Law fit parameters to generates velocities
- applytempcorrection.m: uses data from calibration file to correct voltages for temperatures
- basicstats.m: unknown/unused
- calibparams.m: calculates parameters for yaw (angular) calibration on x-wire probes
- checkIso.m: analyzes x-wire probe data to calculate isotropy coefficient of a data set
- CheckIsoV2.m: analyzes x-wire probe data using improved calibration to calculate isotropy coefficient
- convertDantectemp.m: left over from old temperature conversion attempt, currently unused
- CorrCheck.m: analyzes an angle file and generates plots and statistics for paddle correlations
- CorrCheckAbsoluteValue.m: same as CorrCheck, but defines correlations using absolute values of angles (unused)
- cumulativeaverageenergy.m: plots cumulative average of a data set in time (helps show data convergence)
- Edecfit.m: plots energy decay curves from a set of energy decay tests
- EdecfitClient.m: provides parameters to Edecfit.m so it doesn’t need to be called as a function
- extractangles.m: unknown/unused
- extractcalibdata.m: reads and stores calibration data from a calibration file
- extractexptinfo.m: unknown/unused
- extractfluidinfo.m: reads and stores fluid information from a calibration file (unused)
- extractprobeinfo.m: reads and stores probe information from a calibration file
- extractwiretemps.m: left over from old temperature conversion attempt, currently unused
- findKingslaw.m: finds King’s Law fit parameters from a given calibration data set
- getrmsangle.m: finds actual RMS angles for any number of angle files in a directory
- hwils.m: calculates the integral length scale of a correlation function
- hwspectrum.m: calculates and plots the spectrum of a data set (unused)
- loadLabviewbin.m: reads data from .dat file provided by Labview program
- loadvelocityff.m: reads and converts voltage data into velocity data using calibration curves
- loadvelocityxwire.m: reads and converts voltage data into velocity data using multiangular x-wire calibration
- logfit.m: unknown/unused
- loglogcof.m: calculates flow statistics (e.g. Reynolds number and length scales) and plots correlation function for a given .mat file
- makeallstats.m: calculates statistics for a set of .dat files from Labview
- makeallstats_edec.m: calculates statistics for multiple sets of .dat files (used for energy decay processing)
- makelogtime.m: unknown/unused
- powlawfit.m: unknown/unused
- processcalibrationff.m: processes calibration data from calibration file and returns calibration curve coefficients
- processcalibrationxwire.m: processes calibration data for x-wire probes and returns a matrix of calibration curve coefficients
- sample_calibration.m: sample calibration file, formatted specifically for use with these scripts
- sample_calibration_xwire.m: sample x-wire calibration file with correct formatting
- smoothspectrumnew.m: unknown/unused
- mex files: coded in C/C++ and interface with MATLAB using mex library; used to speed up slow MATLAB algorithms
- xwiretoprobeBruun.m: old calibration processing script for x-wire probes


Questions about architecture or specific programs may be directed to Kevin Griffin (kevinpg@princeton.edu) and Nathan Wei (nwei@princeton.edu).
