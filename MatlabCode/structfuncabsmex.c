
/* 
a mex file for matlab to compute structure functions 
for hot wire data more quickly.  

Sr = structfunc1mex(V, rs, exponent); 

input: 
double V: 1D array, a continuous trace, a component
double rs: 1D array, separations to compute structure functions at
int exponent

output: 
double Sr: 1D array, value for each separation

requires the maximum separation to be the last element in rs!  

Sr = <|dV^exponent|>
dV = V(t+Dt) - V(t)
Dt element of rs

Matlab code: 
cnt = 1; 
for r = rs
	cnt = cnt+1; 
	sr(1:m-r) = V(1+r:m) - V(1:m-r); 
	sr2 = abs(mean(sr(1:n,1:m-r).^2), 2); 
	S2streamtemp(j, k, cnt) = sr2(k); 
end  % for r.rs

Greg Bewley 2010
 */
 
#include "mex.h"
#include "math.h"


void mexFunction(int nlhs, mxArray *plhs[], 
                 int nrhs, const mxArray *prhs[])
{
	  /* inputs: */ 
	double *V; 
	double *rs; 
	int exponent; 
	
	  /* outputs: */
	double *Sr; 
	
	  /* temporary variables: */
	int numseparations, maxseparation, thisseparation, numtimes; 
	double diff, Srsum; 
	int i, j; 
	
	/*mexPrintf("Hello, world!  \n"); */
	
	  /* retrieve the input data */
	V = mxGetPr(prhs[0]); 
	rs = mxGetPr(prhs[1]); 
	exponent = *mxGetPr(prhs[2]); 
	
	  /* Find the dimensions of the data */
	numtimes = mxGetN(prhs[0]); 
	numseparations = mxGetN(prhs[1]); 
	/*maxseparation = max(prhs[2]); */
	maxseparation = (int) rs[numseparations-1]; 
	
	  /* Create an mxArray for the output data */
	plhs[0] = mxCreateDoubleMatrix(numseparations, 1, mxREAL); 
	
	  /* Create a pointer to the output data */
	Sr = mxGetPr(plhs[0]); 
	
	   /* Put data in the output array */
	   /* loop over separations */
	for (i = 0; i < numseparations; i++)
	{
		thisseparation = (int) rs[i]; 
		
		  /* initialize sum */
		Srsum = 0; 
		
		  /* loop over times */
		for (j = 0; j < (numtimes - maxseparation); j++)
		{
			diff = V[j+thisseparation] - V[j]; 
			if (diff < 0) diff = -diff; 
			Srsum = Srsum + pow(diff, exponent); 
		}
		
		  /* turn sum into mean */
		Sr[i] = Srsum / (double) (numtimes - maxseparation); 
	}
	
}
