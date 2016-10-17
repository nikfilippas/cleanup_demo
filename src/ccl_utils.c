#include "ccl_utils.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

/* ------- ROUTINE: ccl_linear spacing ------
INPUTS: [xmin,xmax] of the interval to be divided in bins of width dx
TASK: divide an interval in even bins of width dx; if this doesn't result in even bins, then ?
OUTPUT: bin edges in range [xmin,xmax]
*/

double * ccl_linear_spacing(double xmin, double xmax, double dx, int * N){
    int n = trunc((xmax-xmin)/dx) + 1;
    //is 0.001 this sufficient for our purposes? why not use n as input and determine dx?
    if (fabs(1-(xmax-xmin)/((n-1)*dx))>0.001){ 
        *N = 0;
        fprintf(stderr, "ERROR: Could not evenly divide range [%le, %le] with dx=%le\n", xmin, xmax, dx);
        return NULL;
    }

    double * x = malloc(sizeof(double)*n);
    if (x==NULL){
        fprintf(stderr, "ERROR: Could not allocate memory for linear-spaced array (N=%d)\n", n);
        *N=0;
        return x;
    }
    *N=n;

    for (int i=0; i<n; i++){
        x[i] = xmin + dx*i;
    }

    return x;
}

/* ------- ROUTINE: ccl_log spacing ------
INPUTS: [xmin,xmax] of the interval to be divided logarithmically in N bins
TASK: divide an interval in N logarithmic bins
OUTPUT: bin edges in range [xmin,xmax]
*/

double * ccl_log_spacing(double xmin, double xmax, int N){

    if (N<2){
        fprintf(stderr, "ERROR: Cannot make log-spaced array with %d points - need at least 2\n", N);
        return NULL;
    }

    if (!(xmin>0 && xmax>0)){
        fprintf(stderr, "ERROR: Cannot make log-spaced array xmax or xmax non-positive (had %le, %le)\n", xmin, xmax);
        return NULL;
    }


    double log_xmax = log(xmax);
    double log_xmin = log(xmin);
    double dlog_x = (log_xmax - log_xmin) /  (N-1);


    double * x = malloc(sizeof(double)*N);
    if (x==NULL){
        fprintf(stderr, "ERROR: Could not allocate memory for log-spaced array (N=%d)\n", N);
        return x;
    }


    for (int i=0; i<N; i++){
        x[i] = exp(log_xmin + dlog_x*i);
    }

    return x;
}


