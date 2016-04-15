#include <stdio.h>

/*----------------------------------
DLL for provisional means algorithm.
Exports provmeans.
To be wrapped in Python with ctypes.
Stand alone VS compile with /MT flag.
------------------------------------*/

void provmeans(double *Xs, double *Ws, int NN, int n, double *sw, double *mn, double *cov){
    double w, r;
    int i,j,k;
    double d[400];
    /* loop over observation vectors */
    for(i=0; i<n; i++){
        w = *Ws;
        *sw += *Ws++;
        r = w/(*sw);
        /* mean */
        for(j=0; j<NN; j++){
            d[j] = Xs[i*NN+j]-mn[j];
            mn[j] += d[j]*r;
        }
       /* covariance */
        for(j=0; j<NN; j++)
            for (k=j; k<NN; k++) cov[j*NN+k] += d[j]*d[k]*(1-r)*w;
    }
}
