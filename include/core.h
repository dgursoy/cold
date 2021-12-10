#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>


#define MAX(x, y) (((x) >= (y)) ? (x) : (y))
#define MIN(x, y) (((x) < (y)) ? (x) : (y))
#define PI 3.14159265358979323846


int recposlsqr(float* data, int dx, float* mask, int mx, int corr, int pos, int sig, float tau);
int recsiglsqr(float* data, int dx, float* mask, int mx, int corr, int pos, int sig, float tau);
float rectaulsqr(float* data, int dx, float* mask, int mx, int corr, int pos, int sig, float tau);
void convolve(float* signal, int ns, float* kernel, int nk, float* res);
void gauss(int siz, float mu, float sigma, float* signal);
void lorentz(int siz, float mu, float sigma, float* signal);
void pseudovoigt(int siz, float mu, float sigma, float alpha, float* signal);
void tukey(int size, float tau, float* signal);
void copyarr(float* arr, int size, float* copy);
void maxarr(float* arr, int size, float *maxval);
void minarr(float* arr, int size, float *minval);
void minmaxarr(float* arr, int size, float *minval, float *maxval);
void normalize(float* arr, int size);
