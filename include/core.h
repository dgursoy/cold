#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>


#define MAX(x, y) (((x) >= (y)) ? (x) : (y))
#define MIN(x, y) (((x) < (y)) ? (x) : (y))
#define PI 3.14159265358979323846


int recposlsqr(float* data, int dx, float* mask, int mx, int corr, int pos, int sig, float tau, float dec);
int recsiglsqr(float* data, int dx, float* mask, int mx, int corr, int pos, int sig, float tau, float dec);
float rectaulsqr(float* data, int dx, float* mask, int mx, int corr, int pos, int sig, float tau, float dec);
float recdeclsqr(float* data, int dx, float* mask, int mx, int corr, int pos, int sig, float tau, float dec);
void convolve(float* signal, int ns, float* kernel, int nk, float* res);
void footprint(int size, float tau, float dec, float* signal);
void copyarr(float* arr, int size, float* copy);
void maxarr(float* arr, int size, float *maxval);
void minarr(float* arr, int size, float *minval);
void minmaxarr(float* arr, int size, float *minval, float *maxval);
void normalize(float* arr, int size);
