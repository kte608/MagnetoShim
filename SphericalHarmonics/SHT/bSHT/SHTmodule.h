#ifndef KDEF_SHT_MODULE
#define KDEF_SHT_MODULE
#include <Python.h>
#include <errno.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include "makeweights.h"
#include "fftw3.h"
//#include "makeweights.h"
#include "FST_semi_fly.h"
//#include "csecond.h"

#define max(A, B) ((A) > (B) ? (A) : (B))
double* listtoarray(PyObject* o,int *size);
void SHT_BACKEND(long int bw,double *rdata,double *idata,
		 double **rresult,double **iresult);
double harmcoeff(int n,int m);
#endif
