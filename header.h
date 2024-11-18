#ifndef FUNCTION
#define FUNCTION

#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>

#define eps 1e-15
#define MAX(i, j) (i > j ? y : j)

void f0(double *A, int n, int m, char* filename);

void f1(double *A, int n, int m);

void f2(double *A, int n, int m);
 
void f3(double *A, int n, int m);

void f4(double *A, int n, int m);

void fill(double *A, int n, int m, char *filename, int s);

void fill_right_part(double *A, int n, int m, double *B);





















#endif