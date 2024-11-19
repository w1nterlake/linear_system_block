#ifndef FUNCTION
#define FUNCTION

#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include <stdbool.h>

#define eps 1e-15
#define MAX(i, j) (i > j ? i : j)

bool f0(double *A, int n, int m, char* filename);//

bool f1(double *A, int n, int m);//

bool f2(double *A, int n, int m);//
 
bool f3(double *A, int n, int m);//

bool f4(double *A, int n, int m); //

void fill(double *A, int n, int m, char *filename, int s, int *q);//

void fill_right_part(double *A, int n, int m, double *B);//

// void multiply_vec_on_matrix(double *A, double *B, double *X, int n, int m);

double discrepancy_1(double *A, double *B, double *X, int n, int m);

double discrepancy_2(double *X, int n, int m);

double solve(double *A, double *X, double *B, int n, int m);

void made_vec(double *A, double *X, int n, int m, int i, int j, int q, int p, double *S);

void multiply(double *A, double *X, int n, int m, int i, int j, int p, double *S);

void multiply_right_part(double *X, double *B, int n, int m, int i, int j, int z);

void print_matrix(double *A, int n, int m, int r);







#endif