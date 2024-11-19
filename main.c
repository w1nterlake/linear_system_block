#include "header.h"

int main(){

    int n, m, r;
    double *A, *B, *X;
    char *text = "text.txt";
    scanf("%d", &n);
    scanf("%d", &m);
    
    A = (double*)malloc(n * n * sizeof(double));
    B = (double*)malloc(n * sizeof(double));
    X = (double*)malloc(n * sizeof(double));

    f1(A, n, m);
    print_matrix(A, n, m, n);
    printf("\n\n");
    solve(A, X, B, n, m);
    printf("\n\n\n");
    print_matrix(A, n, m, n);

    free(A);
    free(B);
    free(X);
    return 0;
}