#include "header.h"

bool f0(double *A, int n, int m, char *filename){
    FILE *fin;
    int vert, hor;
    int k = n/m;
    int l = n - k * m;
    fin = fopen(filename, "r");
    int count = 0;
    for(int i = 0; i * m < n; i++){
        vert = i < k ? m : l;
        for(int i1 = 0; i1 < vert; i1++){
            for(int j = 0; j * m < n; j++){
                hor = j < k ? m : l;
                for(int j1 = 0; j1 < hor; j1++){
                    if(fscanf(fin, "%lf", &A[i * n * m + j * vert * m + i1 * hor + j1]) == 1){
                        count++;
                    }
                }
            }
        }
    }
    fclose(fin);
    if(count == n * n){
        return true;
    }
    else{
        return false;
    }
}

bool f1(double *A, int n, int m){
    int vert, hor;
    int k = n/m;
    int l = n - k * m;
    for(int i = 0; i * m < n; i++){
        vert = i < k ? m : l;
        for(int i1 = 0; i1 < vert; i1++){
            for(int j = 0; j * m < n; j++){
                hor = j < k ? m : l;
                for(int j1 = 0; j1 < hor; j1++){
                    A[i * m * n + j * vert * m + i1 * hor + j1] = n - MAX(i * m + i1 + 1, j * m + j1 + 1) + 1;
                }
            }
        }
    }
    return true;
}

bool f2(double *A, int n, int m){
    int vert, hor;
    int k = n/m;
    int l = n - k * m;
    for(int i = 0; i * m < n; i++){
        vert = i < k ? m : l;
        for(int i1 = 0; i1 < vert; i1++){
            for(int j = 0; j * m < n; j++){
                hor = j < k ? m : l;
                for(int j1 = 0; j1 < hor; j1++){
                    A[i * m * n + j * vert * m + i1 * hor + j1] = MAX(i * m + i1 + 1, j * m + j1 + 1);
                }
            }
        }
    }
    return true;
}

bool f3(double *A, int n, int m){
    int vert, hor;
    int k = n/m;
    int l = n - k * m;
    for(int i = 0; i * m < n; i++){
        vert = i < k ? m : l;
        for(int i1 = 0; i1 < vert; i1++){
            for(int j = 0; j * m < n; j++){
                hor = j < k ? m : l;
                for(int j1 = 0; j1 < hor; j1++){
                    A[i * m * n + j * vert * m + i1 * hor + j1] = abs((i * m + i1) - (j * m + j1));
                }
            }
        }
    }
    return true;
}

bool f4(double *A, int n, int m){
    int vert, hor;
    int k = n/m;
    int l = n - k * m;
    for(int i = 0; i * m < n; i++){
        vert = i < k ? m : l;
        for(int i1 = 0; i1 < vert; i1++){
            for(int j = 0; j * m < n; j++){
                hor = j < k ? m : l;
                for(int j1 = 0; j1 < hor; j1++){
                    A[i * m * n + j * vert * m + i1 * hor + j1] = 1 / (double)(i * m + i1 + j * m + j1);
                }
            }
        }
    }
    return true;    
}

void fill(double *A, int n, int m, char *filename, int s, int q){
    switch(s){
        case(0):
            if(!f0(A, n, m, filename)){
                printf("Ошибка при записи с файла");
                q = 1;
                return;
            }
            else{
                return;
            }

        case(1):
            f1(A, n, m);
            return;
        case(2):
            f2(A, n, m);
            return;
        case(3):
            f3(A, n, m);
            return;
        case(4):
            f4(A, n, m);
            return;
        default:
            return;
    }
}

void fill_right_part(double *A, int n, int m, double *B){
    int vert, hor;
    int k = n / m;
    int l = n - k*m;
    for(int i = 0; i * m < n; i++){
        vert = i < k ? m : l;
        for(int y = 0; y < vert; y++){
            B[i * m + y] = 0;
            for(int j = 0; j * m < n; j++){
                hor = j < k ? m : l;
                for(int x = 0; x < hor; x += 2){
                    B[i * m + y] += A[i * m * n + j * vert * m + y * hor + x];
                }
            }
        }
    }
}

// void multiply_vec_on_matrix(double *A, double *B, double *X, int n, int m){
//     int vert, hor, tmp;
//     int k = n / m;
//     int l = n - k * m;
//     double sum = 0;
//     for(int i = 0; i * m < n; i++){
//         vert = i < k ? m : l;
//         for(int y = 0; y < vert; y++){
//             for(int j = 0; j * m < n; j++){
//                 hor = j < k ? m : l;
//                 for(int x = 0; x * m < n; x++){
//                     sum += A[i * m * n + j * vert * m + y * hor + x] * B[i * m + y];
//                 }
//             }
//             X[i * m + y] = sum;
//             sum = 0;
//         }
//     }
// }

double discrepancy_1(double *A, double *B, double *X, int n, int m){
    int vert, hor;
    int k = n / m;
    int l = n - k*m;
    double sum = 0;
    double res1 = 0;
    double res2 = 0;
    for(int i = 0; i * m < n; i++){
        vert = i < k ? m : l;
        for(int y = 0; y < vert; y++){
            for(int j = 0; j * m < n; j++){
                hor = j < k ? m : l;
                for(int x = 0; x < hor; x++){
                    sum += A[i * m * n + j * vert * m + y * hor + x] * X[i * m + y];
                }
            }
            res1 += fabs(sum - B[i * m + y]);
            res2 += fabs(B[i * m + y]);
            sum = 0;
        }
    }
    if(res2 < eps){
        return -1;
    }
    else{
        return res1 / res2;
    }
}

double discrepancy_2(double *X, int n, int m){
    int z;
    int k = n / m;
    int l = n - k * m;
    double sum1 = 0;
    double sum2 = 0;
    for(int i = 0; i * m < n; i++){
        z = i < k ? m : l;
        for(int j = 0; j < z; j++){
            sum1 += fabs(X[i * m + j] - (i * m + j + 1) % 2);
            sum2 += (i * m + j + 1) % 2;
        }
    }
    return sum1 / sum2;
}

void made_vec(double *A, double *X, int n, int m, int i, int j, int q, int p){
    if(i == j){
        int k = n / m;
        int l = n - k * m;
        int vert = i < k ? m : l;
        int hor = i < k ? m : l;
        for(int y = 0; y < vert; y++){
            X[y + i * m] = 0;
        }
        double sum = 0;
        for(int y = q; y < vert; y++){
            X[y + i * m] = A[i * n * m + j * vert * m + y * hor + p];
            sum += X[y + i * m];
        }
        sum = sqrt(sum);
        X[q + i * m] /= sum;
        sum = 0;
        for(int y = q; y < vert; y++){
            sum += X[y + i * m];
        }
        sum = sqrt(sum);
        for(int y = q; y < vert; y++){
            X[y + i * m] /= sum;
        }
    }//
    else{
        int k = n / m;
        int l = n - k * m;
        int vert = i < k ? m : l;
        int hor = j < k ? m : l;
        for(int y = 0; y < vert; y++){
            X[y + i * m] = 0;
        }
        double sum = 0;
        X[p + j * m] = A[j * n * m + j * m * vert + hor * p + p];
        sum += X[j * m + p];
        for(int y = 0; y < vert; y ++){
            X[i * m + y] = A[i * m * n + j * m * vert + y * hor + p];
            sum += X[y + i * m];
        }
        sum = sqrt(sum);
        X[j * m + p] -= sum;
        sum = X[j * m + p];
        for(int y = 0; y < vert; y++){
            sum += X[y + i * m];
        }
        sum = sqrt(sum);
        X[j * m + p] /= sum;
        for(int y = 0; y < vert; y++){
            X[i * m + y] /= sum;
        }
    }
}

void multiply(double *A, double *X, int n, int m, int i, int j, int p){
    if(i == j){
        int k = m / n;
        int l = n - k * m;
        int vert = i < k ? m : l;
        int hor;
        for(int z = j; z * m < n; z++){
            hor = j < k ? m : l;
            for(int x = 0; x < hor; x++){
                double sum = 0;
                for(int y = 0; y < vert; y++){
                    sum += X[i * m + y] * A[i * m * n + z * m * vert + y * hor + x];
                }
                for(int y = 0; y < vert; y++){
                    A[i * m * n + z * m * vert + y * hor + x] -= 2 * sum * X[i * m + y];
                }
                sum = 0;
            }
        }
    }
    else{
        int k = m / n;
        int l = n - k * m;
        int vert = i < k ? m : l;
        int hor;
        for(int z = j; z * m < n; z++){
            hor = j < k ? m : l;
            for(int x = 0; x < hor; x++){
                double sum = 0;
                sum += A[j * m * n + j * m * vert + p * hor + x] * X[j * m + p];
                for(int y = 0; y < vert; y++){
                    sum += A[i * m * n + z * m * vert + y * vert + x] * X[i * m + y];
                }
                A[j * m * n + j * m * vert + p * hor + x] -= 2 * sum * X[j * m + p];
                for(int y = 0; y < vert; y++){
                    A[i * m * n + z * m * vert + y * hor + x] -= 2 * sum * X[i * m + y];
                }
                sum = 0;
            }
        }
    }
}

double solve(double *A, double *B, double *X, int n, int m){

}


