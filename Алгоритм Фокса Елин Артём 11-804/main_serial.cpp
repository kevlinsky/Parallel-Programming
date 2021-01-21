#include <math.h>
#include <cstdio>

int main(){
    int N = 10;
    int A[N][N];
    int B[N][N];
    int C[N][N];

    printf("Matrix A:\n");
    for(int i = 0; i < N; i++){
        for(int j = 0; j < N; j++){
            A[i][j] = rand() % 10;
            printf("%d ", A[i][j]);
        }
        printf("\n");
    }

    printf("Matrix B:\n");
    for(int i = 0; i < N; i++){
        for(int j = 0; j < N; j++){
            B[i][j] = rand() % 10;
            printf("%d ", B[i][j]);
        }
        printf("\n");
    }

    printf("Result matrix C:\n");
    for(int i = 0; i < N; i++){
        for(int j = 0; j < N; j++){
            C[i][j] = 0;
            for (int k = 0; k < N; k++) {
                C[i][j] += A[i][k] * B[k][j];
            }
            printf("%d ", C[i][j]);
        }
        printf("\n");
    }

    return 0;
}
