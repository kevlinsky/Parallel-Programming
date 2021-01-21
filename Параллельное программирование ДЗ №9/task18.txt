#include <iostream>
#include <math.h>
#include "mpi.h"
#include <unistd.h>

int main(int argc, char *argv[]){
    MPI_Init(&argc, &argv);
    int rank, size;
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    int a[8][8], b[4][8], c[4][8];

    MPI_Datatype newType;
    MPI_Type_vector(4, 8, 16, MPI_INT, &newType);
    MPI_Type_commit(&newType);

    if (rank == 0){
        printf("Matrix:\n");
        for (int i = 0; i < 8; ++i) {
            for (int j = 0; j < 8; ++j) {
                a[i][j] = (rand() % 1000);
                printf("%d ", a[i][j]);
            }
            printf("\n");
        }
        MPI_Send(&a[0][0], 1, newType, 1, 1, MPI_COMM_WORLD);
        MPI_Send(&a[1][0], 1, newType, 2, 2, MPI_COMM_WORLD);
    } else if (rank == 1){
        MPI_Status status;
        int size;

        MPI_Probe(0, 1, MPI_COMM_WORLD, &status);
        MPI_Get_count(&status, MPI_INT, &size);

        int *tmp = new int[size];
        MPI_Recv(tmp, size, MPI_INT, 0, 1, MPI_COMM_WORLD, &status);

        int count = 0;
        for (int i = 0; i < 4; i++) {
            for (int j = 0; j < 8; j++) {
                b[i][j] = tmp[count];
                count++;
            }
        }
        printf("Thread 1 received:\n");
        for (int i = 0; i < 4; ++i) {
            for (int j = 0; j < 8; ++j) {
                printf("%d ", b[i][j]);
            }
            printf("\n");
        }
    } else {
        sleep(1);
        MPI_Status status;
        int size;

        MPI_Probe(0, 2, MPI_COMM_WORLD, &status);
        MPI_Get_count(&status, MPI_INT, &size);

        int *tmp = new int[size];
        MPI_Recv(tmp, size, MPI_INT, 0, 2, MPI_COMM_WORLD, &status);

        int count = 0;
        for (int i = 0; i < 4; i++) {
            for (int j = 0; j < 8; j++) {
                c[i][j] = tmp[count];
                count++;
            }
        }
        printf("Thread 2 received:\n");
        for (int i = 0; i < 4; ++i) {
            for (int j = 0; j < 8; ++j) {
                printf("%d ", c[i][j]);
            }
            printf("\n");
        }
    }
    MPI_Type_free(&newType);
    MPI_Finalize();
    return 0;
}