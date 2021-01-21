#include <iostream>
#include <math.h>
#include "mpi.h"
#include <unistd.h>

int main(int argc, char *argv[]){
    MPI_Init(&argc, &argv);
    int rank, size;
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    int a[8][8], d[8][8];

    MPI_Datatype newType;
    int array_of_blockLengths[8] = {8, 7, 6, 5, 4, 3, 2, 1};
    int array_of_displacements[8] = {0, 9, 18, 27, 36, 45, 54, 63};
    MPI_Type_indexed(8, array_of_blockLengths, array_of_displacements, MPI_INT, &newType);
    MPI_Type_commit(&newType);

    if (rank == 0){
        printf("Matrix:\n");
        for (int i = 0; i < 8; ++i) {
            for (int j = 0; j < 8; ++j) {
                if (j >= i){
                    a[i][j] = (rand() % 10) + 1;
                } else {
                    a[i][j] = 0;
                }
                printf("%d ", a[i][j]);
            }
            printf("\n");
        }
        MPI_Bcast(&a, 1, newType, 0, MPI_COMM_WORLD);
    } else {
        sleep(rank);
        MPI_Bcast(&d, 1, newType, 0, MPI_COMM_WORLD);
        printf("Thread %d received:\n", rank);
        for (int i = 0; i < 8; ++i) {
            for (int j = 0; j < 8; ++j) {
                if (j < i){
                    printf("0 ");
                } else {
                    printf("%d ", d[i][j]);
                }
            }
            printf("\n");
        }
    }

    MPI_Type_free(&newType);
    MPI_Finalize();
    return 0;
}