#include <iostream>
#include <math.h>
#include "mpi.h"
#include <unistd.h>

struct MyStruct{
    int d[2][8];
};

int main(int argc, char *argv[]){
    MPI_Init(&argc, &argv);
    int rank, size;
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    MyStruct myStruct;
    int blockLength[1] = {16};
    MPI_Datatype type[1] = {MPI_INT};
    MPI_Aint address[2];
    MPI_Aint disp[1];
    MPI_Datatype newType;
    int a[8][8];

    MPI_Get_address(&myStruct, &address[0]);
    MPI_Get_address(&(myStruct.d[0][0]), &address[1]);
    disp[0] = address[1] - address[0];

    MPI_Type_create_struct(1, blockLength, disp, type, &newType);
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
        int index = 0;
        int prev_index = 0;
        for (int k = 0; k < 4; ++k) {
            for (int i = 0; i < 2; ++i) {
                for (int j = 0; j < 8; ++j) {
                    myStruct.d[i][j] = a[index][j];
                }
                index += 4;
            }
            index = prev_index + 1;
            prev_index++;
            MPI_Send(&myStruct, 1, newType, k + 1, k + 1, MPI_COMM_WORLD);
        }
    } else {
        sleep(rank);
        MPI_Status status;
        MPI_Recv(&myStruct, 1, newType, 0, rank, MPI_COMM_WORLD, &status);

        printf("Thread %d received:\n", rank);
        for (int i = 0; i < 2; ++i) {
            for (int j = 0; j < 8; ++j) {
                printf("%d ", myStruct.d[i][j]);
            }
            printf("\n");
        }
    }
    MPI_Type_free(&newType);
    MPI_Finalize();
    return 0;
}