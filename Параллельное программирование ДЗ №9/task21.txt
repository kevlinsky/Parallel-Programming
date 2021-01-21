#include <iostream>
#include <math.h>
#include "mpi.h"
#include <unistd.h>

struct MyNewStruct{
    int i[10];
    float f[8];
};

int main(int argc, char *argv[]){
    MPI_Init(&argc, &argv);
    int rank, size;
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    MyNewStruct myNewStruct;
    int blockLengths[2] = {10, 8};
    MPI_Datatype type[2] = {MPI_INT, MPI_FLOAT};
    MPI_Aint address[3];
    MPI_Aint disp[2];
    MPI_Datatype newType;

    MPI_Get_address(&myNewStruct, &address[0]);
    MPI_Get_address(&(myNewStruct.i[0]), &address[1]);
    MPI_Get_address(&(myNewStruct.f[0]), &address[2]);
    disp[0] = address[1] - address[0];
    disp[1] = address[2] - address[0];

    MPI_Type_create_struct(1, blockLengths, disp, type, &newType);
    MPI_Type_commit(&newType);

    if (rank == 0){
        printf("Initial structure:\ni[10] = ");
        for (int i = 0; i < 10; ++i) {
            myNewStruct.i[i] = rand() % 10 + 1;
            printf("%d ", myNewStruct.i[i]);
        }
        printf("\nf[8] = ");
        for (int i = 0; i < 8; ++i) {
            myNewStruct.f[i] = ((float) (rand() % 10 + 1)) / 10;
            printf("%f ", myNewStruct.f[i]);
        }
        printf("\n");
        MPI_Bcast(&myNewStruct, 1, newType, 0, MPI_COMM_WORLD);
    } else {
        sleep(rank);
        MyNewStruct myNewStruct1;
        MPI_Bcast(&myNewStruct1, 1, newType, 0, MPI_COMM_WORLD);
        printf("Thread %d received:\ni[10] = ", rank);
        for (int i = 0; i < 10; ++i) {
            printf("%d ", myNewStruct1.i[i]);
        }
        printf("\nf[8] = ");
        for (int i = 0; i < 8; ++i) {
            printf("%f ", myNewStruct1.f[i]);
        }
        printf("\n");
    }

    MPI_Type_free(&newType);
    MPI_Finalize();
    return 0;
}

