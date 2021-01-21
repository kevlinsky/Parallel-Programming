#include <iostream>
#include <mpi.h>
#include <cmath>

int main(int argc, char *argv[]) {
    srand(time(nullptr));
    int size = 10;

    int newGroupOfProcesses[5] = { 8, 3, 9, 1, 6 };

    int world_rank;
    int world_size;

    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);
    MPI_Comm_size(MPI_COMM_WORLD, &world_size);

    MPI_Group worldGroup;
    MPI_Group newGroup;
    MPI_Comm_group(MPI_COMM_WORLD, &worldGroup);
    MPI_Group_incl(worldGroup, 5, newGroupOfProcesses, &newGroup);

    MPI_Comm newComm;
    MPI_Comm_create(MPI_COMM_WORLD, newGroup, &newComm);

    int a[size];

    int newRank = -1;
    if (newComm != MPI_COMM_NULL) {
        MPI_Comm_rank(newComm, &newRank);
    }

    if (newRank == 0 && newComm != MPI_COMM_NULL) {
        for (int i = 0; i < size; i++) {
            a[i] = rand() % 10 + 1;
        }
        MPI_Bcast(&a, size, MPI_INT, 0, newComm);
    }

    if (newRank != -1 && newRank != 0 && newComm != MPI_COMM_NULL) {
        int a1[size];
        MPI_Bcast(&a1, size, MPI_INT, 0, newComm);

        printf("Thread %d from NEW_COMM received a[i] = ", newRank);
        for (int i = 0; i < size; i++) {
            printf("%d ", a1[i]);
        }
        printf("\n");

        if (newRank == 4) {
            MPI_Send(&a1, size, MPI_INT, 0, 0, MPI_COMM_WORLD);
        }
    }

    if (world_rank == 0) {
        int a2[size];
        MPI_Status status;
        MPI_Recv(a2, size, MPI_INT, 6, 0, MPI_COMM_WORLD, &status);

        printf("Thread 0 from WORLD_COMM received from thread 4 from NEW_COMM a[i] = ");
        for (int i = 0; i < size; i++) {
            printf("%d ", a2[i]);
        }
        printf("\n");
    }

    if (newComm != MPI_COMM_NULL) {
        MPI_Comm_free(&newComm);
    }

    MPI_Group_free(&newGroup);
    MPI_Finalize();
    return 0;
}