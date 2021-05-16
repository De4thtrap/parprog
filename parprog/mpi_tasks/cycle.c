// cyclyed resending

#include <stdio.h>
#include <mpi.h>

void main(int argc, char* argv[])
{
    int rank, commsize;
    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &commsize);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    int sbuffer = rank;
    int rbuffer = 0;
    MPI_Status status;

    if (rank == 0) 
    {
        MPI_Send(&sbuffer, 1, MPI_INT, 1, 0, MPI_COMM_WORLD);
        MPI_Recv(&rbuffer, 1, MPI_INT, commsize - 1, 0, MPI_COMM_WORLD, &status);
    }
    else 
    {
        MPI_Recv(&rbuffer, 1, MPI_INT, (rank - 1 < 0 ? commsize - 1 : rank - 1), 0, MPI_COMM_WORLD, &status);  
        MPI_Send(&sbuffer, 1, MPI_INT, (rank + 1 > commsize ? 0 : rank + 1), 0, MPI_COMM_WORLD);
    }

    printf("Process %d sent value : %d\n", rank, sbuffer);  
    printf("Process %d recieved value : %d\n", rank, rbuffer);

    MPI_Finalize();
}