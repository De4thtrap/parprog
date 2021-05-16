// harmonic sum

#include <stdio.h>
#include <mpi.h>

int main(int argc, char* argv[])
{
    int commsize, rank, root = 0;
    double result = 0;
    double partial = 0;

    const int N = 100;    

    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &commsize);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

//    MPI_Bcast(&source, 1, MPI_DOUBLE, root, MPI_COMM_WORLD);

    if (rank != root) 
    {
        for (int n = (rank - 1) * N/commsize; n < rank * N/commsize; n++)
            partial += n > 0 ? 1/n : 0;
    }

    printf("MPI_SUM result on process #%d is %d\n", rank, result);

    MPI_Reduce(&partial, &result, 1, MPI_DOUBLE, MPI_SUM, root, MPI_COMM_WORLD);

    if (rank == root) 
        printf("MPI_SUM result is %d\n", result);

    MPI_Finalize();
    
    return 0;
}
