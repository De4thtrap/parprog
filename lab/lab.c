#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

#define DEFAULT_N 1000000

#define INITIAL_TAG       1
#define READY_TAG         2
#define TRANSMISSION_ADD_TAG  1000000000


struct param_type {
    int start_m; // included
    int end_m;   // excluded
    int K;
};

double f(double x, double t) {
    return t * x;
}

double fi(double x) {
    return x * x - x * x * x / 6;
}

double psi(double t) {
    return t * t;
}

double exact_solution(double x, double t) {
    return (x - t) * (x - t) + 0.5 * t * x * x - x * x * x / 6;
}

double get_cell(int m, int k, double** U, double h, double tau, int start_m) {
    double res = (-U[m-1][k] + U[m-1][k-1] + U[m][k-1]) / 2 / tau;
    res += (U[m-1][k] + U[m-1][k-1] - U[m][k-1]) / 2 / h;
    res += f(h * (m+start_m-1+0.5), tau * (k+0.5));
    res /= 0.5 * (1/h + 1/tau);
    return res;
}

void make_line(int k, int M, double** U, double h, double tau, int start_m) {
    int m;
    for (m = 1; m < M; m++) {
        U[m][k] = get_cell(m, k, U, h, tau, start_m);
    }
}


void slave_task(int rank, int size, double h, double tau) {
    struct param_type param;
    MPI_Recv(&param, sizeof(struct param_type), MPI_INT, 0, INITIAL_TAG, MPI_COMM_WORLD, NULL);

    int K = param.K;
    int end_m   = param.end_m;
    int start_m = param.start_m;
    int M = end_m - start_m + 1;

    double **U = (double**)malloc(M * sizeof(double*));
    int m;
    for (m = 0; m < M; ++m) {
        U[m] = (double*)malloc(K * sizeof(double));
    }

    for (m = 0; m < M; m++) {
        U[m][0] = fi((m + start_m - 1) * h);
    }

    int k;
    for (k = 1; k < K; k++) {
        if (rank == 1) {
            U[0][k] = psi(k * tau);
        } else {
            MPI_Recv(&U[0][k], 1, MPI_DOUBLE, rank - 1, READY_TAG, MPI_COMM_WORLD, NULL);
        }

        make_line(k, M, U, h, tau, start_m);

        if (rank != size - 1) {
            MPI_Send(&U[M-1][k], 1, MPI_DOUBLE, rank + 1, READY_TAG, MPI_COMM_WORLD);
        }
    }
    
    if (rank != 1) {
        free(U[0]);
    }

    for (m = (rank == 1) ? 0 : 1; m < M; ++m) {
        MPI_Send(U[m], K, MPI_DOUBLE, 0, m + start_m - 1 + TRANSMISSION_ADD_TAG, MPI_COMM_WORLD);
        // printf("rank = %2d   m = %7d sent\n", rank, m);
        free(U[m]);
    }

    free(U);
}


double** master_task(int K, int M) {
    int size;
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    int start_m = 1;
    int end_m   = 1;
    int len = M / (size - 1);
    int mod = M % (size - 1);

    // gives tasks to all other nodes
    int i;
    for (i = 1; i < size; i++) {
        start_m = end_m;                       // start x index (included)
        end_m += len + (i <= mod ? 1 : 0);     //   end x index (excluded)

        struct param_type param;
        param.start_m = start_m;
        param.end_m   = end_m;
        param.K = K;

        MPI_Send(&param, 3, MPI_INT, i, INITIAL_TAG, MPI_COMM_WORLD);
    }


    // U(x, t)
    double **U = (double**)malloc((M+1) * sizeof(double*));
    int m;
    for (m = 0; m <= M; ++m) {
        U[m] = (double*)malloc(K * sizeof(double));
    }

    for (m = 0; m <= M; m++){
        MPI_Recv(U[m], K, MPI_DOUBLE, MPI_ANY_SOURCE, m + TRANSMISSION_ADD_TAG, MPI_COMM_WORLD, NULL);
        // printf("m = %7d recieved\n", m);
    }

    return U;
}


double** single_node_version(int K, int M, double h, double tau) {
    double **U = (double**)malloc(M * sizeof(double*));
    int m;
    for (m = 0; m < M; ++m) {
        U[m] = (double*)malloc(K * sizeof(double));
    }

    for (m = 0; m < M; m++) {
        U[m][0] = fi(m * h);
    }

    int k;
    for (k = 1; k < K; k++) {
        U[0][k] = psi(k * tau);
        make_line(k, M, U, h, tau, 1);
    }

    return U;
}






int main(int argc, char** argv) {
    MPI_Init(&argc, &argv);
    

    int K = strtol(argv[1], NULL, 10); // max t number
    int M = strtol(argv[2], NULL, 10); // max x number
    double tau = strtod(argv[3], NULL);
    double h   = strtod(argv[4], NULL);
    char if_print = (argc >= 6 && !strcmp(argv[5], "print"));
    

    int rank;
    int size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    // U(x, t)
    double** U;

    if (size == 1) {
        double start = MPI_Wtime();
        U = single_node_version(K + 1, M + 1, h, tau);
        double duration = MPI_Wtime() - start;
        printf("duration is %lf\n", duration);
    } else if (rank == 0) { // master branch
        double start = MPI_Wtime();
        U = master_task(K + 1, M);
        double duration = MPI_Wtime() - start;
        printf("process %2d (master) duration is %lf (same as full program duration)\n", rank, duration);
    } else { // slave branch
        double start = MPI_Wtime();
        slave_task(rank, size, h, tau);
        double duration = MPI_Wtime() - start;
        printf("process %2d duration is %lf\n", rank, duration);
    }
    
    
    // printing results and clearing memory
    if (rank == 0) {
        if (if_print) {
            printf("\nnumerical solution:\n");
            int k, m;
            for (k = K; k >= 0; --k) {
                for (m = 0; m <= M; ++m)
                    printf("%10lf, ", U[m][k]);
                printf("\n");
            }
            printf("\nexact solution:\n");
            for (k = K; k >= 0; --k) {
                for (m = 0; m <= M; ++m)
                    printf("%10lf, ", exact_solution(m * h, k * tau));
                printf("\n");
            }
        }

        int m;
        for (m = 0; m <= M; ++m)
            free(U[m]);
        free(U);
    }


    MPI_Finalize();
    return 0;
}
