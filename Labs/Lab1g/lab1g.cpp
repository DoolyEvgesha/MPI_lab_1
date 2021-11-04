#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <mpi.h>
#include <algorithm>

#define ISIZE 10
#define JSIZE 10
#define MASTER_RANK 0
#define IN 3
#define JN 3

#define TAG_PLACE 1
#define TAG_ARRAY_ELEMENT 2
#define TAG_STOP 3
#define TAG_ARRAY_NUM 4

void printInFile(FILE *ff, double **a) {
    for (int i = 0; i < ISIZE; i++) {
        for (int j = 0; j < JSIZE; j++) {
            fprintf(ff, "%f ", a[i][j]);
        }
        fprintf(ff, "\n");
    }
}

double *create1DArray(int size) {
    auto *a = static_cast<double *>(calloc(size, sizeof(double)));
    return a;
}

double **create2DArray() {
    auto **a = static_cast<double **>(calloc(JSIZE, sizeof(double *)));
    for (int i = 0; i < JSIZE; i++) {
        a[i] = create1DArray(ISIZE);
    }
    return a;
}

void free1DArray(double *a) {
    free(a);
}

void free2DArray(double **a) {
    for (int i = 0; i < ISIZE; i++) {
        free1DArray(a[i]);
    }
    free(a);
}

void computeCycleMaster(double **a, int numThreads) {
    MPI_Status status;
    int rowNum, place;
    for (int i = 0; i < numThreads - 1; i++) {
        // send to all slave threads
        MPI_Send(&i, 1, MPI_INT, i + 1, TAG_PLACE, MPI_COMM_WORLD);
        MPI_Send(a[i], JSIZE, MPI_DOUBLE, i + 1, TAG_ARRAY_ELEMENT, MPI_COMM_WORLD);
    }

    for (int i = numThreads - 1; i < ISIZE; i++) {
        MPI_Recv(&rowNum, 1, MPI_INT, MPI_ANY_SOURCE, TAG_PLACE, MPI_COMM_WORLD, &status);
        MPI_Recv(a[rowNum], JSIZE, MPI_DOUBLE, status.MPI_SOURCE, TAG_ARRAY_ELEMENT, MPI_COMM_WORLD, &status);

        MPI_Send(&i, 1, MPI_INT, status.MPI_SOURCE, TAG_PLACE, MPI_COMM_WORLD);
        MPI_Send(a[i], JSIZE, MPI_DOUBLE, status.MPI_SOURCE, TAG_ARRAY_ELEMENT, MPI_COMM_WORLD);
    }

    for (int i = 0; i < numThreads - 1; i++) {
        MPI_Recv(&rowNum, 1, MPI_INT, MPI_ANY_SOURCE, TAG_PLACE, MPI_COMM_WORLD, &status);
        MPI_Recv(a[rowNum], JSIZE, MPI_DOUBLE, status.MPI_SOURCE, TAG_ARRAY_ELEMENT, MPI_COMM_WORLD, &status);

        MPI_Send(&rowNum, 1, MPI_INT, status.MPI_SOURCE, TAG_STOP, MPI_COMM_WORLD);
    }

    printf("[MASTER] finished 3rd cycle\n");
    fflush(stdout);

    double *a1 = create1DArray(JSIZE);
    double *a2 = create1DArray(JSIZE);
    double *a3 = create1DArray(JSIZE);

    int jMax = std::min(numThreads - 1, 3);

    for (int i = 1; i < ISIZE; i++) {
        for (int j = 0; j < jMax; j++) {
            MPI_Send(&j, 1, MPI_INT, j + 1, TAG_ARRAY_NUM, MPI_COMM_WORLD);
            MPI_Send(&i, 1, MPI_INT, j + 1, TAG_PLACE, MPI_COMM_WORLD);
            MPI_Send(a[i], JSIZE, MPI_DOUBLE, j + 1, TAG_ARRAY_ELEMENT, MPI_COMM_WORLD);
            MPI_Send(a[i - 1], JSIZE, MPI_DOUBLE, j + 1, TAG_ARRAY_ELEMENT, MPI_COMM_WORLD);
        }
        for (int j = jMax; j < 3; j++) {
            MPI_Recv(&place, 1, MPI_INT, MPI_ANY_SOURCE, TAG_ARRAY_NUM, MPI_COMM_WORLD, &status);
            MPI_Recv(&rowNum, 1, MPI_INT, status.MPI_SOURCE, TAG_PLACE, MPI_COMM_WORLD, &status);
            if (place == 0) {
                MPI_Recv(a1, JSIZE, MPI_DOUBLE, status.MPI_SOURCE, TAG_ARRAY_ELEMENT, MPI_COMM_WORLD, &status);
            } else if (place == 1) {
                MPI_Recv(a2, JSIZE, MPI_DOUBLE, status.MPI_SOURCE, TAG_ARRAY_ELEMENT, MPI_COMM_WORLD, &status);
            } else if (place == 2) {
                MPI_Recv(a3, JSIZE, MPI_DOUBLE, status.MPI_SOURCE, TAG_ARRAY_ELEMENT, MPI_COMM_WORLD, &status);
            }
            MPI_Send(&j, 1, MPI_INT, j + 1, TAG_ARRAY_NUM, MPI_COMM_WORLD);
            MPI_Send(&i, 1, MPI_INT, status.MPI_SOURCE, TAG_PLACE, MPI_COMM_WORLD);
            MPI_Send(a[i], JSIZE, MPI_DOUBLE, status.MPI_SOURCE, TAG_ARRAY_ELEMENT, MPI_COMM_WORLD);
            MPI_Send(a[i - 1], JSIZE, MPI_DOUBLE, status.MPI_SOURCE, TAG_ARRAY_ELEMENT, MPI_COMM_WORLD);
        }
        for (int j = 0; j < jMax; j++) {
            MPI_Recv(&place, 1, MPI_INT, MPI_ANY_SOURCE, TAG_ARRAY_NUM, MPI_COMM_WORLD, &status);
            MPI_Recv(&rowNum, 1, MPI_INT, status.MPI_SOURCE, TAG_PLACE, MPI_COMM_WORLD, &status);
            if (place == 0) {
                MPI_Recv(a1, JSIZE, MPI_DOUBLE, status.MPI_SOURCE, TAG_ARRAY_ELEMENT, MPI_COMM_WORLD, &status);
            } else if (place == 1) {
                MPI_Recv(a2, JSIZE, MPI_DOUBLE, status.MPI_SOURCE, TAG_ARRAY_ELEMENT, MPI_COMM_WORLD, &status);
            } else if (place == 2) {
                MPI_Recv(a3, JSIZE, MPI_DOUBLE, status.MPI_SOURCE, TAG_ARRAY_ELEMENT, MPI_COMM_WORLD, &status);
            }

            MPI_Send(&rowNum, 1, MPI_INT, status.MPI_SOURCE, TAG_STOP, MPI_COMM_WORLD);
        }
        for (int j = 3; j < JSIZE; j++) {
            if (j % 3 == 0) {
                a[i][j] = a1[j];
            } else if (j % 3 == 1) {
                a[i][j] = a2[j];
            } else if (j % 3 == 2) {
                a[i][j] = a3[j];
            }
        }
    }

    printf("[MASTER] finished 4th cycle\n");
    fflush(stdout);

//    for (int i = 0; i < numThreads - 1; i++) {
//        MPI_Send(&place, 1, MPI_INT, i + 1, TAG_STOP, MPI_COMM_WORLD);
//    }
}

void computeCycleSlave(double *aRows, double *aRows_1) {
    int place;
    MPI_Status status;
    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    int begin = rank + JN - 1;
    int end = JSIZE - (JSIZE - begin) % JN;

    printf("[SLAVE #%d] begin = %d, end = %d\n", rank, begin, end);
    fflush(stdout);

    while (true) {
        MPI_Recv(&place, 1, MPI_INT, MASTER_RANK, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
        if (status.MPI_TAG == TAG_STOP) {
            break;
        }
        MPI_Recv(aRows, JSIZE, MPI_DOUBLE, MASTER_RANK, TAG_ARRAY_ELEMENT, MPI_COMM_WORLD, &status);
        for (int j = 0; j < JSIZE; j++) {
            aRows[j] = 10 * place + j;
        }
        MPI_Send(&place, 1, MPI_INT, MASTER_RANK, TAG_PLACE, MPI_COMM_WORLD);
        MPI_Send(aRows, JSIZE, MPI_DOUBLE, MASTER_RANK, TAG_ARRAY_ELEMENT, MPI_COMM_WORLD);
    }
    printf("[SLAVE #%d] finished 1st cycle\n", rank);
    fflush(stdout);

    while (true) {
        MPI_Recv(&place, 1, MPI_INT, MASTER_RANK, MPI_ANY_TAG, MPI_COMM_WORLD, &status);

        printf("[SLAVE #%d] received place %d, status = %d\n", rank, place, status.MPI_TAG);
        fflush(stdout);

        if (status.MPI_TAG == TAG_STOP) {
            break;
        }
        MPI_Recv(aRows, JSIZE, MPI_DOUBLE, MASTER_RANK, TAG_ARRAY_ELEMENT, MPI_COMM_WORLD, &status);
        MPI_Recv(aRows_1, JSIZE, MPI_DOUBLE, MASTER_RANK, TAG_ARRAY_ELEMENT, MPI_COMM_WORLD, &status);
        for (int j = begin; j < end; j += 3) {
            aRows[j] = sin(0.00001 * aRows_1[j]);
        }
        MPI_Send(&place, 1, MPI_INT, MASTER_RANK, TAG_PLACE, MPI_COMM_WORLD);
        MPI_Send(aRows, JSIZE, MPI_DOUBLE, MASTER_RANK, TAG_ARRAY_ELEMENT, MPI_COMM_WORLD);
    }

    printf("[SLAVE #%d] finished 2nd cycle\n", rank);
    fflush(stdout);
}

// computed by only one thread
// as it is written in PDF
void computeSolo(double **a) {
    int i, j;
    for (i = 0; i < ISIZE; i++) {
        for (j = 0; j < JSIZE; j++) {
            a[i][j] = 10 * i + j;
        }
    }
    for (i = 1; i < ISIZE; i++) {
        for (j = 3; j < JSIZE; j++) {
            a[i][j] = sin(0.00001 * a[i - 1][j - 3]);
        }
    }
}

void compareTwoFiles() {

}

// for (i = 1; i < ISIZE; i++) {
// for (j = 3; j < JSIZE; j++) {
// a[i][j] = sin(0.00001 * a[i-1][j-3];
int main(int argc, char **argv) {
    int numThreads;
    int rank;
    FILE *ff;

    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &numThreads);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    if (numThreads < 2) {
        exit(EXIT_FAILURE);
    }

    // do master things
    if (rank == MASTER_RANK) {
        printf("Number of threads = %d\n", numThreads);
        auto a = create2DArray();

        double startTime = MPI_Wtime();
        computeCycleMaster(a, numThreads);
        double endTime = MPI_Wtime();
        double timeDiff = endTime - startTime;
        printf("time parallel: %0.16f\n", timeDiff);

        ff = fopen("result_parallel.txt", "w");
        printInFile(ff, a);
        fclose(ff);

        startTime = MPI_Wtime();
        computeSolo(a);
        endTime = MPI_Wtime();
        timeDiff = endTime - startTime;
        printf("time solo: %0.16f\n", timeDiff);

        ff = fopen("result_solo.txt", "w");
        printInFile(ff, a);
        fclose(ff);

        free2DArray(a);
        compareTwoFiles();

    } else { // if you are a slave thread
        double *aRows = create1DArray(JSIZE);
        double *aRows_1 = create1DArray(JSIZE);

        computeCycleSlave(aRows, aRows_1);

        free1DArray(aRows);
        free1DArray(aRows_1);
    }

    MPI_Barrier(MPI_COMM_WORLD);
    MPI_Finalize();
    return 0;
}