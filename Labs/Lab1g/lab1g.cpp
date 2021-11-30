#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <mpi.h>
#include <algorithm>

#define ISIZE 1000
#define JSIZE 1000
#define MASTER_RANK 0
#define IN 3
#define JN 3

#define TAG_ROW_NUM 1
#define TAG_ARRAY_ELEMENT 2
#define TAG_STOP 3
#define TAG_ARRAY_NUM 4

void printInFile(FILE *ff, double **a) {
    for (int i = 0; i < ISIZE; i++) {
        for (int j = 0; j < JSIZE; j++) {
            fprintf(ff, "%e ", a[i][j]);
        }
        fprintf(ff, "\n");
    }
}

double *create1DArray(int size) {
    auto *a = static_cast<double *>(calloc(size, sizeof(double)));
    return a;
}

double **create2DArray() {
    auto **a = static_cast<double **>(calloc(ISIZE, sizeof(double *)));
    for (int i = 0; i < ISIZE; i++) {
        a[i] = create1DArray(JSIZE);
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
        MPI_Send(&i, 1, MPI_INT, i + 1, TAG_ROW_NUM, MPI_COMM_WORLD);
        MPI_Send(a[i], JSIZE, MPI_DOUBLE, i + 1, TAG_ARRAY_ELEMENT, MPI_COMM_WORLD);
    }

    for (int i = numThreads - 1; i < ISIZE; i++) {
        MPI_Recv(&rowNum, 1, MPI_INT, MPI_ANY_SOURCE, TAG_ROW_NUM, MPI_COMM_WORLD, &status);
        MPI_Recv(a[rowNum], JSIZE, MPI_DOUBLE, status.MPI_SOURCE, TAG_ARRAY_ELEMENT, MPI_COMM_WORLD, &status);

        MPI_Send(&i, 1, MPI_INT, status.MPI_SOURCE, TAG_ROW_NUM, MPI_COMM_WORLD);
        MPI_Send(a[i], JSIZE, MPI_DOUBLE, status.MPI_SOURCE, TAG_ARRAY_ELEMENT, MPI_COMM_WORLD);
    }

    for (int i = 0; i < numThreads - 1; i++) {
        MPI_Recv(&rowNum, 1, MPI_INT, MPI_ANY_SOURCE, TAG_ROW_NUM, MPI_COMM_WORLD, &status);
        MPI_Recv(a[rowNum], JSIZE, MPI_DOUBLE, status.MPI_SOURCE, TAG_ARRAY_ELEMENT, MPI_COMM_WORLD, &status);

        MPI_Send(&rowNum, 1, MPI_INT, status.MPI_SOURCE, TAG_STOP, MPI_COMM_WORLD);
    }

    double *a1 = create1DArray(JSIZE);
    double *a2 = create1DArray(JSIZE);
    double *a3 = create1DArray(JSIZE);

    int jMax = std::min(numThreads - 1, JN);

    for (int i = 1; i < ISIZE; i++) {
        for (int j = 0; j < jMax; j++) {
            MPI_Send(&j, 1, MPI_INT, j + 1, TAG_ARRAY_NUM, MPI_COMM_WORLD);
            MPI_Send(a[i - 1], JSIZE, MPI_DOUBLE, j + 1, TAG_ARRAY_ELEMENT, MPI_COMM_WORLD);
        }
        for (int j = jMax; j < JN; j++) {
            MPI_Recv(&place, 1, MPI_INT, MPI_ANY_SOURCE, TAG_ARRAY_NUM, MPI_COMM_WORLD, &status);
            if (place == 0) {
                MPI_Recv(a1, JSIZE, MPI_DOUBLE, status.MPI_SOURCE, TAG_ARRAY_ELEMENT, MPI_COMM_WORLD, &status);
            } else if (place == 1) {
                MPI_Recv(a2, JSIZE, MPI_DOUBLE, status.MPI_SOURCE, TAG_ARRAY_ELEMENT, MPI_COMM_WORLD, &status);
            } else if (place == 2) {
                MPI_Recv(a3, JSIZE, MPI_DOUBLE, status.MPI_SOURCE, TAG_ARRAY_ELEMENT, MPI_COMM_WORLD, &status);
            }
            MPI_Send(&j, 1, MPI_INT, status.MPI_SOURCE, TAG_ARRAY_NUM, MPI_COMM_WORLD);
            MPI_Send(a[i - 1], JSIZE, MPI_DOUBLE, status.MPI_SOURCE, TAG_ARRAY_ELEMENT, MPI_COMM_WORLD);
        }
        for (int j = 0; j < jMax; j++) {
            MPI_Recv(&place, 1, MPI_INT, MPI_ANY_SOURCE, TAG_ARRAY_NUM, MPI_COMM_WORLD, &status);
            if (place == 0) {
                MPI_Recv(a1, JSIZE, MPI_DOUBLE, status.MPI_SOURCE, TAG_ARRAY_ELEMENT, MPI_COMM_WORLD, &status);
            } else if (place == 1) {
                MPI_Recv(a2, JSIZE, MPI_DOUBLE, status.MPI_SOURCE, TAG_ARRAY_ELEMENT, MPI_COMM_WORLD, &status);
            } else if (place == 2) {
                MPI_Recv(a3, JSIZE, MPI_DOUBLE, status.MPI_SOURCE, TAG_ARRAY_ELEMENT, MPI_COMM_WORLD, &status);
            }

        }

        for (int j = JN; j < JSIZE; j++) {
            if (j % JN == 0) {
                a[i][j] = a1[j];
            } else if (j % JN == 1) {
                a[i][j] = a2[j];
            } else if (j % JN == 2) {
                a[i][j] = a3[j];
            }
        }
    }

    for (int i = 0; i < numThreads - 1; i++) {
        MPI_Send(&place, 1, MPI_INT, i + 1, TAG_STOP, MPI_COMM_WORLD);
    }
    free(a1);
    free(a2);
    free(a3);
}

void computeCycleSlave(double *aRows, double *aRows_1) {
    int place;
    int rowNum;
    MPI_Status status;

    while (true) {
        MPI_Recv(&rowNum, 1, MPI_INT, MASTER_RANK, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
        if (status.MPI_TAG == TAG_STOP) {
            break;
        }
        MPI_Recv(aRows, JSIZE, MPI_DOUBLE, MASTER_RANK, TAG_ARRAY_ELEMENT, MPI_COMM_WORLD, &status); // aRows?
        for (int j = 0; j < JSIZE; j++) {
            aRows[j] = 10 * rowNum + j;
        }
        MPI_Send(&rowNum, 1, MPI_INT, MASTER_RANK, TAG_ROW_NUM, MPI_COMM_WORLD);
        MPI_Send(aRows, JSIZE, MPI_DOUBLE, MASTER_RANK, TAG_ARRAY_ELEMENT, MPI_COMM_WORLD);
    }

    while (true) {
        MPI_Recv(&place, 1, MPI_INT, MASTER_RANK, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
        if (status.MPI_TAG == TAG_STOP) {
            break;
        }

        int begin = place + JN;

        MPI_Recv(aRows_1, JSIZE, MPI_DOUBLE, MASTER_RANK, TAG_ARRAY_ELEMENT, MPI_COMM_WORLD, &status);
        for (int j = begin; j < JSIZE; j += JN) {
            aRows[j] = sin(0.00001 * aRows_1[j - JN]);
        }
        MPI_Send(&place, 1, MPI_INT, MASTER_RANK, TAG_ARRAY_NUM, MPI_COMM_WORLD);
        MPI_Send(aRows, JSIZE, MPI_DOUBLE, MASTER_RANK, TAG_ARRAY_ELEMENT, MPI_COMM_WORLD);
    }
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
        for (j = JN; j < JSIZE; j++) {
            a[i][j] = sin(0.00001 * a[i - 1][j - JN]);
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
