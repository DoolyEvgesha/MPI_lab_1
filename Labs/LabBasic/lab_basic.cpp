#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <mpi.h>

#define ISIZE 1000
#define JSIZE 1000
#define MASTER_RANK 0

#define TAG_ROW_NUM 1
#define TAG_ARRAY_ELEMENT 2
#define TAG_STOP 3

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
    int rowNum;
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
    // for sinus
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
}

void computeCyclesSlave(double *a) {
    MPI_Status status;
    int numRow;

    while (true) {
        MPI_Recv(&numRow, 1, MPI_INT, MASTER_RANK, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
        if (status.MPI_TAG == TAG_STOP) {
            break;
        }
        MPI_Recv(a, JSIZE, MPI_DOUBLE, MASTER_RANK, TAG_ARRAY_ELEMENT, MPI_COMM_WORLD, &status);
        for (int j = 0; j < JSIZE; j++) {
            a[j] = 10 * numRow + j;
        }
        MPI_Send(&numRow, 1, MPI_INT, MASTER_RANK, TAG_ROW_NUM, MPI_COMM_WORLD);
        MPI_Send(a, JSIZE, MPI_DOUBLE, MASTER_RANK, TAG_ARRAY_ELEMENT, MPI_COMM_WORLD);
    }

    while (true) {
        MPI_Recv(&numRow, 1, MPI_INT, MASTER_RANK, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
        if (status.MPI_TAG == TAG_STOP) {
            break;
        }
        MPI_Recv(a, JSIZE, MPI_DOUBLE, MASTER_RANK, TAG_ARRAY_ELEMENT, MPI_COMM_WORLD, &status);
        for (int j = 0; j < JSIZE; j++) {
            a[j] = sin(0.00001 * a[j]);
        }
        MPI_Send(&numRow, 1, MPI_INT, MASTER_RANK, TAG_ROW_NUM, MPI_COMM_WORLD);
        MPI_Send(a, JSIZE, MPI_DOUBLE, MASTER_RANK, TAG_ARRAY_ELEMENT, MPI_COMM_WORLD);
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
    for (i = 0; i < ISIZE; i++) {
        for (j = 0; j < JSIZE; j++) {
            a[i][j] = sin(0.00001 * a[i][j]);
        }
    }
}

// mpiexec -np 4 LabBasic
int main(int argc, char **argv) {
    int numThreads;
    int rank;
    FILE *ff;

    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &numThreads);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    fflush(stdout);
    // do master deeds
    if (rank == MASTER_RANK) {
        auto a = create2DArray();

        double startTime = MPI_Wtime();
        //two cycles as in compute solo ones
        computeCycleMaster(a, numThreads);
        double endTime = MPI_Wtime();
        double timeDiff = endTime - startTime;
        printf("time parallel: %0.16f\n", timeDiff);

        ff = fopen("result.txt", "w");
        printInFile(ff, a);
        fclose(ff);

        startTime = MPI_Wtime();
        computeSolo(a);
        endTime = MPI_Wtime();
        timeDiff = endTime - startTime;
        printf("time solo: %0.16f\n", timeDiff);

        free2DArray(a);

    } else { // if you are a slave thread
        double *aRows = create1DArray(JSIZE);
        computeCyclesSlave(aRows);

        free1DArray(aRows);
    }

    MPI_Barrier(MPI_COMM_WORLD);
    MPI_Finalize();
    return 0;
}

