#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <mpi.h>

#define ISIZE 1000
#define JSIZE 1000
#define MASTER_RANK 0

#define TAG_PLACE 1
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
}

void computeCyclesSlave(double *a) {
    MPI_Status status;
    int numRow;

    // receive from MASTER thread
    // for simple computation of 10 * i + j
    while (true) {
        MPI_Recv(&numRow, 1, MPI_INT, MASTER_RANK, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
        if (status.MPI_TAG == TAG_STOP) {
            break;
        }
        MPI_Recv(a, JSIZE, MPI_DOUBLE, MASTER_RANK, TAG_ARRAY_ELEMENT, MPI_COMM_WORLD, &status);
        for (int j = 0; j < JSIZE; j++) {
            a[j] = 10 * numRow + j;
        }
        MPI_Send(&numRow, 1, MPI_INT, MASTER_RANK, TAG_PLACE, MPI_COMM_WORLD);
        MPI_Send(a, JSIZE, MPI_DOUBLE, MASTER_RANK, TAG_ARRAY_ELEMENT, MPI_COMM_WORLD);
    }

    // for sinus
    while (true) {
        MPI_Recv(&numRow, 1, MPI_INT, MASTER_RANK, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
        if (status.MPI_TAG == TAG_STOP) {
            break;
        }
        MPI_Recv(a, JSIZE, MPI_DOUBLE, MASTER_RANK, TAG_ARRAY_ELEMENT, MPI_COMM_WORLD, &status);
        for (int j = 0; j < JSIZE; j++) {
            a[j] = sin(0.00001 * a[j]);
        }
        MPI_Send(&numRow, 1, MPI_INT, MASTER_RANK, TAG_PLACE, MPI_COMM_WORLD);
        MPI_Send(a, JSIZE, MPI_DOUBLE, MASTER_RANK, TAG_ARRAY_ELEMENT, MPI_COMM_WORLD);
    }

    // for sinus and (10 * i + j)
    while (true) {
        MPI_Recv(&numRow, 1, MPI_INT, MASTER_RANK, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
        if (status.MPI_TAG == TAG_STOP) {
            break;
        }
        MPI_Recv(a, JSIZE, MPI_DOUBLE, MASTER_RANK, TAG_ARRAY_ELEMENT, MPI_COMM_WORLD, &status);
        for (int j = 0; j < JSIZE; j++) {
            a[j] = sin(0.00001 * (10 * numRow + j));
        }
        MPI_Send(&numRow, 1, MPI_INT, MASTER_RANK, TAG_PLACE, MPI_COMM_WORLD);
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

// mpiexec -n 4 LabBasic
//C:\Users\darya\CLionProjects\Lab1MPIBasicTask\Labs\Lab1r\lab_1r.cpp
int main(int argc, char **argv) {
    int numThreads;
    int rank;
    FILE *ff;

    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &numThreads);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    printf("Number of threads = %d, rank = %d\n", numThreads, rank);
    fflush(stdout);
    // do master deeds
    if (rank == MASTER_RANK) {
        auto a = create2DArray();

        double startTime = MPI_Wtime();
        //two cycles as in compute solo ones
        computeCycleMaster(a, numThreads);
        computeCycleMaster(a, numThreads);
        double endTime = MPI_Wtime();
        double timeDiff = endTime - startTime;
        printf("time without changing graph: %0.16f\n", timeDiff);

        startTime = MPI_Wtime();
        //one cycle for both 10 * i + j and sin
        computeCycleMaster(a, numThreads);
        endTime = MPI_Wtime();
        timeDiff = endTime - startTime;
        printf("time with changing graph: %0.16f\n", timeDiff);

        ff = fopen("result.txt", "w");
        printInFile(ff, a);
        fclose(ff);

        // start region solo computational area
        startTime = MPI_Wtime();
        computeSolo(a);
        endTime = MPI_Wtime();
        timeDiff = endTime - startTime;
        printf("time solo: %0.16f\n", timeDiff);
        // end region

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

//Number of threads = 2
//time without changing graph: 0.0308983000000000
//time with changing graph: 0.0286849000000000
//time solo: 0.1230460000000000

//Number of threads = 3
//time without changing graph: 0.0279267000000000
//time with changing graph: 0.0129701000000000
//time solo: 0.1011228999999999

//Number of threads = 4
//time without changing graph: 0.0313012000000000
//time with changing graph: 0.0127284000000000
//time solo: 0.1231338000000001

//Number of threads = 5
//time without changing graph: 0.0397578000000000
//time with changing graph: 0.0141880000000000
//time solo: 0.1052784999999999

//Number of threads = 6
//time without changing graph: 0.0489594000000000
//time with changing graph: 0.0142033000000000
//time solo: 0.1237226000000000

//Number of threads = 7, rank = 6
//time without changing graph: 0.1680818000000000
//time with changing graph: 0.0097112000000000
//time solo: 0.0828555999999999

//Number of threads = 8
//time without changing graph: 0.1288778000000000
//time with changing graph: 0.0093467000000000
//time solo: 0.0180776000000000
