#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <mpi.h>

#define ISIZE 10000
#define JSIZE 10000
#define MASTER_RANK 0
#define JN 5
#define IN 4

#define TAG_ROW_NUM 1
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
    double *a4 = create1DArray(JSIZE);

    double *a5 = create1DArray(JSIZE);
    double *a6 = create1DArray(JSIZE);
    double *a7 = create1DArray(JSIZE);
    double *a8 = create1DArray(JSIZE);

    double *a9 = create1DArray(JSIZE);
    double *a10 = create1DArray(JSIZE);
    double *a11 = create1DArray(JSIZE);
    double *a12 = create1DArray(JSIZE);

    double *a13 = create1DArray(JSIZE);
    double *a14 = create1DArray(JSIZE);
    double *a15 = create1DArray(JSIZE);
    double *a16 = create1DArray(JSIZE);

    double *a17 = create1DArray(JSIZE);
    double *a18 = create1DArray(JSIZE);
    double *a19 = create1DArray(JSIZE);
    double *a20 = create1DArray(JSIZE);

    //sending to count sinus
    for (int i = 0; i < ISIZE; i += 4) {
        if (i + 7 > ISIZE )
            break;
        for (int j = 0; j < JN; j++) {
            int i_4 = i + 4;
            MPI_Send(&j, 1, MPI_INT, j + 1, TAG_ARRAY_NUM, MPI_COMM_WORLD);

            MPI_Send(a[i_4], JSIZE, MPI_DOUBLE, j + 1, TAG_ARRAY_ELEMENT, MPI_COMM_WORLD);
            MPI_Send(a[i_4 + 1], JSIZE, MPI_DOUBLE, j + 1, TAG_ARRAY_ELEMENT, MPI_COMM_WORLD);
            MPI_Send(a[i_4 + 2], JSIZE, MPI_DOUBLE, j + 1, TAG_ARRAY_ELEMENT, MPI_COMM_WORLD);
            MPI_Send(a[i_4 + 3], JSIZE, MPI_DOUBLE, j + 1, TAG_ARRAY_ELEMENT, MPI_COMM_WORLD);
        }
//        printf("[MASTER] complete i = %d\n", i);

        for (int j = 0; j < JN; j++) {
            MPI_Recv(&place, 1, MPI_INT, MPI_ANY_SOURCE, TAG_ARRAY_NUM, MPI_COMM_WORLD, &status);
            if (place == 0) {
                MPI_Recv(a1, JSIZE, MPI_DOUBLE, status.MPI_SOURCE, TAG_ARRAY_ELEMENT, MPI_COMM_WORLD, &status);
                MPI_Recv(a2, JSIZE, MPI_DOUBLE, status.MPI_SOURCE, TAG_ARRAY_ELEMENT, MPI_COMM_WORLD, &status);
                MPI_Recv(a3, JSIZE, MPI_DOUBLE, status.MPI_SOURCE, TAG_ARRAY_ELEMENT, MPI_COMM_WORLD, &status);
                MPI_Recv(a4, JSIZE, MPI_DOUBLE, status.MPI_SOURCE, TAG_ARRAY_ELEMENT, MPI_COMM_WORLD, &status);
            } else if (place == 1) {
                MPI_Recv(a5, JSIZE, MPI_DOUBLE, status.MPI_SOURCE, TAG_ARRAY_ELEMENT, MPI_COMM_WORLD, &status);
                MPI_Recv(a6, JSIZE, MPI_DOUBLE, status.MPI_SOURCE, TAG_ARRAY_ELEMENT, MPI_COMM_WORLD, &status);
                MPI_Recv(a7, JSIZE, MPI_DOUBLE, status.MPI_SOURCE, TAG_ARRAY_ELEMENT, MPI_COMM_WORLD, &status);
                MPI_Recv(a8, JSIZE, MPI_DOUBLE, status.MPI_SOURCE, TAG_ARRAY_ELEMENT, MPI_COMM_WORLD, &status);
            } else if (place == 2) {
                MPI_Recv(a9, JSIZE, MPI_DOUBLE, status.MPI_SOURCE, TAG_ARRAY_ELEMENT, MPI_COMM_WORLD, &status);
                MPI_Recv(a10, JSIZE, MPI_DOUBLE, status.MPI_SOURCE, TAG_ARRAY_ELEMENT, MPI_COMM_WORLD, &status);
                MPI_Recv(a11, JSIZE, MPI_DOUBLE, status.MPI_SOURCE, TAG_ARRAY_ELEMENT, MPI_COMM_WORLD, &status);
                MPI_Recv(a12, JSIZE, MPI_DOUBLE, status.MPI_SOURCE, TAG_ARRAY_ELEMENT, MPI_COMM_WORLD, &status);
            } else if (place == 3) {
                MPI_Recv(a13, JSIZE, MPI_DOUBLE, status.MPI_SOURCE, TAG_ARRAY_ELEMENT, MPI_COMM_WORLD, &status);
                MPI_Recv(a14, JSIZE, MPI_DOUBLE, status.MPI_SOURCE, TAG_ARRAY_ELEMENT, MPI_COMM_WORLD, &status);
                MPI_Recv(a15, JSIZE, MPI_DOUBLE, status.MPI_SOURCE, TAG_ARRAY_ELEMENT, MPI_COMM_WORLD, &status);
                MPI_Recv(a16, JSIZE, MPI_DOUBLE, status.MPI_SOURCE, TAG_ARRAY_ELEMENT, MPI_COMM_WORLD, &status);
            } else if (place == 4) {
                MPI_Recv(a17, JSIZE, MPI_DOUBLE, status.MPI_SOURCE, TAG_ARRAY_ELEMENT, MPI_COMM_WORLD, &status);
                MPI_Recv(a18, JSIZE, MPI_DOUBLE, status.MPI_SOURCE, TAG_ARRAY_ELEMENT, MPI_COMM_WORLD, &status);
                MPI_Recv(a19, JSIZE, MPI_DOUBLE, status.MPI_SOURCE, TAG_ARRAY_ELEMENT, MPI_COMM_WORLD, &status);
                MPI_Recv(a20, JSIZE, MPI_DOUBLE, status.MPI_SOURCE, TAG_ARRAY_ELEMENT, MPI_COMM_WORLD, &status);
            }
        }
//        printf("[MASTER] 2nd cycle i = %d\n", i);
//        fflush(stdout);

        for (int j = JN; j < JSIZE; j++) {
            switch (j % JN) {
                case 0: {
                    a[i][j] = a1[j];
                    a[i + 1][j] = a2[j];
                    a[i + 2][j] = a3[j];
                    a[i + 3][j] = a4[j];
                }
                case 1: {
                    a[i][j] = a5[j];
                    a[i + 1][j] = a6[j];
                    a[i + 2][j] = a7[j];
                    a[i + 3][j] = a8[j];
                }
                case 2: {
                    a[i][j] = a9[j];
                    a[i + 1][j] = a10[j];
                    a[i + 2][j] = a11[j];
                    a[i + 3][j] = a12[j];
                }
                case 3: {
                    a[i][j] = a13[j];
                    a[i + 1][j] = a14[j];
                    a[i + 2][j] = a15[j];
                    a[i + 3][j] = a16[j];
                }
                case 4: {
                    a[i][j] = a17[j];
                    a[i + 1][j] = a18[j];
                    a[i + 2][j] = a19[j];
                    a[i + 3][j] = a20[j];
                }
            }
        }
    }
//    printf("[MASTER] DONE\n");

    for (int j = 0; j < JN; j++) {
        MPI_Send(&rowNum, 1, MPI_INT, j + 1, TAG_STOP, MPI_COMM_WORLD);
    }

    free(a1);
    free(a2);
    free(a3);
    free(a4);

    free(a5);
    free(a6);
    free(a7);
    free(a8);

    free(a9);
    free(a10);
    free(a11);
    free(a12);

    free(a13);
    free(a14);
    free(a15);
    free(a16);

    free(a17);
    free(a18);
    free(a19);
    free(a20);
}

void computeCyclesSlave(double *aRows) {
    double *aRowsPlus4 = create1DArray(JSIZE);
    double *aRowsPlus5 = create1DArray(JSIZE);
    double *aRowsPlus6 = create1DArray(JSIZE);
    double *aRowsPlus7 = create1DArray(JSIZE);

    double *aRows1 = create1DArray(JSIZE);
    double *aRows2 = create1DArray(JSIZE);
    double *aRows3 = create1DArray(JSIZE);
    double *aRows4 = create1DArray(JSIZE);

    MPI_Status status;
    int rowNum, place;

    while (true) {
        MPI_Recv(&rowNum, 1, MPI_INT, MASTER_RANK, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
        if (status.MPI_TAG == TAG_STOP) {
            break;
        }
        MPI_Recv(aRows, JSIZE, MPI_DOUBLE, MASTER_RANK, TAG_ARRAY_ELEMENT, MPI_COMM_WORLD, &status);
        for (int j = 0; j < JSIZE; j++) {
            aRows[j] = 10 * rowNum + j;
        }
        MPI_Send(&rowNum, 1, MPI_INT, MASTER_RANK, TAG_ROW_NUM, MPI_COMM_WORLD);
        MPI_Send(aRows, JSIZE, MPI_DOUBLE, MASTER_RANK, TAG_ARRAY_ELEMENT, MPI_COMM_WORLD);
    }

    while (true) {
        MPI_Recv(&place, 1, MPI_INT, MASTER_RANK, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
        if (status.MPI_TAG == TAG_STOP) {
//            printf("[SLAVE] received stop\n");
            break;
        }
        MPI_Recv(aRowsPlus4, JSIZE, MPI_DOUBLE, MASTER_RANK, TAG_ARRAY_ELEMENT, MPI_COMM_WORLD, &status);
        MPI_Recv(aRowsPlus5, JSIZE, MPI_DOUBLE, MASTER_RANK, TAG_ARRAY_ELEMENT, MPI_COMM_WORLD, &status);
        MPI_Recv(aRowsPlus6, JSIZE, MPI_DOUBLE, MASTER_RANK, TAG_ARRAY_ELEMENT, MPI_COMM_WORLD, &status);
        MPI_Recv(aRowsPlus7, JSIZE, MPI_DOUBLE, MASTER_RANK, TAG_ARRAY_ELEMENT, MPI_COMM_WORLD, &status);

        for (int j = place; j < JSIZE; j += JN) {
            aRows1[j] = sin(0.00001 * aRowsPlus4[j - JN]);
            aRows2[j] = sin(0.00001 * aRowsPlus5[j - JN]);
            aRows3[j] = sin(0.00001 * aRowsPlus6[j - JN]);
            aRows4[j] = sin(0.00001 * aRowsPlus7[j - JN]);
        }
        MPI_Send(&place, 1, MPI_INT, MASTER_RANK, TAG_ARRAY_NUM, MPI_COMM_WORLD);

        MPI_Send(aRows1, JSIZE, MPI_DOUBLE, MASTER_RANK, TAG_ARRAY_ELEMENT, MPI_COMM_WORLD);
        MPI_Send(aRows2, JSIZE, MPI_DOUBLE, MASTER_RANK, TAG_ARRAY_ELEMENT, MPI_COMM_WORLD);
        MPI_Send(aRows3, JSIZE, MPI_DOUBLE, MASTER_RANK, TAG_ARRAY_ELEMENT, MPI_COMM_WORLD);
        MPI_Send(aRows4, JSIZE, MPI_DOUBLE, MASTER_RANK, TAG_ARRAY_ELEMENT, MPI_COMM_WORLD);
    }
    free1DArray(aRowsPlus4);
    free1DArray(aRowsPlus5);
    free1DArray(aRowsPlus6);
    free1DArray(aRowsPlus7);

    free1DArray(aRows1);
    free1DArray(aRows2);
    free1DArray(aRows3);
    free1DArray(aRows4);
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
    for (i = 0; i < ISIZE - 4; i++) {
        for (j = JN; j < JSIZE; j++) {
            a[i][j] = sin(0.00001 * a[i + 4][j - JN]);
        }
    }
}

// mpiexec -n 4 Lab2g
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
        computeCycleMaster(a, numThreads);
        double endTime = MPI_Wtime();
        double timeDiff = endTime - startTime;
        printf("time parallel: %0.16f\n", timeDiff);

        ff = fopen("result_parallel_2.txt", "w");
        printInFile(ff, a);
        fclose(ff);

        startTime = MPI_Wtime();
        computeSolo(a);
        endTime = MPI_Wtime();
        timeDiff = endTime - startTime;
        printf("time solo: %0.16f\n", timeDiff);

        ff = fopen("result_solo_2.txt", "w");
        printInFile(ff, a);
        fclose(ff);

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