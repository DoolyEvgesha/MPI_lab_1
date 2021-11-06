Start with:

    mpiexec -np <num_threads> [LabBasic/Lab1g]

## Basic Task:
    #define ISIZE 1000
    #define JSIZE 1000
    // bla bla bla main etc
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

# Task 1G:

    #define ISIZE 1000
    #define JSIZE 1000
    // bla bla bla main etc
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


# Task 2G:

    #define ISIZE 1000
    #define JSIZE 1000
    // bla bla bla main etc
    int i, j;

    for (i = 0; i < ISIZE; i++) {

        for (j = 0; j < JSIZE; j++) {

            a[i][j] = 10 * i + j;

        }

    }

    for (i = 0; i < ISIZE - 4; i++) {

        for (j = 5; j < JSIZE; j++) {

            a[i][j] = sin(0.00001 * a[i + 4][j - 5]);

        }

    }
    // write in the file

Source: https://drive.google.com/drive/folders/1n1AwSxdBu9_r_iQixB8yDj1UvbOkRyDf

Лабораторная_работа.pdf

P.S. in case I complete it and it helps you, you owe me a beer :)