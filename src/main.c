#include <stdio.h>
#include <mpi.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

const int MAIN_PID = 0;
const int DEFAULT_TAG = 0;

struct HeatSource {
    int upperLeftRow;
    int upperLeftColumn;
    int lowerRightRow;
    int lowerRightColumn;
    double temp;
};

void inputHeatSource(struct HeatSource *heatSource) {
    int upperLeftRow;
    int upperLeftColumn;
    int lowerRightRow;
    int lowerRightColumn;
    double temp;

    printf("Upper left row:\n");
    scanf("%i", &upperLeftRow);

    printf("Upper left column:\n");
    scanf("%i", &upperLeftColumn);

    printf("Lower right row:\n");
    scanf("%i", &lowerRightRow);

    printf("Lower right column:\n");
    scanf("%i", &lowerRightColumn);

    printf("Temp (0-1):\n");
    scanf("%lf", &temp);

    heatSource->upperLeftRow = upperLeftRow;
    heatSource->upperLeftColumn = upperLeftColumn;
    heatSource->lowerRightRow = lowerRightRow;
    heatSource->lowerRightColumn = lowerRightColumn;
    heatSource->temp = temp;
}

void inputData(int *m, int *n, int *numberOfHeatSources, struct HeatSource **heatSources) {

    printf("M:\n");
    scanf("%i", m);

    printf("N:\n");
    scanf("%i", n);

    printf("Number of heat sources:\n");
    scanf("%i", numberOfHeatSources);

    *heatSources = calloc((size_t) *numberOfHeatSources, sizeof(struct HeatSource));

    int i;
    for (i = 0; i < *numberOfHeatSources; i++) {
        printf("Heat source %i:\n", i);
        struct HeatSource heatSource;
        inputHeatSource(&heatSource);
        (*heatSources)[i] = heatSource;
    }
}

int getPid() {
    int pid;
    MPI_Comm_rank(MPI_COMM_WORLD, &pid);
    return pid;
}

int getNumberOfProcesses() {
    int numberOfProcesses;
    MPI_Comm_size(MPI_COMM_WORLD, &numberOfProcesses);
    return numberOfProcesses;
}

void fillSheetWithHeatSources(double *sheet, int m, int n, int numberOfHeatSources, struct HeatSource *heatSources) {
    int i;
    for (i = 0; i < numberOfHeatSources; i++) {
//        printf("Heat source %i\n", i);

        struct HeatSource heatSource = heatSources[i];

        int row;
        int column;
        for (row = heatSource.upperLeftRow; row <= heatSource.lowerRightRow; row++) {
            for (column = heatSource.upperLeftColumn; column <= heatSource.lowerRightColumn; column++) {
                *(sheet + column * n + row) = heatSource.temp;
//                printf("r: %i, c: %i, temp: %lf\n", row, column, heatSource.temp);
            }
        }

        for (row = 1; row < m - 1; row++) {
            for (column = 1; column < n - 1; column++) {
                *(sheet + column * n + row) = 0;
            }
        }
    }
}

void printSheet(double *sheet, int numberOfCells, int n) {
    int i;
    for (i = 0; i < numberOfCells; i++) {
        if (i % n == 0) {
            printf("\n");
        }
        printf("%lf ", sheet[i]);
    }
    printf("\n");
}

void printHeatSources(struct HeatSource *heatSources, int numberOfHeatSources) {
    int i;
    for (i = 0; i < numberOfHeatSources; i++) {
        printf("ulr: %i, ulc: %i, lrr: %i, lrc: %o, temp: %lf\n", heatSources[i].upperLeftRow,
               heatSources[i].upperLeftColumn, heatSources[i].lowerRightRow, heatSources[i].lowerRightColumn,
               heatSources[i].temp);
    }
}

void printParams(int m, int n, int numberOfHeatSources, struct HeatSource *heatSources) {
    printf("m: %i, n: %i, h: %i\n", m, n, numberOfHeatSources);

    printHeatSources(heatSources, numberOfHeatSources);

}

void run(int pid, int numberOfProcesses) {
    int m;
    int n;

    int mWithAir;
    int nWithAir;

    int numberOfHeatSources;
    struct HeatSource *heatSources;

    int numberOfCellsPerProcess;
    int numberOfRowsPerProcess;
    int numberOfCells;

    double *sheet;

    double *mySheet;
    double *myOldSheet;

    if (pid == MAIN_PID) {
        inputData(&m, &n, &numberOfHeatSources, &heatSources);
    }

    MPI_Bcast(&n, 1, MPI_INT, MAIN_PID, MPI_COMM_WORLD);
    MPI_Bcast(&m, 1, MPI_INT, MAIN_PID, MPI_COMM_WORLD);

    mWithAir = m + 2;
    nWithAir = n + 2;

    numberOfRowsPerProcess = m / numberOfProcesses;
    numberOfCellsPerProcess = (numberOfRowsPerProcess + 2) * nWithAir;

    if (pid == MAIN_PID) {
        numberOfCells = mWithAir * nWithAir;
        sheet = calloc((size_t) numberOfCells, sizeof(double));

        fillSheetWithHeatSources(sheet, mWithAir, nWithAir, numberOfHeatSources, heatSources);

//        printf("norpp: %i\n", numberOfRowsPerProcess);

//        printParams(m, n, numberOfHeatSources, heatSources);
        printSheet(sheet, numberOfCells, nWithAir);

        int i;
        for (i = 0; i < numberOfProcesses; i++) {
            MPI_Send((sheet + ((i * numberOfRowsPerProcess) * nWithAir)), numberOfCellsPerProcess, MPI_DOUBLE, i,
                     DEFAULT_TAG,
                     MPI_COMM_WORLD);
        }
    }

    mySheet = calloc((size_t) numberOfCellsPerProcess, sizeof(double));
    myOldSheet = calloc((size_t) numberOfCellsPerProcess, sizeof(double));
    MPI_Recv(mySheet, numberOfCellsPerProcess, MPI_DOUBLE, MAIN_PID, DEFAULT_TAG, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

    int i = 0;
    int hasChanged = 1;
    while (hasChanged == 1) {
        memcpy(myOldSheet, mySheet, numberOfCellsPerProcess * sizeof(double));

        printSheet(myOldSheet, numberOfCellsPerProcess, nWithAir);

        int row;
        int col;
        for (row = 1; row < numberOfRowsPerProcess + 1; row++) {
            for (col = 1; col < nWithAir - 1; col++) {
                double rowAbove = *(myOldSheet + (((row + 1) * nWithAir) + col));
                double rowBelow = *(myOldSheet + (((row - 1) * nWithAir) + col));
                double colRight = *(myOldSheet + ((row * nWithAir) + (col + 1)));
                double colLeft = *(myOldSheet + ((row * nWithAir) + (col - 1)));
                double curr = *(myOldSheet + ((row * nWithAir) + col));
                double sum = (rowAbove + rowBelow + colRight + colLeft);
                double new =  sum / 4;
                printf("i %i, r %i, c %i, curr %lf, sum %lf, new %lf, ra %lf, rb %lf, cr %lf, cl %lf\n", i, row, col, curr, sum, new, rowAbove, rowBelow, colRight,
                       colLeft);
                *(mySheet + ((row * nWithAir) + col)) = new;
                if (fabs(curr - new) < .000001) {
                    hasChanged = 0;
                }
            }
        }

        if (pid == 0) {
            // Send below
            MPI_Send(&mySheet);
        } else if (pid == numberOfProcesses - 1) {
            // Send above
            MPI_Send();
        } else {
            // Send above and below
            MPI_Send();
            MPI_Send();
        }

//        printSheet(mySheet, numberOfCellsPerProcess, nWithAir);
        i += 1;
    }

//    printf("Printing mysheet %i\n", pid);

//    int i;
//    for (i = 0; i < numberOfProcesses; i++) {
//        if (pid == i) {
//            printf("pid %i\n", pid);
//            printSheet(mySheet, numberOfCellsPerProcess, nWithAir);
//        }
//        MPI_Barrier(MPI_COMM_WORLD);
//    }
}

// https://www.geeksforgeeks.org/dynamically-allocate-2d-array-c/
int main(int argc, char **argv) {
    MPI_Init(&argc, &argv);

    int pid = getPid();
    int numberOfProcesses = getNumberOfProcesses();

    run(pid, numberOfProcesses);

    MPI_Finalize();
}
