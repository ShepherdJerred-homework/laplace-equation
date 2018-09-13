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

//    printf("Upper left row:\n");
    scanf("%i", &upperLeftColumn);

//    printf("Upper left column:\n");
    scanf("%i", &upperLeftRow);

//    printf("Lower right row:\n");
    scanf("%i", &lowerRightColumn);

//    printf("Lower right column:\n");
    scanf("%i", &lowerRightRow);

//    printf("Temp (0-1):\n");
    scanf("%lf", &temp);

    heatSource->upperLeftRow = upperLeftRow;
    heatSource->upperLeftColumn = upperLeftColumn;
    heatSource->lowerRightRow = lowerRightRow;
    heatSource->lowerRightColumn = lowerRightColumn;
    heatSource->temp = temp;
}

void inputData(int *m, int *n, int *numberOfHeatSources, struct HeatSource **heatSources) {

//    printf("M:\n");
    scanf("%i", m);

//    printf("N:\n");
    scanf("%i", n);

//    printf("Number of heat sources:\n");
    scanf("%i", numberOfHeatSources);

    *heatSources = calloc((size_t) *numberOfHeatSources, sizeof(struct HeatSource));

    int i;
    for (i = 0; i < *numberOfHeatSources; i++) {
//        printf("Heat source %i:\n", i);
        struct HeatSource heatSource;
        inputHeatSource(&heatSource);
        (*heatSources)[i] = heatSource;
    }

//    printf("Input finished! Executing..\n");
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

void printMyRows(double *sheet, int numberOfRowsPerProcess, int nWithAir) {
    int row;
    int col;
    for (row = 1; row < numberOfRowsPerProcess + 1; row++) {
        for (col = 0; col < nWithAir; col++) {
            printf("%lf ", *(sheet + ((row * nWithAir) + col)));
        }
        printf("\n");
    }
}

void printFirstRow(double *sheet, int nWithAir) {
    int i;
    for (i = 0; i < nWithAir; i++) {
        printf("%lf ", *(sheet + i));
    }
    printf("\n");
}

void printLastRow(double *sheet, int numberOfRowsPerProcess, int pid, int nWithAir) {
    double *lastRow = (sheet + (((numberOfRowsPerProcess + 1) * nWithAir)));
    int i;
    for (i = 0; i < nWithAir; i++) {
        printf("%lf ", *(lastRow + i));
    }
    printf("\n");
}

void synchronize(int pid, int numberOfProcesses, int numberOfRowsPerProcess, double *sheet, int nWithAir) {
    double *lastRowBelongingToMe = (sheet + (numberOfRowsPerProcess * nWithAir));
    double *firstRowBelongingToMe = (sheet + nWithAir);

    double *firstRowInMySheet = (sheet);
    double *lastRowInMySheet = (sheet + (numberOfRowsPerProcess * nWithAir) + (nWithAir));

    if (pid == 0) {
        // Send to pid + 1, receive from pid + 1
        MPI_Send(lastRowBelongingToMe, nWithAir, MPI_DOUBLE, pid + 1, DEFAULT_TAG, MPI_COMM_WORLD);
        MPI_Recv(lastRowInMySheet, nWithAir, MPI_DOUBLE, pid + 1, DEFAULT_TAG, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    } else if (pid == numberOfProcesses - 1) {
        // Send to p - 1, receive from pid - 1
        MPI_Send(firstRowBelongingToMe, nWithAir, MPI_DOUBLE, pid - 1, DEFAULT_TAG, MPI_COMM_WORLD);
        MPI_Recv(firstRowInMySheet, nWithAir, MPI_DOUBLE, pid - 1, DEFAULT_TAG, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    } else {
        // Send to pid + 1, receive from pid + 1
        MPI_Send(lastRowBelongingToMe, nWithAir, MPI_DOUBLE, pid + 1, DEFAULT_TAG, MPI_COMM_WORLD);
        MPI_Recv(lastRowInMySheet, nWithAir, MPI_DOUBLE, pid + 1, DEFAULT_TAG, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        // Send to p - 1, receive from pid - 1
        MPI_Send(firstRowBelongingToMe, nWithAir, MPI_DOUBLE, pid - 1, DEFAULT_TAG, MPI_COMM_WORLD);
        MPI_Recv(firstRowInMySheet, nWithAir, MPI_DOUBLE, pid - 1, DEFAULT_TAG, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    }
}

void prettyPrintMySheet(int numberOfProcesses, int pid, double* mySheet, int nWithAir, int numberOfRowsPerProcess) {
    int i;
    for (i = 0; i < numberOfProcesses; i++) {
        if (pid == i && pid == 0) {
            printFirstRow(mySheet, nWithAir);
        }
        if (pid == i) {
            printMyRows(mySheet, numberOfRowsPerProcess, nWithAir);
        }
        if (pid == i && pid == numberOfProcesses - 1) {
            printLastRow(mySheet, numberOfRowsPerProcess, pid, nWithAir);
        }
        MPI_Barrier(MPI_COMM_WORLD);
    }
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

    memcpy(myOldSheet, mySheet, numberOfCellsPerProcess * sizeof(double));

    int i = 0;
    int shouldContinue = 1;
    int hasChanged;

    while (shouldContinue) {
        i += 1;
        hasChanged = 0;
        shouldContinue = 0;

        double *c;
        c = myOldSheet;
        myOldSheet = mySheet;
        mySheet = c;

        int row;
        int col;
        for (row = 1; row < numberOfRowsPerProcess + 1; row++) {
            for (col = 1; col < nWithAir - 1; col++) {
                double rowAbove = *(myOldSheet + (((row + 1) * nWithAir) + col));
                double rowBelow = *(myOldSheet + (((row - 1) * nWithAir) + col));
                double columnRight = *(myOldSheet + ((row * nWithAir) + (col + 1)));
                double columnLeft = *(myOldSheet + ((row * nWithAir) + (col - 1)));

                double currentValue = *(myOldSheet + ((row * nWithAir) + col));
                double sumOfAdjacent = (rowAbove + rowBelow + columnRight + columnLeft);
                double newValue = sumOfAdjacent / 4;
                double currentAndNewDifference = fabs(currentValue - newValue);
//                printf("i %i, r %i, c %i, currentValue %lf, sumOfAdjacent %lf, newValue %lf, currentAndNewDifference %lf, ra %lf, rb %lf, cr %lf, cl %lf\n", i,
//                       row, col, currentValue, sumOfAdjacent, newValue, currentAndNewDifference, rowAbove, rowBelow, columnRight,
//                       columnLeft);
                *(mySheet + ((row * nWithAir) + col)) = newValue;
                if (currentAndNewDifference > .000001) {
                    hasChanged = 1;
                }
            }
        }

        synchronize(pid, numberOfProcesses, numberOfRowsPerProcess, mySheet, nWithAir);
        MPI_Allreduce(&hasChanged, &shouldContinue, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
    }
    prettyPrintMySheet(numberOfProcesses, pid, mySheet, nWithAir, numberOfRowsPerProcess);
}

// https://www.geeksforgeeks.org/dynamically-allocate-2d-array-c/
int main(int argc, char **argv) {
    MPI_Init(&argc, &argv);

    int pid = getPid();
    int numberOfProcesses = getNumberOfProcesses();

    run(pid, numberOfProcesses);

    MPI_Finalize();
}
