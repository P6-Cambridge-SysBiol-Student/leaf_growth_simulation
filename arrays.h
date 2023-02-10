//
// Created by finley on 10/02/23.
//

#ifndef FRAP_ARRAYS_H
#define FRAP_ARRAYS_H
#include <stdlib.h>
#include <string.h>
#endif //FRAP_ARRAYS_H

int** create2Darray(int rows, int cols) {  /// here rows  = nbo and cols = NAW

    /// mallocing the memory for the pointers to each row
    int **arr = (int **) malloc(rows * sizeof(int *));
    /// mallocing memory for each column
    for (int i = 0; i < rows; i++) {
        arr[i] = (int *) malloc(cols * sizeof(int));
    }
    return arr;
}

int* create1Darray(int size) {
    int* array = (int*) malloc(size * sizeof(int));
    return array;
}

/// fill neighbourhood array with -1's
void init2DArray(int** inputArrayPointer, int rows, int cols, int inputFillValue){
    for (int i = 0; i < rows; i++) {
        memset(inputArrayPointer[i], inputFillValue, cols * sizeof(int));
    }
}

void init1DArray(int* input1DArrayPointer, int size, int inputFillValue){
    memset(input1DArrayPointer, inputFillValue, size * sizeof(int));
}


