//
// Created by finley on 29/03/23.
//

#include <stdio.h>
#include <stdlib.h>

#ifndef FRAP_WRITING_H
#define FRAP_WRITING_H

#endif //FRAP_WRITING_H

// TODO these functions need to output the coefficients into the

void exampleWritingFunction(){
    FILE *fptr;
    fptr = fopen("example_param_file", "w");
    fprintf(fptr, "This is text that will overwrite what is in this file.\n");
    fclose(fptr);
}




