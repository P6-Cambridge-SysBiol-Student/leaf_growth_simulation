//
// Created by finley on 13/02/23.
//

#ifndef FRAP_SIGMOID_H
#define FRAP_SIGMOID_H

#endif //FRAP_SIGMOID_H

#ifndef DEF_SIGMOID
#define DEF_SIGMOID

#define EULER_NUMBER 2.71828
#define EULER_NUMBER_F 2.71828182846
#define EULER_NUMBER_L 2.71828182845904523536

#endif
#include <math.h>

double sigmoid(double n) {
    return (1 / (1 + pow(EULER_NUMBER, -n)));
}

float sigmoidf(float n) {
    return (1 / (1 + powf(EULER_NUMBER_F, -n)));
}

long double sigmoidl(long double n) {
    return (1 / (1 + powl(EULER_NUMBER_L, -n)));
}