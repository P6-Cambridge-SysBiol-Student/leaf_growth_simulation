//
// Created by finley on 07/04/23.
//

#ifndef FRAP_FITNESS_H
#define FRAP_FITNESS_H

#endif //FRAP_FITNESS_H


double** computeDeltaFourierCoeffs(int desiredNumFourierCoeffs) { /// "nonuniform discrete Fourier transform of type II (NUDFT-II)"
    /// mallocing the memory for the pointers to each row
    double **FourierCoeffs = (double **) malloc(desiredNumFourierCoeffs * sizeof(double *));
    for (int i = 0; i < desiredNumFourierCoeffs; i++) {
        FourierCoeffs[i] = (double *) malloc(2 * sizeof(double));
    }
    double polarCoords[nbo][2];

    for (int i = 0; i < nbo; i++) {
        Point &cell = pointsArray[i];
        polarCoords[i][0] = cell.disVec.magnitude(); ///the radius value
        polarCoords[i][1] = atan2(cell.disVec.yy, cell.disVec.xx); ///the theta value
    }

    double maxRadiusValue = 0; /// find max radius value to scale radius_n values to between 0 and 1
    for (int jj = 0; jj<nbo; jj++){
        if (polarCoords[jj][0] > maxRadiusValue){
            maxRadiusValue = polarCoords[jj][0];
        }
    }

    /// calculate the sin and cos components of the fourier coefficients
    for (int k = 0; k < desiredNumFourierCoeffs; k++) {
        double &realComp = FourierCoeffs[k][0];
        double &imgComp = FourierCoeffs[k][1];
        realComp = 0;
        imgComp = 0;

        if (k == 0) {
            for (int n = 0; n < nbo; n++) {
                double &radiusN = polarCoords[n][0];
                double &thetaN = polarCoords[n][1];

                realComp += 1.0/nbo * radiusN;

            }
        } else {
            for (int n = 0; n < nbo; n++) {
                double &radiusN = polarCoords[n][0];
                double &thetaN = polarCoords[n][1];
                realComp += 1.0/nbo * radiusN * cos(k * thetaN);
                imgComp += 1.0/nbo * radiusN * sin(k * thetaN);
            }
        }
    }
    return FourierCoeffs;
}

void printDeltaFourierCoeffs(double** inputFourierArray, int desiredNumOfFourierCoeffs){
    for (int m = 0; m<desiredNumOfFourierCoeffs; m++) {
        double &realValue = inputFourierArray[m][0];
        double &imgValue = inputFourierArray[m][1];
        {
            printf("Magnitude/Phase of coefficient %d: %f   %f\n",
                   m, sqrt(realValue*realValue + imgValue*imgValue), atan2(imgValue, realValue));
        }
    }
}

void reconstructShape(double** inputFourierArray, int desiredNumOfFourierCoeffs){
    double x, y;
    int numPoints = 3000; /// higher = smoother curve
    double dTheta = 2*M_PI / numPoints; /// stepsize for the curve, should be sized to loop once
    double currentTheta = 0;

    double &a_0_real = inputFourierArray[0][0];
    double a_0 = sqrt(a_0_real*a_0_real);

    double reconstructedPolarCoords[numPoints][2];


    /// f(t) = a_0 + Σ(a_n*cos(2*pi*n*t/T) + b_n*sin(2*pi*n*t/T))
    for (int n = 0; n < numPoints; n++) {
        double reconstructedRadiusN = 2*a_0; /// doubling first coeff to account for how this is the average

        for (int k = 1; k < desiredNumOfFourierCoeffs; k++) {
            double &realComp = inputFourierArray[k][0];
            double &imgComp = inputFourierArray[k][1];

            reconstructedRadiusN += realComp * cos(k * currentTheta) + imgComp * sin(k * currentTheta);
        }
        reconstructedPolarCoords[n][0] = reconstructedRadiusN;
        reconstructedPolarCoords[n][1] = currentTheta;
        currentTheta += dTheta;
    }

    drawSquare(xBound, yBound);
    glLineWidth(3);
    glBegin(GL_LINE_LOOP);
    for (int i = 0; i < numPoints; i++) {
        double theta = i * dTheta; /// Calculate the angle t for the current step
        x = reconstructedPolarCoords[i][0] * cos(theta);
        y = reconstructedPolarCoords[i][0] * sin(theta);
        glColor3f(1.0,1.0,1.0);
        glVertex2f(x, y);
    }
    glEnd();
}

void outputFourierToFile(double** inputFourierArray, int desiredNumOfFourierCoeffs, const char* filename) {
    FILE *file = fopen(filename, "w");
    if (file == NULL) {
        printf("Error opening file `%s'!\n", filename);
        exit(1);
    }

    fprintf(file, "Index,Real,Imaginary,Magnitude,Phase\n");
    for (int m = 0; m < desiredNumOfFourierCoeffs; m++) {
        double &realValue = inputFourierArray[m][0];
        double &imgValue = inputFourierArray[m][1];
        double magnitude = sqrt(realValue * realValue + imgValue * imgValue);
        double phase = atan2(imgValue, realValue);

        fprintf(file, "%d,%f,%f,%f,%f\n", m, realValue, imgValue, magnitude, phase);
    }

    fclose(file);
}



