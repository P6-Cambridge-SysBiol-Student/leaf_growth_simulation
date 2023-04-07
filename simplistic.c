/*
 Simplistic C-program to read parameters from a file
 FJN 27.03.2023, Copyright 2023 Cambridge University
 */

#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <ctype.h>

double alpha = 0;
double beta = 0;
double zeta = 0;

//------------------------------------------------------------------------------
int readParameter(const char arg[], const char name[], double* ptr)
{
    int res = 0;
    char * dup = strdup(arg);
    char * val = dup;
    char * key = strsep(&val, "=");
    char * end = NULL;
    double d = 0;

    while ( isspace(*key) ) ++key;
    if ( *key )
    {
        if ( key == strstr(key, name) )
        {
            while ( isspace(*val) ) ++val;
            //printf("    found `%s' in [%s] : %s\n", name, key, val);
            d = strtod(val, &end);
            if ( end > val )
            {
                *ptr = d;
                res = 1;
            }
        }
    }
    free(dup);
    return res;
}

int readLine(const char arg[])
{
    if ( readParameter(arg, "alpha", &alpha) ) return 1;
    if ( readParameter(arg, "beta", &beta) ) return 1;
    if ( readParameter(arg, "zeta", &zeta) ) return 1;
    return 0;
}

void readFile(const char path[])
{
    FILE * file = NULL;
    char * line = NULL;
    size_t len = 0;
    ssize_t read = 0;
    
    file = fopen(path, "r");
    if ( !file ) {
        printf("Error: file `%s' cannot be found\n", path);
        return;
    }
    if ( ferror(file) ) {
        fclose(file);
        printf("Error: file `%s' cannot be read\n", path);
        return;
    }
    //printf("reading file [%s]\n", path);
    while ((read = getline(&line, &len, file)) != -1 )
    {
        line[read-1] = 0;
        //printf("  reading [%s]:\n", line);
        readLine(line);
    }
    free(line);
    fclose(file);
}


void writeFile(const char path[], double data)
{
    FILE * file = fopen(path, "w");
    if ( !file ) {
        printf("Error: file `%s' cannot be created\n", path);
        return;
    }
    if ( ferror(file) ) {
        fclose(file);
        printf("Error: file `%s' cannot be written\n", path);
        return;
    }
    fprintf(file, "%f\n", data);
    fclose(file);
}


/* program entry */
int main(int argc, char *argv[])
{
    for ( int i = 1; i < argc; ++i )
    {
        const char * arg = argv[i];
        size_t n = strlen(arg);
        if ( n > 4 && !strcmp(arg+n-4, ".cym") )
            readFile(arg);
        else if ( n > 2 && 0 == readLine(arg) )
            printf("warning: argument '%s' was ignored\n", arg);
    }
    
    double result = alpha + beta - zeta;
    
    writeFile("result.txt", result);
}
