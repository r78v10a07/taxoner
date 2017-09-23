/*
 * File:   main.c
 * Author: lorinc and roberto
 *
 * Created on August 15, 2014, 11:27 AM
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <getopt.h>
#include <time.h>
#include <stdbool.h>
#include <pthread.h>
#include <unistd.h>
#include <ctype.h>

#include "inputs.h"
#include "main.h"
#include "nodes.h"

char * program_name;
struct params * parms;
struct database * db;
struct TaxonNodes * nodes;

/*
 *  Main function
 */
int main(int argc, char** argv)
{
    struct timespec start, stop;

    clock_gettime(CLOCK_MONOTONIC, &start);
    program_name = argv[0];

    readInputs(argc, argv);
    ImportNodes();
    if(parms->neighborOnly == 0) AlignReads();
    GetNeighbors();
    FreeInputs();
    clock_gettime(CLOCK_MONOTONIC, &stop);
    printf("\n\tThe total time was %lu sec\n\n", timespecDiff(&stop, &start) / 1000000000);

    return 0;
}
