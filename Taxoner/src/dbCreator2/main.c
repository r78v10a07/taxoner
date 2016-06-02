/* 
 * File:   main.c
 * Author: roberto and lorinc
 *
 * Created on August 15, 2014, 11:27 AM
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <getopt.h>
#include <time.h>
#include <stdbool.h>

#include "utils.h"
#include "importGI.h"
#include "importNodes.h"
#include "IncludeExclude.h"
#include "parsefasta.h"
char *program_name;
char *nt, *gi, *nodes, *skip, *include;
extern int MaxGb;
/*
 *  Main function
 */
int main(int argc, char** argv)
{
    int next_option, verbose, write;
    const char* const short_options = "hvn:g:o:s:i:";

    char *giIndex = NULL;
    char *includeIndex = NULL;
    char *skipIndex = NULL;

    program_name = argv[0];

    const struct option long_options[] = {
        { "help", 0, NULL, 'h'},
        { "verbose", 0, NULL, 'v'},
        { "nt", 1, NULL, 'n'},
        { "gi", 1, NULL, 'g'},
        { "nodes", 1, NULL, 'o'},
        { "skip", 1, NULL, 's'},
        { "include", 1, NULL, 'i'},
        { "dbSize", 1, NULL, 'd'},
        { NULL, 0, NULL, 0} /* Required at end of array.  */
    };

    write = verbose = 0;
    nt = gi = nodes = skip = include = NULL;
    do {
        next_option = getopt_long(argc, argv, short_options, long_options, NULL);

        switch (next_option) {
            case 'h':
                print_usage(stdout, 0);

            case 'v':
                verbose = 1;
                break;

            case 'n':
                nt = (char *) malloc(sizeof (char) * (strlen(optarg) + 1));
                strcpy(nt, optarg);
                break;

            case 'g':
                gi = (char *) malloc(sizeof (char) * (strlen(optarg) + 1));
                strcpy(gi, optarg);
                break;

            case 'o':
                nodes = (char *) malloc(sizeof (char) * (strlen(optarg) + 1));
                strcpy(nodes, optarg);
                break;

            case 's':
                skip = (char *) malloc(sizeof (char) * (strlen(optarg) + 1));
                strcpy(skip, optarg);
                break;

            case 'i':
                include = (char *) malloc(sizeof (char) * (strlen(optarg) + 1));
                strcpy(include, optarg);
                break;

            case 'd':
                MaxGb = atoi(optarg);
                break;
        }
    } while (next_option != -1);

    InputCheck(nt, gi, nodes, skip, include);

    if(verbose == 1)
        print_parameters();

    if(verbose == 1)
        printf("Getting biggest GI id from %s\n", gi);
    GetMaxGiNuclDmp(gi);
    if(verbose == 1)
        printf("Importing GI and taxon ID pairs into memory from %s\n", gi);
    ImportGiNuclDmp(gi);

    if(verbose == 1)
        printf("Getting biggest taxon ID from %s\n", nodes);
    GetMaxNodeslDmp(nodes);
    if(verbose == 1)
        printf("Importing taxon parent and child IDs into memory from %s\n", nodes);
    ImportNodeslDmp(nodes);

    if(include != NULL) {
        if(verbose == 1)
            printf("Importing include taxon IDs from %s\n", include);
        ImportInclude(include);
    }

    if(skip != NULL) {
        if(verbose == 1)
            printf("Importing exclude taxon IDs from %s\n", skip);
        ImportExclude(skip);
    }
    if(verbose == 1)
        printf("Reading and parsing fasta file %s\n", nt);
    ReadFasta(nt);

    if (nt) free(nt);
    if (gi) free(gi);
    if (nodes) free(nodes);
    if (skip) free(skip);
    if (include) free(include);
    FreeIncludeExclude();
    FreeGiId();
    FreeTaxons();

    return 0;
}
