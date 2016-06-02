/* 
 * File:   main.c
 * Author: roberto and lorinc
 *
 * Created on April 5, 2014, 11:27 AM
 */

#define _GNU_SOURCE
#define _FILE_OFFSET_BITS 64
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <getopt.h>
#include <stdbool.h>
#include <time.h>

#include "utils.h"
#include "files.h"
#include "btree.h"
#include "index.h"
#include "fasta.h"

char *program_name;

/*
 *  Main function
 */
int main(int argc, char** argv) {

    int next_option, verbose, write;
    const char* const short_options = "hvn:g:o:s:i:";
    char *nt, *gi, *nodes, *skip, *include;
    
    node *giIndex = NULL;
    node *includeIndex = NULL;
    node *skipIndex = NULL;

    program_name = argv[0];

    const struct option long_options[] = {
        { "help", 0, NULL, 'h'},
        { "verbose", 0, NULL, 'v'},
        { "nt", 1, NULL, 'n'},
        { "gi", 1, NULL, 'g'},
        { "nodes", 1, NULL, 'o'},
        { "skip", 1, NULL, 's'},
        { "include", 1, NULL, 'i'},
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
        }
    } while (next_option != -1);

    InputCheck(nt, gi, nodes, skip, include);
    
    printf("\nReading GIs and taxIds for Btree creation\n");
    giIndex = createGiIndex(gi);
    
    if (include != NULL){
        printf("\nReading TaxIds for the include list for Btree creation\n");
        includeIndex = createTaxIndex(include);
    }
    
    if (skip != NULL){
        printf("\nReading TaxIds for the skip list for Btree creation\n");
        skipIndex = createTaxIndex(skip);
    }
    
    printf("\nParsing fasta file\n");
    ReadFasta(nt, giIndex, includeIndex, skipIndex);
    
    destroy_tree(giIndex, freeBtreeData);
    destroy_tree(includeIndex, freeBtreeData);
    destroy_tree(skipIndex, freeBtreeData);

    if (nt) free(nt);
    if (gi) free(gi);
    if (nodes) free(nodes);
    if (skip) free(skip);
    if (include) free(include);
    return (EXIT_SUCCESS);
}
