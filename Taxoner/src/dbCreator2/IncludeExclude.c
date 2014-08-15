#define _GNU_SOURCE
#define _FILE_OFFSET_BITS 64
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <string.h>
#include "importGI.h"

extern int max_gi;
extern int *gi_id;

extern int * taxons;
extern int max_tax;

extern char *skip, *include;

int *includeTaxon = NULL;
int *excludeTaxon = NULL;

void FreeIncludeExclude(void) {
    if(includeTaxon)
        free(includeTaxon);

    if(excludeTaxon)
        free(excludeTaxon);
}

/*
 *  Check if taxon ID or its lineage is in the include list
 */
int CheckInclude(int source, int status) {
    if(includeTaxon == NULL)
        return 1;

    if(source < max_tax) {
        if(taxons[source] != -1) {
            if(includeTaxon[source] == 1) {
                status = 1;
                return status;
            }

            if(taxons[source] == 1)
                return status;

            else
                return CheckInclude(taxons[source], status);
        }
    }

    return -1;
}

/*
 *  Check if taxon ID or its lineage is in the exclude list
 */
int CheckExclude(int source, int status) {
    if(skip == NULL)
        return 1;

    if(source < max_tax) {
        if(taxons[source] != -1) {
            if(excludeTaxon[source] == 1) {
                status = -1;
                return status;
            }

            if(taxons[source] == 1)
                return status;

            else
                return CheckExclude(taxons[source], status);
        }
    }

    return -1;
}

/*
 *  Import taxons that have to be included from the include file
 */
void ImportInclude(char * infile) {
    FILE * handle = fopen(infile, "r+");
    char line[BUFSIZ];
    int token = 0;
    int tempTok;

    includeTaxon = (int *)calloc(max_tax + 1, sizeof(int));

    for(tempTok = 0; tempTok < max_tax; tempTok++)
        includeTaxon[tempTok] = -1;

    while(fgets(line,sizeof(line), handle)) {
        tempTok = atoi(line);

        if(tempTok > max_tax || tempTok <= 0) {
            printf("Taxon \"%d\" not found (larger than max taxon in nodes.dmp)\n", tempTok);
        }

        else {
            if(taxons[tempTok] == -1) {
                printf("Taxon \"%d\" has no parent taxon in lineage\n", tempTok);
            }

            else {
                includeTaxon[tempTok] = 1;
            }
        }
    }

    fclose(handle);
}

/*
 *  Import taxons that have to be excluded from the exclude file
 */
void ImportExclude(char * infile) {
    FILE * handle = fopen(infile, "r+");
    char line[BUFSIZ];
    int token = 0;
    int tempTok;

    excludeTaxon = (int *)calloc(max_tax + 1, sizeof(int));

    for(tempTok = 0; tempTok < max_tax; tempTok++)
        excludeTaxon[tempTok] = -1;

    while(fgets(line,sizeof(line), handle)) {
        tempTok = atoi(line);

        if(tempTok > max_tax || tempTok <= 0) {
            printf("Taxon \"%d\" not found (larger than max taxon in nodes.dmp)\n", tempTok);
        }

        else {
            if(taxons[tempTok] == -1) {
                printf("Taxone \"%d\" has no parent taxon in lineage\n", tempTok);
            }

            else {
                excludeTaxon[tempTok] = 1;
            }
        }
    }

    fclose(handle);
}
