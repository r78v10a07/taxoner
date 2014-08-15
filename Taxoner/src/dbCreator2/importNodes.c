#define _GNU_SOURCE
#define _FILE_OFFSET_BITS 64
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <string.h>
#include "importGI.h"

int *taxons;
int max_tax = 0;

extern int max_gi;
extern int *gi_id;

void FreeTaxons(void) {
    if(taxons)
        free(taxons);
}

/*
 *  Finds the biggest taxon ID in the nodes.dmp file
 */
void GetMaxNodeslDmp(char * infile) {
    FILE * handle = fopen(infile, "r+");
    char line[BUFSIZ];
    char * ptr;
    int token = 0;
    int tempTok;

    while(fgets(line,sizeof(line), handle)) {
        ptr = strtok(line, "\t");
        token = 0;

        while(ptr != NULL) {
            if(token == 0) {
                tempTok = atoi(ptr);
                if(max_tax > max_gi)
                    max_tax = tempTok;
            }

            if(token == 2) {
                tempTok = atoi(ptr);
                if(tempTok > max_tax)
                    max_tax = tempTok;
            }

            token++;
            ptr = strtok(NULL, "\t");
        }
    }

    fclose(handle);
}

/*
 *  Imports taxon ID and parent taxon ID into memory from nodes.dmp file.
 *  The child taxon ID is the index, the parent taxon ID is the value in 'int *taxon array'
 */
void ImportNodeslDmp(char * infile) {
    FILE * handle = fopen(infile, "r+");
    char line[BUFSIZ];
    char * ptr;
    int token = 0;
    int tempTok;
    int i;
    int parent;
    int child;

    taxons = (int *)calloc(max_tax + 1, sizeof(int));
    for(i = 0; i < max_tax; i++)
        taxons[i] = -1;

    while(fgets(line,sizeof(line), handle)) {
        ptr = strtok(line, "\t");
        token = 0;
        parent = -1;
        child = -1;

        while(ptr != NULL) {
            if(token == 0) {
                tempTok = atoi(ptr);
                if(max_tax > tempTok)
                    parent = tempTok;
            }

            if(token == 2) {
                tempTok = atoi(ptr);
                if(tempTok < max_tax)
                    child = tempTok;
            }

            token++;
            ptr = strtok(NULL, "\t");
        }

        if(parent != -1 && child != -1)
            taxons[parent] = child;
    }

    fclose(handle);
}
