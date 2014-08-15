#define _GNU_SOURCE
#define _FILE_OFFSET_BITS 64
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <string.h>

extern int * taxons;
extern int max_tax;

int * gi_id = NULL;
int max_gi = 0;

void FreeGiId(void) {
    if(gi_id)
        free(gi_id);
}

/*
 *  Finds the biggest GI and taxon ID in the gi_taxid_nucl.dmp file
 */
void GetMaxGiNuclDmp(char * infile) {
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
                if(tempTok > max_gi)
                    max_gi = tempTok;
            }

            if(token == 1) {
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
 *  Parses input line read from gi_taxid_nucl.dmp file
 *  (from: ImportGiNuclDmp function)
 */
void ParseGiNuclLine(char * input) {
    char * ptr = strtok(input, "\t");
    int token = 0;
    int tgi = -1; //gi id (first column
    int tti = -1; //tax is (second column)

    while(ptr != NULL) {
        if(token == 0)
            tgi = atoi(ptr);

        if(token == 1)
            tti = atoi(ptr);

        token++;
        ptr = strtok(NULL, "\t");
    }

    //if both ids are non-negative -> store data
    if(tgi > 0 && tti > 0) {
        if(tgi < max_gi)
            gi_id[tgi] = tti;
    }
}

/*
 *  Imports GI and taxon ID pairs into memory from gi_taxid_nucl.dmp file.
 *  The GI ID is the index, the taxon ID is the value in 'int *gi_id array'
 */
void ImportGiNuclDmp(char * infile) {
    FILE * handle = fopen(infile, "r+");
    char line[BUFSIZ];
    int i = 0;

    if(max_gi <= 0 || max_tax <= 0) {
        printf("No ids detected from %s\n", infile);
        return;
    }


    //allocate and empty gi_id array
    gi_id = (int *)calloc(max_gi + 1, sizeof(int));
    for(i = 0; i < max_gi; i++)
       gi_id[i] = -1;

    while(fgets(line,sizeof(line), handle)) {
        ParseGiNuclLine(line);
    }

    fclose(handle);
}
