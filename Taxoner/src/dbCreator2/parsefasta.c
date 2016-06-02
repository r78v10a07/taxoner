#define _GNU_SOURCE
#define _FILE_OFFSET_BITS 64
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <string.h>

#include "parsefasta.h"
#include "IncludeExclude.h"

extern int max_gi;
extern int *gi_id;

extern int max_tax;
extern int *taxons;

extern char *skip, *include, *nt;

int MbSize = 1048576;
int GbSize = 1024;

int MaxGb = 4;

int tempChars;
int tempMb;
int tempGb;

int FastaNum = 0;

FILE * OutFasta = NULL;

int printData = 0; //check to print fasta file

/*
 *  Open new fasta file to print results from NT file
 */
void OpenFastaFile(void) {
    char buffer[50];

    if (OutFasta != NULL)
        fclose(OutFasta);

    sprintf(buffer, "%d.fasta", FastaNum);
    printf("Index: %s\n", buffer);
    OutFasta = fopen(buffer, "w+");

    FastaNum++;
}

/*
 *  Count the current amaount of chars printed to last fasta file
 */
void CountData(int a) {
    int i;

    if (tempChars >= MbSize) {
        while (tempChars > MbSize) {
            tempChars -= MbSize;
            tempMb++;
        }

        if (tempMb >= GbSize) {
            while (tempMb >= GbSize) {
                tempGb++;
                tempMb -= GbSize;
            }

            if (tempGb >= MaxGb) {
                OpenFastaFile();
                tempGb = 0;
                tempMb = 0;
                tempChars = 0;
            }
        }
    }
}

/*
 *  Get taxon id from GI id
 */
int FindTaxonWithGi(int gids) {
    if (gids > max_gi || gids <= 0) {
        return -1;
    }

    return gi_id[gids];
}

/*
 *  Compare strings by length and chars
 */
int compareStrings(char * a, char * b) {
    if (strlen(a) == strlen(b)) {
        if (strncmp(a, b, strlen(a)) == 0) {
            return 1;
        }
    }

    return -1;
}

/*
 *  Parse NT fasta file and get GI id
 */
void ParseFastaTitle(char * input) {
    char * ptr = strtok(input, "|");
    int token = 0;
    int taxa = -1;
    char buffer[BUFSIZ];

    while (ptr != NULL) {
        if (token == 1) {
            token++;
            taxa = FindTaxonWithGi(atoi(ptr));
            if (taxa != -1 && taxa > 0 && taxa < max_gi && CheckInclude(taxa, -1) == 1 && CheckExclude(taxa, 1) == 1) {
                sprintf(buffer, ">%d;%d\n", atoi(ptr), taxa);
                printData = 1;
                CountData(strlen(buffer));
                fprintf(OutFasta, "%s", buffer);
            }
        }

        if (compareStrings("gi", ptr) == 1 || compareStrings(">gi", ptr) == 1)
            token++;

        ptr = strtok(NULL, "|");
    }
}

/*
 *  Read NT fasta file
 */
void ReadFasta(char * inFile) {
    FILE * handle;

    if (strcmp(inFile, "-") != 0) {
        handle = fopen(inFile, "r+");
    } else {
        handle = stdin;
    }

    char line[1000000]; //large buffer for very long fasta names

    OpenFastaFile();

    while (fgets(line, sizeof (line), handle)) {
        if (line[0] == '>') {
            printData = 0;
            ParseFastaTitle(line);
        } else {
            if (printData == 1) {
                tempChars += strlen(line);
                fprintf(OutFasta, "%s", line);
            }
        }
    }

    if (strcmp(inFile, "-") != 0) {
        fclose(handle);
    }
}
