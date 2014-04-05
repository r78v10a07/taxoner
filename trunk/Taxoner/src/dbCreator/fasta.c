

#define _GNU_SOURCE
#define _FILE_OFFSET_BITS 64
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <stdbool.h>
#include "btree.h"
#include "utils.h"

FILE * OutFasta = NULL;

int maxID;
int MbSize = 1048576;
int GbSize = 1024;

int MaxGb = 4;

int tempChars;
int tempMb;
int tempGb;

int FastaNum = 0;

void EmptyFileData(void) {
    tempMb = 0;
    tempGb = 0;
    tempChars = 0;
}

void OpenFastaFile(void) {
    char buffer[50];

    if (OutFasta != NULL)
        fclose(OutFasta);

    sprintf(buffer, "%d.fasta", FastaNum);
    printf("\nIndex: %s\n", buffer);
    OutFasta = fopen(buffer, "w+");

    FastaNum++;
}

void CountData(int a) {
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

void PrintTitle(int a, int b) {
    char buffer[BUFSIZ];

    sprintf(buffer, ">%d;%d", a, b);
    CountData(strlen(buffer));

    fprintf(OutFasta, "%s\n", buffer);
}

int checkIfBlacklisted(int taxId, node *skipIndex) {
    if (skipIndex != NULL) {
        /**
         * Not NULL if the TaxId is include into the include list
         */
        record *rec = find(skipIndex, taxId, false);

        if (rec != NULL) {
            return 1;
        }
    }
    return 0;
}

int searchId(int gi, node *giIndex, node *includeIndex) {
    int taxon;

    /**
     * Get the TaxId from the Gi index
     */
    record *rec = find(giIndex, gi, false);

    if (rec == NULL)
        return -1;

    taxon = *((int *) rec->value);

    if (includeIndex != NULL) {
        /**
         * Not NULL if the TaxId is include into the include list
         */
        rec = find(includeIndex, taxon, false);

        if (rec == NULL)
            return -1;
    }

    if (taxon > 0)
        return taxon;

    return -1;
}

int FormatFastaTitle(char * title, node *giIndex, node *includeIndex, node *skipIndex) {
    char * ptr = strtok(title, "|");
    int nTok = 0;
    char * temp = NULL;
    record *rec;

    while (ptr != NULL) {
        if (temp != NULL)
            free(temp);

        if (ptr[0] == '>') {
            ptr += 1;
            temp = CopyString(ptr);
            ptr -= 1;
        } else
            temp = CopyString(ptr);

        if (CompareStrings(temp, "gi") == 0)
            nTok = 1;

        if (nTok == 2) {
            /*
             * If the index of the include exist then it check if the taxId for 
             * that GI is included.
             * If the index of the include does not exist the just return the 
             * TaxId
             */
            nTok = searchId(atoi(ptr), giIndex, includeIndex);

            if (nTok == -1) {
                free(temp);
                return 0;
            } else {
                free(temp);
                /*
                 * Check if the TaxId is included into the blacklist
                 * 
                 * Cero if it not included
                 * 
                 */
                if (checkIfBlacklisted(nTok, skipIndex) == 0) {
                    PrintTitle(atoi(ptr), nTok);
                    return 1;
                }

                return 0;
            }
        }

        if (nTok == 1)
            nTok++;

        ptr = strtok(NULL, "|");
    }

    return 0;
}

void ReadFasta(char * filename, node *giIndex, node *includeIndex, node *skipIndex) {
    struct timespec start, stop;
    FILE *fb;
    off_t fileLen;
    off_t pos;
    char *line = NULL;
    size_t len = 0;
    ssize_t read;
    int count = 0;
    float percent = 0;
    float estimated;
    int SeqPrint = 0;

    clock_gettime(CLOCK_MONOTONIC, &start);
    fb = fopen(filename, "r");
    if (!fb) {
        fprintf(stdout, "Unable to open file %s\n", filename);
        exit(-1);
    }

    fseeko(fb, 0, SEEK_END);
    fileLen = ftello(fb);
    fseeko(fb, 0, SEEK_SET);
    
    EmptyFileData();
    OpenFastaFile();

    while ((read = getline(&line, &len, fb)) != -1) {
        if (line[0] == '>') {
            count++;
            SeqPrint = FormatFastaTitle(line, giIndex, includeIndex, skipIndex);
        } else {
            if (SeqPrint != 0) {
                fprintf(OutFasta, "%s", line);
                tempChars += strlen(line);
            }
        }
        pos = ftello(fb);
        percent = (pos * 100) / fileLen;
        if (count % 1000 == 0) {
            clock_gettime(CLOCK_MONOTONIC, &stop);
            estimated = (fileLen * (timespecDiff(&stop, &start) / 1000000000)) / pos;
            printf("Reading Fasta entries: Total: %10d\t\tPercent: %6.2f%%\t\tEstimated time: %10.1f s   \r", count, percent, estimated);
        }
    }

    if (line) free(line);
    fclose(fb);
    clock_gettime(CLOCK_MONOTONIC, &stop);
    printf("Reading  Fasta entries: Total: %10d\t\tPercent: %6.2f%%\t\tEstimated time: %10.1f s   \r", count, percent, estimated);
    fflush(NULL);
}