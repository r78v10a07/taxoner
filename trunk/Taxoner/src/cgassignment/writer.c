/* 
 * File:   writer.c
 * Author: roberto
 *
 * Created on February 4, 2014, 3:09 PM
 */
#define _GNU_SOURCE
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdbool.h>
#include "btree.h"
#include "genbank.h"
#include "index.h"
#include "writer.h"

#define _STATIC_SIZE_ 10000

void writer(char *text, char *bin, char *index, int verbose) {
    int i;
    FILE *fp;
    FILE *fb;
    FILE *fo;
    off_t pos = 0;
    int count = 0;
    char *line = NULL;
    size_t len = 0;
    ssize_t read;
    int prev = 0;
    genBank_t *g = NULL;

    int gi, taxId, proteinGi, pFrom, pTo;
    char *cog = malloc(sizeof (char) * _STATIC_SIZE_);
    char *protClust = malloc(sizeof (char) * _STATIC_SIZE_);
    char *prot = malloc(sizeof (char) * _STATIC_SIZE_);
    char *locus = malloc(sizeof (char) * _STATIC_SIZE_);
    char *glocus = malloc(sizeof (char) * _STATIC_SIZE_);

    fp = fopen(text, "r");
    if (!fp) {
        printf("Unable to open file!");
        return;
    }

    fb = fopen(bin, "wb");
    if (!fb) {
        printf("Unable to open file!");
        return;
    }

    fo = fopen(index, "wb");
    if (!fo) {
        printf("Unable to open file!");
        return;
    }

    while ((read = getline(&line, &len, fp)) != -1) {
        memset(cog, 0, sizeof (char) * _STATIC_SIZE_);
        memset(protClust, 0, sizeof (char) * _STATIC_SIZE_);
        memset(prot, 0, sizeof (char) * _STATIC_SIZE_);
        memset(locus, 0, sizeof (char) * _STATIC_SIZE_);
        sscanf(line, "%d\t%s\t%d\t%s\t%s\t%d\t%d\t%s\t%s", &gi, glocus, &taxId, prot, locus, &pFrom, &pTo, cog, protClust);

        if (prev != gi) {
            if (g != NULL) {
                count++;
                if (verbose == 1) printGenBank(g);
                fwrite(&(g->gi), sizeof (int), 1, fo);
                fwrite(&(pos), sizeof (off_t), 1, fo);
                printGenBankBinary(fb, g);
                freeGenBank(g);
                pos = ftello(fb);
            }
            g = initGenBank(gi, glocus, taxId, prot, locus, pFrom, pTo, cog, protClust);
            prev = gi;
        } else {
            addCDS(&(g->cds), &(g->cds_number), prot, locus, pFrom, pTo, cog, protClust);
        }
    }
    count++;
    if (verbose == 1) printGenBank(g);
    pos = ftello(fb);
    fwrite(&(g->gi), sizeof (int), 1, fo);
    fwrite(&(pos), sizeof (long int), 1, fo);
    printGenBankBinary(fb, g);

    printf("There are %d GI in the index file\n", count);

    fclose(fb);
    fclose(fo);
    fclose(fp);
    freeGenBank(g);
    if (line) free(line);
    if (cog) free(cog);
    if (protClust) free(protClust);
    if (locus) free(locus);
    if (glocus) free(glocus);
    if (prot) free(prot);
}