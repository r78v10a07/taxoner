
#define _GNU_SOURCE
#define _FILE_OFFSET_BITS 64
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <stdbool.h>
#include "btree.h"
#include "utils.h"

void freeBtreeData(void *d) {
    int *g = (int *) d;
    free(g);
}

node *createGiIndex(char * filename) {
    FILE *fb;
    off_t fileLen;
    off_t pos;
    char *line = NULL;
    size_t len = 0;
    ssize_t read;
    node *root = NULL;
    int gi;
    int *taxid;
    int count = 0;
    float percent = 0;
    float estimated;


    fb = fopen(filename, "r");
    if (!fb) {
        fprintf(stdout, "Unable to open file %s\n", filename);
        exit(-1);
    }

    fseeko(fb, 0, SEEK_END);
    fileLen = ftello(fb);
    fseeko(fb, 0, SEEK_SET);

    while ((read = getline(&line, &len, fb)) != -1) {
        taxid = (int *) malloc(sizeof (int));
        sscanf(line, "%d\t%d\n", &gi, taxid);
        root = insert(root, gi, taxid);
        pos = ftello(fb);
        percent = (pos * 100) / fileLen;
        if (count % 100000 == 0) {
            printf("Reading GIs: Total: %10d\t\tPercent: %6.2f%%\t\t\r", count, percent);
        }
        count++;
    }

    if (line) free(line);
    fclose(fb);
    printf("Reading GIs: Total: %10d\t\tPercent: %6.2f%%\t\t\n", count, percent);
    printf("\n\tThere are %d GIs into the B+Tree.\n\n", count);
    fflush(NULL);
    return root;
}

node *createTaxIndex(char * filename) {
    struct timespec start, stop;
    FILE *fb;
    off_t fileLen;
    off_t pos;
    char *line = NULL;
    size_t len = 0;
    ssize_t read;
    node *root = NULL;
    int *taxid;
    int count = 0;
    float percent = 0;
    float estimated;

    fb = fopen(filename, "r");
    if (!fb) {
        fprintf(stdout, "Unable to open file %s\n", filename);
        exit(-1);
    }

    fseeko(fb, 0, SEEK_END);
    fileLen = ftello(fb);
    fseeko(fb, 0, SEEK_SET);

    while ((read = getline(&line, &len, fb)) != -1) {
        taxid = (int *) malloc(sizeof (int));
        sscanf(line, "%d\n", taxid);
        root = insert(root, *taxid, taxid);
        pos = ftello(fb);
        percent = (pos * 100) / fileLen;
        if (count % 100000 == 0) {
            printf("Reading TaxIds: Total: %10d\t\tPercent: %6.2f%%\t\t\r", count, percent);
        }
        count++;
    }

    if (line) free(line);
    fclose(fb);
    printf("Reading TaxIds: Total: %10d\t\tPercent: %6.2f%%\t\t\n", count, percent);
    printf("\n\tThere are %d TaxIds into the B+Tree\n\n", count);
    fflush(NULL);
    return root;
}