
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
    struct timespec start, stop;
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


    clock_gettime(CLOCK_MONOTONIC, &start);
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
        if (count % 1000 == 0) {
            clock_gettime(CLOCK_MONOTONIC, &stop);
            estimated = (fileLen * (timespecDiff(&stop, &start) / 1000000000)) / pos;
            printf("Reading GIs: Total: %10d\t\tPercent: %6.2f%%\t\tEstimated time: %10.1f s   \r", count, percent, estimated);
        }
        count++;
    }

    if (line) free(line);
    fclose(fb);
    clock_gettime(CLOCK_MONOTONIC, &stop);
    printf("Reading GIs: Total: %10d\t\tPercent: %6.2f%%\t\tEstimated time: %10.1f s   \r", count, percent, estimated);
    printf("\n\tThere are %d GIs into the B+Tree. Elapsed time: %lu sec\n\n", count, timespecDiff(&stop, &start) / 1000000000);
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


    clock_gettime(CLOCK_MONOTONIC, &start);
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
        if (count % 10000 == 0) {
            clock_gettime(CLOCK_MONOTONIC, &stop);
            estimated = (fileLen * (timespecDiff(&stop, &start) / 1000000000)) / pos;
            printf("Reading TaxIds: Total: %10d\t\tPercent: %6.2f%%\t\tEstimated time: %10.1f s   \r", count, percent,estimated);
        }
        count++;
    }

    if (line) free(line);
    fclose(fb);
    clock_gettime(CLOCK_MONOTONIC, &stop);
    printf("Reading TaxIds: Total: %10d\t\tPercent: %6.2f%%\t\tEstimated time: %10.1f s   \r", count, percent, estimated);
    printf("\n\tThere are %d TaxIds into the B+Tree. Elapsed time: %lu sec\n\n", count, timespecDiff(&stop, &start) / 1000000000);
    fflush(NULL);
    return root;
}