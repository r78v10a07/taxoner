/* 
 * File:   index.c
 * Author: roberto
 *
 * Created on February 6, 2014, 10:20 AM
 */

#define _GNU_SOURCE
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <stdbool.h>
#include "index.h"
#include "btree.h"
#include "genbank.h"

int cmpfun(void *d1, void *d2) {
    giOffset_t *g1 = (giOffset_t *) d1;
    giOffset_t *g2 = (giOffset_t *) d2;

    return g1->gi - g2->gi;
}

void freeData(void *d) {
    off_t *g = (off_t *) d;
    free(g);
}

void freeDataAssign(void *d) {
    int *o = (int *) d;
    free(o);
}

void printData(void *d) {
    giOffset_t *g = (giOffset_t *) d;
    printf("%d", g->gi);
}

giOffset_t *readIndex(int *count, char *index) {
    FILE *fi;
    off_t fileLen;
    off_t pos;
    size_t len;
    giOffset_t *list = NULL;

    fi = fopen(index, "rb");
    if (!fi) {
        printf("Unable to open file!");
        return;
    }

    //Get file length
    fseeko(fi, 0, SEEK_END);
    fileLen = ftello(fi);
    fseeko(fi, 0, SEEK_SET);

    pos = 0;
    len = 0;
    while (pos < fileLen) {
        list = realloc(list, sizeof (giOffset_t) * (len + 1));
        fread(&(list[len].gi), sizeof (int), 1, fi);
        fread(&(list[len].offset), sizeof (long int), 1, fi);
        pos = ftell(fi);
        len++;
    }

    printf("There are %d gi\n", len);

    *count = len;
    return list;
}

node *readIndexBtree(char *index) {
    struct timespec start, stop;
    node *root = NULL;
    FILE *fi;
    off_t fileLen;
    off_t pos;
    off_t *value;
    int gi;
    int count = 0;

    clock_gettime(CLOCK_MONOTONIC, &start);
    fi = fopen(index, "rb");
    if (!fi) {
        printf("Unable to open file!");
        return;
    }

    fseeko(fi, 0, SEEK_END);
    fileLen = ftello(fi);
    fseeko(fi, 0, SEEK_SET);

    pos = 0;
    while (pos < fileLen) {
        value = malloc(sizeof (off_t));
        fread(&gi, sizeof (int), 1, fi);
        fread(value, sizeof (off_t), 1, fi);

        root = insert(root, gi, value);
        pos = ftell(fi);
        count++;
    }
    fclose(fi);
    clock_gettime(CLOCK_MONOTONIC, &stop);
    printf("\n\tThere are %d GIs into the B+Tree. Elapsed time: %lu sec\n\n", count, timespecDiff(&stop, &start) / 1000000000);
    fflush(NULL);
    return root;
}

int64_t timespecDiff(struct timespec *timeA_p, struct timespec *timeB_p) {
    return ((timeA_p->tv_sec * 1000000000) + timeA_p->tv_nsec) -
            ((timeB_p->tv_sec * 1000000000) + timeB_p->tv_nsec);
}

void createBTreeIndex(char *text, char *bin, char *index, char *output, int verbose) {
    int i, j, k;
    struct timespec start, stop;
    node *root;
    node *assign = NULL;
    FILE *fb;
    FILE *ft;
    FILE *fo;
    record *rec;
    genBank_t *g = NULL;
    int count = 0;
    char *line = NULL;
    size_t len = 0;
    ssize_t read;
    off_t oTotal, oCurr;
    int taxId;
    int gi;
    float score;
    int pFrom;
    int pTo;
    int noGene = 0;

    fb = fopen(bin, "rb");
    if (!fb) {
        fprintf(stdout, "Unable to open file %s\n", bin);
        exit(-1);
    }

    ft = fopen(text, "r");
    if (!ft) {
        fprintf(stdout, "Unable to open file %s\n", text);
        exit(-1);
    }

    fo = fopen(output, "w");
    if (!fo) {
        fprintf(stdout, "Unable to open file %s\n", output);
        exit(-1);
    }

    root = readIndexBtree(index);

    fseeko(ft, 0, SEEK_END);
    oTotal = ftello(ft);
    fseeko(ft, 0, SEEK_SET);

    k = 1;
    clock_gettime(CLOCK_MONOTONIC, &start);
    while ((read = getline(&line, &len, ft)) != -1) {
        oCurr = ftello(ft);
        if ((i = sscanf(line, "%*s\t%d\t%d\t%f\t%d\t%d", &taxId, &gi, &score, &pFrom, &pTo)) == 5) {
            clock_gettime(CLOCK_MONOTONIC, &stop);
            j = timespecDiff(&stop, &start) / 1000000000;
            if (verbose) {
                printf("\tLines reads %d. Reads without genes %d. GIs with genes %d. Elapsed time %lu sec. Estimated time %lu sec.\r", k, noGene, count, j, (j * oTotal / oCurr));
            }
            rec = find(root, gi, false);
            if (rec != NULL) {
                if (!assignGenes(&assign, &g, &count, fb, *((off_t *) rec->value), pFrom, pTo)) {
                    noGene++;
                }
            } else {
                noGene++;
            }
            k++;
        } else {
            printf("Can't parse the line: %d\t", i);
            printf("[[%s]]\n", line);
            exit(-1);
        }
    }
    printf("\tLines reads %d. Reads without genes %d. GIs with genes %d. Elapsed time %lu sec. Estimated time %lu sec.\n", k++, noGene, count, j, (j * oTotal / oCurr));
    fflush(NULL);
    if (line) free(line);
    line = malloc(sizeof (char) * 1);
    len = 1;

    if (g != NULL) {
        qsort(g, count, sizeof (genBank_t), cmpGenBank);
        for (i = 0; i < count; i++) {
            if (g[i].cds != NULL) {
                qsort(g[i].cds, g[i].cds_number, sizeof (cds_t), cmpCDSSum);
                for (j = 0; j < g[i].cds_number; j++) {

                    line[0] = '\0';
                    if (g[i].cds[j].cog_number != 0) {
                        for (k = 0; k < g[i].cds[j].cog_number; k++) {
                            if (len <= strlen(line) + strlen(g[i].cds[j].cog[k]) + 2) {
                                len += strlen(g[i].cds[j].cog[k]) + 2;
                                if (k < g[i].cds[j].cog_number - 1) len++;
                                line = (char *) realloc(line, sizeof (char) * len);
                            }
                            strcat(line, g[i].cds[j].cog[k]);
                            if (k < g[i].cds[j].cog_number - 1) {
                                strcat(line, ",");
                            }
                        }
                        if (strlen(line) + 3 < len) {
                            len += 3;
                            line = (char *) realloc(line, sizeof (char) * len);
                        }
                        strcat(line, "\t");
                    } else {
                        if (len < 3) {
                            len = 3;
                            line = (char *) realloc(line, sizeof (char) * len);
                        }
                        strcat(line, "-\t");
                    }
                    if (g[i].cds[j].protClust_number != 0) {
                        for (k = 0; k < g[i].cds[j].protClust_number; k++) {
                            if (len <= strlen(line) + strlen(g[i].cds[j].protClust[k]) + 2) {
                                len += strlen(g[i].cds[j].protClust[k]) + 2;
                                if (k < g[i].cds[j].protClust_number - 1) len++;
                                line = (char *) realloc(line, sizeof (char) * len);
                            }
                            strcat(line, g[i].cds[j].protClust[k]);
                            if (k < g[i].cds[j].protClust_number - 1) {
                                strcat(line, ",");
                            }
                        }
                    } else {
                        if (strlen(line) + 4 > len) {
                            len += 4;
                            line = (char *) realloc(line, sizeof (char) * len);
                        }
                        strcat(line, "-");
                    }

                    fprintf(fo, "%s\t%d\t%s\t%s\t%d\t%s\n",
                            g[i].locusName, g[i].taxId,
                            g[i].cds[j].proteinId,
                            g[i].cds[j].locusName,
                            g[i].cds[j].hits,
                            line);
                    freeCDS(g[i].cds + j);
                }
                free(g[i].cds);
            }
            if (g[i].locusName) free(g[i].locusName);
        }
        free(g);
    }

    destroy_tree(assign, freeDataAssign);
    destroy_tree(root, freeData);
    if (line) free(line);
    fclose(fb);
    fclose(ft);
    fclose(fo);
}