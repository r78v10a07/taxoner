
#define _GNU_SOURCE
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdbool.h>
#include "btree.h"
#include "genbank.h"

cds_t *allocCDS() {
    cds_t *cds = (cds_t *) malloc(sizeof (cds_t));
    cds->proteinId = NULL;
    cds->locusName = NULL;
    cds->cog = NULL;
    cds->cog_number = 0;
    cds->protClust = NULL;
    cds->protClust_number = 0;
    cds->pFrom = 0;
    cds->pTo = 0;
    cds->hits = 0;
    return cds;
}

void freeCDS(cds_t *c) {
    int i;
    if (c != NULL) {
        if (c->proteinId) free(c->proteinId);
        if (c->locusName) free(c->locusName);
        if (c->cog != NULL) {
            for (i = 0; i < c->cog_number; i++) {
                free(c->cog[i]);
            }
            free(c->cog);
        }
        if (c->protClust != NULL) {
            for (i = 0; i < c->protClust_number; i++) {
                free(c->protClust[i]);
            }
            free(c->protClust);
        }
    }
}

int str_split(char*** result, char* a_str, const char a_delim) {
    size_t count = 0;
    char delim[2];
    delim[0] = a_delim;
    delim[1] = 0;

    *result = NULL;
    if (a_str == NULL || a_str[0] == '\0' || a_str[0] == '-') {
        return 0;
    }

    char* token = strtok(a_str, delim);
    while (token) {
        *result = realloc(*result, sizeof (char*) * (count + 1));
        *(*result + count++) = strdup(token);
        token = strtok(0, delim);
    }
    return count;
}

int cmpCDSSort(const void *a, const void * b) {
    cds_t *c1 = (cds_t *) a;
    cds_t *c2 = (cds_t *) b;
    int comp = c1->pFrom - c2->pFrom;
    if (comp == 0) {
        return (c1->pTo - c2->pTo);
    }
    return comp;
}

int cmpCDS(const void *a, const void * b) {
    cds_t *c1 = (cds_t *) a;
    cds_t *c2 = (cds_t *) b;
    int comp = strcmp(c1->proteinId, c2->proteinId);
    if (comp == 0) {
        comp = strcmp(c1->locusName, c2->locusName);
        if (comp == 0) {
            return cmpCDSSort(a, b);
        }
    }

    return comp;
}

int cmpCDSSum(const void *a, const void * b) {
    cds_t *c1 = (cds_t *) a;
    cds_t *c2 = (cds_t *) b;
    int comp = c2->hits - c1->hits;
    if (comp == 0) {
        return strcmp(c1->proteinId, c2->proteinId);
    }
    return comp;
}

int cmpGenBank(const void *a, const void * b) {
    genBank_t *g1 = (genBank_t *) a;
    genBank_t *g2 = (genBank_t *) b;
    int comp = g1->taxId - g2->taxId;
    if (comp == 0) {
        return strcmp(g1->locusName, g2->locusName);
    }
    return comp;
}

cds_t *initCDS(char *proteinId, char *locusName, int pFrom, int pTo, char *cog, char *protClust) {
    cds_t *cds = allocCDS();
    cds->proteinId = strdup(proteinId);
    cds->locusName = strdup(locusName);
    cds->cog_number = str_split(&(cds->cog), cog, ',');
    cds->protClust_number = str_split(&(cds->protClust), protClust, ',');
    cds->pFrom = pFrom;
    cds->pTo = pTo;
    return cds;
}

void cdsdup(cds_t **des, int *index, cds_t src) {
    int i;

    *des = (cds_t *) realloc((*des), sizeof (cds_t) * (*index + 1));

    (*des)[*index].proteinId = strdup(src.proteinId);
    (*des)[*index].locusName = strdup(src.locusName);
    (*des)[*index].cog_number = src.cog_number;
    (*des)[*index].protClust_number = src.protClust_number;
    (*des)[*index].pFrom = src.pFrom;
    (*des)[*index].pTo = src.pTo;
    (*des)[*index].hits = src.pTo;

    if ((*des)[*index].cog_number == 0) {
        (*des)[*index].cog = NULL;
    } else {
        (*des)[*index].cog = (char **) malloc(sizeof (char **) * (*des)[*index].cog_number);
        for (i = 0; i < (*des)[*index].cog_number; i++) {
            (*des)[*index].cog[i] = strdup(src.cog[i]);
        }
    }

    if ((*des)[*index].protClust_number == 0) {
        (*des)[*index].protClust = NULL;
    } else {
        (*des)[*index].protClust = (char **) malloc(sizeof (char **) * (*des)[*index].protClust_number);
        for (i = 0; i < (*des)[*index].protClust_number; i++) {
            (*des)[*index].protClust[i] = strdup(src.protClust[i]);
        }
    }

    (*index)++;
}

void addCDS(cds_t **cds, int *index, char *proteinId, char *locusName, int pFrom, int pTo, char *cog, char *protClust) {
    *cds = (cds_t *) realloc(*cds, sizeof (cds_t) * (*index + 1));
    (*cds)[*index].proteinId = strdup(proteinId);
    (*cds)[*index].locusName = strdup(locusName);
    (*cds)[*index].cog_number = str_split(&((*cds)[*index].cog), cog, ',');
    (*cds)[*index].protClust_number = str_split(&((*cds)[*index].protClust), protClust, ',');
    (*cds)[*index].pFrom = pFrom;
    (*cds)[*index].pTo = pTo;
    (*cds)[*index].hits = 0;
    (*index)++;
}

genBank_t *allocGenBank() {
    genBank_t *g = (genBank_t *) malloc(sizeof (genBank_t));
    g->locusName = NULL;
    g->cds_number = 0;
    g->cds = NULL;
    return g;
}

void freeGenBank(genBank_t *g) {
    int i;
    if (g != NULL) {
        if (g->locusName) free(g->locusName);
        if (g->cds != NULL) {
            for (i = 0; i < g->cds_number; i++) {
                freeCDS(g->cds + i);
            }
            free(g->cds);
        }
        free(g);
    }
}

genBank_t *initGenBank(int gi, char *glocusName, int taxId, char *proteinId, char *locusName, int pFrom, int pTo, char *cog, char *protClust) {
    genBank_t *g = (genBank_t *) allocGenBank();
    g->gi = gi;
    g->locusName = strdup(glocusName);
    g->taxId = taxId;
    g->cds_number = 1;
    g->cds = initCDS(proteinId, locusName, pFrom, pTo, cog, protClust);
    return g;
}

void printCDS(cds_t * c) {
    int i;
    if (c != NULL) {
        printf("\tProteinId: %s\n", c->proteinId);
        printf("\tLocus: %s\n", c->locusName);
        printf("\tHits: %d\n", c->hits);
        if (c->cog != NULL) {
            for (i = 0; i < c->cog_number; i++) {
                printf("\t\tCOG: %s\n", c->cog[i]);
            }
        }
        if (c->protClust != NULL) {
            for (i = 0; i < c->protClust_number; i++) {
                printf("\t\tPRO: %s\n", c->protClust[i]);
            }
        }
        printf("\t\tLocation: %d -- %d\n", c->pFrom, c->pTo);
        fflush(stdout);
    }
}

void printGenBank(genBank_t *g) {
    int i;
    if (g != NULL) {
        printf("Gi: %d\n", g->gi);
        fflush(NULL);
        printf("Locus: %s\n", g->locusName);
        fflush(NULL);
        printf("TaxId: %d\n", g->taxId);
        fflush(NULL);
        if (g->cds != NULL) {
            if (g->cds_number > 1) {
                qsort(g->cds, g->cds_number, sizeof (cds_t), cmpCDSSort);
            }
            for (i = 0; i < g->cds_number; i++) {
                printCDS(&(g->cds[i]));
            }
        }
    }
}

void printCDSBinary(FILE *fb, cds_t c) {
    int i, size;
    
    size = strlen(c.proteinId) + 1;
    fwrite(&(size), sizeof (int), 1, fb);
    fwrite(c.proteinId, sizeof (char), size, fb);
    size = strlen(c.locusName) + 1;
    fwrite(&(size), sizeof (int), 1, fb);
    fwrite(c.locusName, sizeof (char), size, fb);
    fwrite(&(c.cog_number), sizeof (int), 1, fb);
    fwrite(&(c.protClust_number), sizeof (int), 1, fb);
    fwrite(&(c.pFrom), sizeof (int), 1, fb);
    fwrite(&(c.pTo), sizeof (int), 1, fb);
    fwrite(&(c.hits), sizeof (int), 1, fb);
    for (i = 0; i < c.cog_number; i++) {
        size = strlen(c.cog[i]) + 1;
        fwrite(&(size), sizeof (int), 1, fb);
        fwrite(c.cog[i], sizeof (char), size, fb);
    }
    for (i = 0; i < c.protClust_number; i++) {
        size = strlen(c.protClust[i]) + 1;
        fwrite(&(size), sizeof (int), 1, fb);
        fwrite(c.protClust[i], sizeof (char), size, fb);
    }
}

cds_t *readCDSBinary(FILE *fb, int number) {
    int i, j, size;

    cds_t *c = malloc(sizeof (cds_t) * number);

    for (i = 0; i < number; i++) {
        fread(&size, sizeof (int), 1, fb);
        c[i].proteinId = malloc(sizeof (char) * size);
        fread(c[i].proteinId, sizeof (char), size, fb);
        fread(&size, sizeof (int), 1, fb);
        c[i].locusName = malloc(sizeof (char) * size);
        fread(c[i].locusName, sizeof (char), size, fb);
        fread(&(c[i].cog_number), sizeof (int), 1, fb);
        fread(&(c[i].protClust_number), sizeof (int), 1, fb);
        fread(&(c[i].pFrom), sizeof (int), 1, fb);
        fread(&(c[i].pTo), sizeof (int), 1, fb);
        fread(&(c[i].hits), sizeof (int), 1, fb);
        if (c[i].cog_number == 0) {
            c[i].cog = NULL;
        } else {
            c[i].cog = malloc(sizeof (char **) * c[i].cog_number);
            for (j = 0; j < c[i].cog_number; j++) {
                fread(&size, sizeof (int), 1, fb);
                c[i].cog[j] = malloc(sizeof (char) * size);
                fread(c[i].cog[j], sizeof (char), size, fb);
            }
        }
        if (c[i].protClust_number == 0) {
            c[i].protClust = NULL;
        } else {
            c[i].protClust = malloc(sizeof (char **) * c[i].protClust_number);
            for (j = 0; j < c[i].protClust_number; j++) {
                fread(&size, sizeof (int), 1, fb);
                c[i].protClust[j] = malloc(sizeof (char) * size);
                fread(c[i].protClust[j], sizeof (char), size, fb);
            }
        }
    }
    return c;
}

void printGenBankBinary(FILE *fb, genBank_t *g) {
    int i, size;
    if (g != NULL) {
        fwrite(&(g->gi), sizeof (int), 1, fb);
        size = strlen(g->locusName) + 1;
        fwrite(&(size), sizeof (int), 1, fb);
        fwrite(g->locusName, sizeof (char), size, fb);
        fwrite(&(g->taxId), sizeof (int), 1, fb);
        fwrite(&(g->cds_number), sizeof (int), 1, fb);
        if (g->cds != NULL) {
            if (g->cds_number > 1) {
                qsort(g->cds, g->cds_number, sizeof (cds_t), cmpCDSSort);
            }
            for (i = 0; i < g->cds_number; i++) {
                printCDSBinary(fb, g->cds[i]);
            }
        }
    }
}

genBank_t *readGenBankBinary(FILE *fb) {
    int size;

    genBank_t *g = allocGenBank();
    fread(&(g->gi), sizeof (int), 1, fb);
    fread(&size, sizeof (int), 1, fb);
    g->locusName = malloc(sizeof (char) * size);
    fread(g->locusName, sizeof (char), size, fb);
    fread(&(g->taxId), sizeof (int), 1, fb);
    fread(&(g->cds_number), sizeof (int), 1, fb);
    g->cds = readCDSBinary(fb, g->cds_number);
    return g;
}

genBank_t *readGenBankBinaryOffSet(FILE *fb, off_t offset) {
    fseeko(fb, offset, SEEK_SET);
    return readGenBankBinary(fb);
}

int assignGenes(node **root, genBank_t **gs, int *count, FILE *fb, off_t offset, int pFrom, int pTo) {
    int i, j, k, *out, *index;
    record *rec;
    genBank_t *g = readGenBankBinaryOffSet(fb, offset);

    if (pFrom > pTo) {
        i = pFrom;
        pFrom = pTo;
        pTo = i;
    }

    j = 0;
    if (g != NULL) {
        cds_t *c = g->cds;
        for (i = 0; i < g->cds_number && c[i].pFrom <= pTo; i++) {
            if ((c[i].pFrom <= pFrom && c[i].pTo >= pFrom) ||
                    (c[i].pFrom <= pTo && c[i].pTo >= pTo) ||
                    (c[i].pFrom >= pFrom && c[i].pTo <= pTo) ||
                    (c[i].pFrom <= pFrom && c[i].pTo >= pTo)) {

                if ((rec = find(*root, g->gi, false)) != NULL) {
                    index = (int *) rec->value;
                    for (k = 0; k < (*gs)[*index].cds_number; k++) {
                        if (cmpCDS(&((*gs)[*index].cds[k]), &(c[i])) == 0) {
                            break;
                        }
                    }
                    if (k >= (*gs)[*index].cds_number) {
                        cdsdup(&((*gs)[*index].cds), &((*gs)[*index].cds_number), c[i]);
                        (*gs)[(*index)].cds[(*gs)[*index].cds_number - 1].hits = 1;
                        j++;
                    } else {
                        (*gs)[*index].cds[k].hits++;
                    }
                } else {
                    out = (int *) malloc(sizeof (int));
                    *out = *count;
                    *root = insert(*root, g->gi, out);

                    (*gs) = (genBank_t *) realloc((*gs), sizeof (genBank_t) * (*out + 1));
                    (*gs)[*out].gi = g->gi;
                    (*gs)[*out].taxId = g->taxId;
                    (*gs)[*out].locusName = strdup(g->locusName);
                    (*gs)[*out].cds_number = 0;
                    (*gs)[*out].cds = NULL;
                    cdsdup(&((*gs)[*out].cds), &((*gs)[*out].cds_number), c[i]);
                    (*gs)[*out].cds[0].hits = 1;
                    (*count)++;
                    j++;
                }
            }
        }
        freeGenBank(g);
    }
    
    return j;
}

