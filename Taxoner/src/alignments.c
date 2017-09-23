#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <getopt.h>
#include <time.h>
#include <stdbool.h>
#include <pthread.h>
#include <unistd.h>
#include "main.h"
#include <ctype.h>

extern struct params * parms;
extern struct database * db;
extern struct TaxonNodes * nodes;

/*
 *  Created bowtie2 command, and returns pointer
 */
char * ReturnBowtieCommand(char * InputIndex, char * indexName) {
    char * command = (char *)calloc(10000, sizeof(char));
    char temp[1000];
    char hits[20];
    if(parms->bowtieExec) strcpy(temp, parms->bowtieExec);
    else strcpy(temp, "bowtie2");

    if(parms->fastaInput == 1) strcat(temp, " -f");

    if(parms->bowtieAllHits == 1) strcpy(hits, "-a");
    else sprintf(hits, "-k %d", parms->bowtieMaxHits);

    if(parms->bowtie2_params) {
        strcat(temp, " ");
        strcat(temp, parms->bowtie2_params);
    }

    if(parms->pairedReads)
        sprintf(command, "%s -p %d %s --no-mixed --no-discordant -I %d -X %d -x %s -1 %s -2 %s > %s/%s.sam", temp, parms->threads, hits, parms->minInsert, parms->maxInsert, InputIndex, parms->input, parms->pairedReads, parms->alignOut, indexName);

    else
        sprintf(command, "%s -p %d %s -x %s %s > %s/%s.sam", temp, parms->threads, hits,InputIndex, parms->input, parms->alignOut, indexName);

    return command;
}

char * ReturnBowtieCommandHost(char * InputIndex, char * indexName) {
    char * command = (char *)calloc(10000, sizeof(char));
    char temp[1000];
    char hits[20];
    if(parms->bowtieExec) strcpy(temp, parms->bowtieExec);
    else strcpy(temp, "bowtie2");

    if(parms->fastaInput == 1) strcat(temp, " -f");

    if(parms->bowtieAllHits == 1) strcpy(hits, "-a");
    else sprintf(hits, "-k %d", parms->bowtieMaxHits);

    if(parms->bowtie2_params) {
        strcat(temp, " ");
        strcat(temp, parms->bowtie2_params);
    }

    if(parms->pairedReads) {
        sprintf(command, "%s -p %d --un-conc %s/unaligned --no-mixed --no-discordant -I %d -X %d -x %s -1 %s -2 %s > %s/%s.sam", temp, parms->threads, parms->alignOut, parms->minInsert, parms->maxInsert, InputIndex, parms->input, parms->pairedReads, parms->alignOut, indexName);
    }

    else {
        sprintf(command, "%s -p %d --un %s/unaligned -x %s %s > %s/%s.sam", temp, parms->threads, parms->alignOut, InputIndex, parms->input, parms->alignOut, indexName);
    }

    return command;
}


/*
 *  Aligns reads
 */
void StartAlignment(void) {
    int i = 0;
    char * bowtieCommand = NULL;

    if(parms->hostFilter == 0 && parms->hostindex != NULL) {
        if(bowtieCommand) free(bowtieCommand);
        bowtieCommand = ReturnBowtieCommandHost(parms->hostindex, "host");
        if(parms->verbose == 1) printf("Aligning reads to host %s\n", parms->hostindex);
        if(parms->verbose == 1) printf("Bowtie2 command: %s\n", bowtieCommand);
        if(parms->verbose == 1) printf("Alignment may take a long time.\n");
        if(parms->verbose == 1) fflush(NULL);
        printf("Bowtie2 output:\n");
        system(bowtieCommand);
        printf("\n");
        exit(EXIT_FAILURE);
    }

    for(i = 0; i < db->indexes; i++) {
        if(bowtieCommand) free(bowtieCommand);
        bowtieCommand = ReturnBowtieCommand(db->paths[i], db->names[i]);
        if(parms->verbose == 1) printf("\nIndex [ %d / %d ]: %s\n", i + 1 , db->indexes, db->names[i]);
        if(parms->verbose == 1) printf("Bowtie2 command: %s\n", bowtieCommand);
        if(parms->verbose == 1) printf("Alignment may take a long time.\n");
        if(parms->verbose == 1) fflush(NULL);
        printf("Bowtie2 output:\n");
        system(bowtieCommand);
        printf("\n");
    }

    if(bowtieCommand) free(bowtieCommand);
}

void AlignReads(void) {
    if(parms->verbose == 1) printf("\n***************************************\n");
    if(parms->verbose == 1) printf("\tAligning reads");
    if(parms->verbose == 1) printf("\n***************************************\n");
    StartAlignment();
    if(parms->verbose == 1) printf("\n***************************************\n");
}
