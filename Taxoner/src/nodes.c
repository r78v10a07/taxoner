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
 *  Created new nodes for the taxon IDs stored in TaxonNodes * nodes
 */
void AllocateNodes(void) {
    int i = 0;
    nodes = (struct TaxonNodes *)calloc(db->maxNode + 1, sizeof(struct TaxonNodes));

    if(parms->verbose == 1) printf("Allocating %d nodes into memory\n", db->maxNode);

    for(i = 0; i < db->maxNode; i++) {
        nodes[i].curr = 0;
        nodes[i].next = 0;
        nodes[i].ok = 0;
    }
}

/*
 *  Looks for largest node
 */
void GetMaxNodes(char * infile)
{
    FILE * handle = fopen(infile, "r+");
    char line[BUFSIZ];
    int tempnode = 0;

    if(parms->verbose == 1) printf("Searching for highes node in %s\n", parms->names);

    while(fgets(line, sizeof(line), handle))
    {
        tempnode = atoi(line);

        if(tempnode > db->maxNode)
            db->maxNode = tempnode + 1;
    }

    if(parms->verbose == 1)
        printf("Highest taxon (node): %d\n", db->maxNode - 1);

    if(db->maxNode <= 0) {
        printf("ERROR: No taxon IDs found in %s\n", parms->names);
        FreeInputs();
        exit(EXIT_SUCCESS);
    }

    fclose(handle);
}

/*
 *  Fills in taxon node data
 */
void FillNodes(char * infile)
{
    FILE * handle = fopen(infile, "r+");
    char line[BUFSIZ];
    int tempnode;
    int token = 0;
    char * ptr;

    if(parms->verbose == 1) printf("Filling nodes into memory\n");

    while(fgets(line, sizeof(line), handle))
    {
        ptr = strtok(line, "|");
        token = 0;

        while(ptr != NULL)
        {
            if(token == 0)
                tempnode = atoi(ptr);

            if(token == 1)
            {
                nodes[tempnode].curr = tempnode;
                nodes[tempnode].next = atoi(ptr);
                nodes[tempnode].ok = 1;
            }

            ptr = strtok(NULL, "|");
            token++;
        }
    }

    fclose(handle);
}

/*
 *  Removes unwanted taxon IDs
 */
void RemoveUnwantedTaxons(void) {
    FILE * handle = fopen(parms->removeTaxons, "r+");
    char line[BUFSIZ];

    if(parms->verbose == 1) printf("Importing taxons that should be skipped from %s\n", parms->removeTaxons);

    while(fgets(line, sizeof(line), handle)) {
        if(atoi(line) > 0 && atoi(line) < db->maxNode)
            nodes[atoi(line)].ok = 0;
    }

    fclose(handle);
}

/*
 *  Nodes import
 */
void ImportNodes(void) {
    db->maxNode = 0;
    if(parms->verbose == 1) printf("\n***************************************\n");
    if(parms->verbose == 1) printf("\tImporting nodes");
    if(parms->verbose == 1) printf("\n***************************************\n");
    GetMaxNodes(parms->names);
    AllocateNodes();
    FillNodes(parms->names);
    if(parms->filterVirus == 1){
        nodes[10239].ok = 0;
        nodes[28384].ok = 0;
    }

    if(parms->removeTaxons)RemoveUnwantedTaxons();
    if(parms->verbose == 1) printf("\n***************************************\n");
}
