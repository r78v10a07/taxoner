#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <getopt.h>
#include <time.h>
#include <stdbool.h>
#include <pthread.h>
#include <unistd.h>
#include "main.h"

extern struct params * parms;
extern struct database * db;
extern struct TaxonNodes * nodes;

struct ReadChunk {
    char * data; //SAM line
    int datalength; //SAM length
    int accession; //accession id fom ncbi
    int seqstart; //starting position of read
    int seqlen; //length of read
    float as; //alignment score
    int taxon; //taxon id (from ncbi)
    int iter; //iterator at multithread option
    int element;
    int pos;
    char * chel; //pointer position in *data
    int tempPos;
    int nameElement;
    struct ReadChunk * node;
    struct ReadChunk * next;
    struct ReadChunk * prev;
    struct ReadChunk * prevNode;
};

struct ReadChunk * RMaster = NULL;
struct ReadChunk * Rhead = NULL;
struct ReadChunk * Rcurr = NULL;

int chunkSize = 40000; //number of reads for each thread
int AllTax = 100; //max taxons while multithread reading

int CheckIfAligned(char * input) {
    if (strstr(input, "AS:i") != NULL)
        return 0;

    return 1;
}

int CompareStrings(char * s1, char * s2){
    if(strlen(s1) == strlen(s2)) {
        if(strncmp(s1, s2, strlen(s1)) == 0)
            return 1;
    }

    return 0;
}

char * ReturnName(char * input, char * dest) {
    int j = 0;
    int len = strlen(input) - 1;

    if (dest == NULL)
        free(dest);

    dest = (char *) calloc(len + 1, sizeof (char));

    while (input[j] != '\t' && input[j] != ' ' && input[j] != '\n' && j < len) {
        dest[j] = input[j];
        j++;
    }

    return dest;
}

struct ReadChunk * UninitializeStruct(struct ReadChunk * psr) {
    psr->accession = -1;
    psr->as = -1;
    psr->chel = NULL;
    psr->datalength = -1;
    psr->element = -1;
    psr->nameElement = -1;
    psr->seqlen = -1;
    psr->seqstart = -1;
    psr->taxon = -1;

    return psr;
}
//FreeReadStruct(struct ReadChunk * abc, FILE * output)
//Reads ReadStruct linked list, and prints results to file

void FreeReadStruct(struct ReadChunk * abc, FILE * output) {
    struct ReadChunk * temp = NULL;
    struct ReadChunk * Direct = NULL;
    char * tempOut = NULL;

    while (abc != NULL) {
        Direct = abc;
        abc = abc->node;

        while (Direct != NULL) {
            temp = Direct;
            Direct = Direct->next;
            tempOut = ReturnName(temp->data, tempOut);

            if (temp->as > 0) //check if alignment score is ok
                fprintf(output, "%s\t%d\t%d\t%.3f\t%d\t%d\n", tempOut, temp->taxon, temp->accession, temp->as, temp->seqstart, (temp->seqlen + temp->seqstart));

            if (temp->data != NULL)
                free(temp->data);

            if (temp != NULL)
                free(temp);
        }
    }
    if (tempOut) {
        free(tempOut);
    }
}

void NewReadStruct(void) {
    RMaster = (struct ReadChunk *) calloc(1, sizeof (struct ReadChunk));
    Rhead = RMaster;
    Rcurr = RMaster;
    Rhead->prevNode = NULL;
    return;
}

void AddReadStructNode(char * input, int nnum) {
    int slen = strlen(input);
    if (RMaster == NULL)
        NewReadStruct();

    else {
        Rhead->node = (struct ReadChunk *) calloc(1, sizeof (struct ReadChunk));
        Rhead->node->prevNode = Rhead;
        Rhead = Rhead->node;
        Rcurr = Rhead;
    }

    Rhead = UninitializeStruct(Rhead);

    Rhead->next = NULL;
    Rhead->prev = NULL;
    Rhead->data = NULL;
    Rhead->node = NULL;

    Rcurr->data = (char *) calloc(slen, sizeof (char));
    strncpy(Rcurr->data, input, slen - 1);

    Rcurr->nameElement = nnum;
}

void AddReadElement(char * input, int nnum) {
    int slen = strlen(input);
    Rcurr->next = (struct ReadChunk *) calloc(1, sizeof (struct ReadChunk));
    Rcurr = Rcurr->next;

    Rcurr->data = (char *) calloc(slen, sizeof (char));
    strncpy(Rcurr->data, input, slen - 1);

    Rcurr = UninitializeStruct(Rcurr);
    Rcurr->nameElement = nnum;
}

//GetToken(char * source, char delim, int length, int prev)
//looks for token specified by 'delim' in 'source'. 'Prev' tells function to continue
//from specified position in string.
//Returns position where 'delim' can be found in string

int GetToken(char * source, char delim, int length, int prev) {
    int i;

    if (prev >= length)
        return -1;

    for (i = prev; i < length; i++) {
        if (source[i] == delim)
            return i + 1;
    }

    return -1;
}

struct ReadChunk * GoThroughData(struct ReadChunk * psr) {
    struct ReadChunk * tempStruct = psr;
    int deleteElement = 0;

    if (tempStruct->taxon == -1 || tempStruct->taxon >= db->maxNode)
        deleteElement = 1;

    else {
        if (nodes[tempStruct->taxon].next == 0)
            deleteElement = 1;
    }

    if (deleteElement == 1) {
        psr->nameElement = -1;
    }


    return psr;
}


void * MultiParser(void * voidA) {
    struct ReadChunk * psr = (struct ReadChunk *) voidA;

    while (psr != NULL) {
        psr->pos = 0;
        psr->iter = 0;
        psr->datalength = strlen(psr->data);
        psr->element = 0;

        while (psr->pos > -1) {
            if (psr->iter == 2) {
                psr->data += psr->pos;
                psr->accession = atoi(psr->data);
                psr->data -= psr->pos;

                psr->pos = GetToken(psr->data, ';', psr->datalength, psr->pos);

                if (psr->pos > 0) {
                    psr->data += psr->pos;
                    psr->taxon = atoi(psr->data);
                    psr->data -= psr->pos;
                } else
                    psr->taxon = -1;
            } else if (psr->iter == 3) {
                if (psr->pos > 0) {
                    psr->data += psr->pos;
                    psr->seqstart = atoi(psr->data);
                    psr->data -= psr->pos;
                } else
                    psr->seqstart = -1;
            } else if (psr->iter == 9) {
                if (psr->pos > 0) {
                    psr->data += psr->pos;
                    psr->seqlen = GetToken(psr->data, '\t', psr->datalength - psr->pos, 0) - 1;
                    psr->data -= psr->pos;
                }
            } else if (psr->iter > 10) {
                psr->data += psr->pos;
                psr->chel = strstr(psr->data, "AS:i:");

                if (psr->chel != NULL) {
                    psr->chel += 5;
                    psr->as = atof(psr->chel);
                    psr->as = ((float) psr->seqlen + psr->as) / ((float) psr->seqlen);
                    psr->chel -= 5;
                }

                psr->data -= psr->pos;
                psr->pos = -1;
            }

            if (psr->pos != -1)
                psr->pos = GetToken(psr->data, '\t', psr->datalength, psr->pos);

            psr->iter++;
        }

        psr = GoThroughData(psr);
        psr = psr->next;
    }

    return NULL;
}


void ParseDataMulti(void) {
    struct ReadChunk * tempStruct = RMaster;
    pthread_t thread_id[parms->threads];
    int i = 0;

    while (tempStruct != NULL) { //calls parser for each thread
        pthread_create(&thread_id[i], NULL, &MultiParser, tempStruct);
        tempStruct = tempStruct->node;
        i++;
    }

    i = 0;
    tempStruct = RMaster;

    while (tempStruct != NULL) { //joins threads
        pthread_join(thread_id[i], NULL);
        tempStruct = tempStruct->node;
        i++;
    }
}

//ReturnNeighbor(int a, int b)
//looks for nearest neighbor for taxons 'a' and 'b', and returns
//nearest neighbor id (integer)

int ReturnNeighbor(int a, int b) {
    int * s1, * s2;
    int l1, l2, i, j;

    int temptax = a;

    if (a > db->maxNode || b > db->maxNode)
        return -1;

    if (b == 0 || nodes[b].next == 0)
        return a;

    s1 = (int *) calloc(AllTax * 2, sizeof (int));
    s2 = (int *) calloc(AllTax * 2, sizeof (int));

    l1 = 0;
    l2 = 0;

    while (nodes[temptax].curr != nodes[temptax].next && l1 < AllTax) { //store taxID lineage in s1
        s1[l1] = temptax;
        temptax = nodes[temptax].next;
        l1++;
    }
    temptax = b;
    while (nodes[temptax].curr != nodes[temptax].curr && l2 < AllTax) { //store taxID lineage in s2
        s2[l2] = temptax;
        temptax = nodes[temptax].next;
        l2++;
    }

    if (l1 > l2) {
        j = l1 - 1;

        for (i = l2 - 1; i > -1; i--) {
            if (s2[i] == s1[j])
                temptax = s2[i];

            else
                i = -1;

            j--;
        }
    } else {
        j = l2 - 1;

        for (i = l1 - 1; i > -1; i--) {
            if (s2[j] == s1[i])
                temptax = s1[i];

            else
                i = -1;

            j--;
        }
    }

    if (s1) free(s1);
    if (s2) free(s2);

    return temptax;
}

/*isVirus: Check if taxon id is a virus or other sequences, and exclude from results*/
int isVirus(int a) {
    int temptax = a;
    int l1 = 0;

    if (a == -1)
        return 1;

    if (parms->filterVirus == 0)
        return 0;

    while (nodes[temptax].curr != nodes[temptax].next && l1 < AllTax) {
        if (nodes[temptax].ok == 0)
            return 1;
        temptax = nodes[temptax].next;
        l1++;
    }

    return 0;
}

//ReadChunk * RemoveVirus(struct ReadChunk * psr)
//Remove virus reads form ReadChunk linked list

struct ReadChunk * RemoveVirus(struct ReadChunk * psr) {
    struct ReadChunk * ptr = psr;
    if(parms->removeTaxons == NULL && parms->filterVirus == 0) return psr;

    while (ptr != NULL) {
        if (isVirus(ptr->taxon) == 1) {
            ptr->as = -1;
        }

        ptr = ptr->next;
    }

    return psr;
}

//ReadChunk * FreeNextBranch(struct ReadChunk * psr)
//Frees found branch

struct ReadChunk * FreeNextBranch(struct ReadChunk * psr) {
    struct ReadChunk * tempStruct = psr->next;
    psr->next = NULL;
    psr->next = tempStruct->next;

    if (tempStruct->data) free(tempStruct->data);
    if (tempStruct) free(tempStruct);

    return psr;
}

//ParallelNeighbor(void * voidA)
//looks for nearest neighbors in reads stored in '*voidA'

void * ParallelNeighbor(void * voidA) {
    struct ReadChunk * psr = (struct ReadChunk *) voidA;

    psr = RemoveVirus(psr);
    while (psr != NULL) {
        if (psr->next != NULL && psr->iter > 3) {
            if (psr->next->nameElement == psr->nameElement && psr->as != -1) {
                if (psr->next->next != NULL) {
                    if (psr->next->as > parms->neighborScore * psr->as)
                        psr->taxon = ReturnNeighbor(psr->taxon, psr->next->taxon);

                    psr = FreeNextBranch(psr);

                } else {
                    if (psr->next->as > parms->neighborScore * psr->as)
                        psr->taxon = ReturnNeighbor(psr->taxon, psr->next->taxon);

                    psr = FreeNextBranch(psr);
                }

            } else
                psr = psr->next;
        } else
            psr = psr->next;
    }

    return NULL;
}

//InitializeNearestMulti(struct ReadChunk * psr)
//calls parallel nearest neighbor identification

void InitializeNearestMulti(struct ReadChunk * psr) {
    struct ReadChunk * tempStruct = psr;
    pthread_t thread_id[parms->threads];
    int i = 0;

    while (tempStruct != NULL) { //calls threads
        pthread_create(&thread_id[i], NULL, &ParallelNeighbor, tempStruct);
        tempStruct = tempStruct->node;
        i++;
    }

    i = 0;
    tempStruct = psr;

    while (tempStruct != NULL) { //joins threads
        pthread_join(thread_id[i], NULL);
        tempStruct = tempStruct->node;
        i++;
    }
}

void SingleThread(void) {
    struct ReadChunk * psr = RMaster;

    psr = RemoveVirus(psr);
    while (psr != NULL) {
        if (psr->next != NULL && psr->iter > 3) {
            if (psr->next->nameElement == psr->nameElement && psr->as != -1) {
                if (psr->next->next != NULL) {
                    if (psr->next->as > parms->neighborScore * psr->as)
                        psr->taxon = ReturnNeighbor(psr->taxon, psr->next->taxon);

                    psr = FreeNextBranch(psr);

                } else {
                    if (psr->next->as > parms->neighborScore * psr->as)
                        psr->taxon = ReturnNeighbor(psr->taxon, psr->next->taxon);

                    psr = FreeNextBranch(psr);
                }

            } else
                psr = psr->next;
        } else
            psr = psr->next;
    }
}



void ReadSamFile(char * infile, char * outfile) {
    FILE * handler = fopen(infile, "r");
    FILE * out = fopen(outfile, "w");
    char line[BUFSIZ];
    int EqName = 0;
    char * currName = NULL;
    char * prevName = NULL;
    int i;
    int num = 0;

    while (fgets(line, sizeof (line), handler)) {
        if (line[0] != '@') {
            if (CheckIfAligned(line) == 0) {
                currName = ReturnName(line, currName); //get name of read

                if (RMaster == NULL) { //crate struct to store data
                    AddReadStructNode(line, EqName);
                } else if (CompareStrings(currName, prevName) == 1) {
                    AddReadElement(line, EqName); //add elemnent to struct
                } else {
                    EqName++;
                    if (i >= chunkSize) { //if read count reaches chunckSize
                        i = 0;
                        num++; //add thread

                        if (num == parms->threads) { //if thread number equals to number of threads
                            i = 0;
                            num = 0;
                            EqName = 0;

                            ParseDataMulti(); //parse data
                            InitializeNearestMulti(RMaster); //get neighbor
                            SingleThread();
                            FreeReadStruct(RMaster, out); //print results and free data
                            RMaster = NULL;
                        }

                        AddReadStructNode(line, EqName); //create new struct
                    } else if (RMaster == NULL)
                        AddReadStructNode(line, EqName); //add new therad struct

                    else
                        AddReadElement(line, EqName);
                }

                prevName = ReturnName(line, prevName);
                i++;
            }
        }
    }

    InitializeNearestMulti(RMaster); //get neighbor
    fclose(handler);
    FreeReadStruct(RMaster, out);
    fclose(out);

    RMaster = NULL;
    Rcurr = NULL;
    Rhead = NULL;

    if (currName) {
        free(currName);
    }
    if (prevName) {
        free(prevName);
    }
}

void StartNeighbor(void) {
    int i = 0;
    char * inFile = NULL;
    char * outFile = NULL;
    for(i = 0; i < db->indexes; i++) {
        if(inFile) free(inFile);

        inFile = (char *)calloc(strlen(db->names[i]) + strlen(parms->alignOut) + 6, sizeof(char));
        strcpy(inFile, parms->alignOut);
        strcat(inFile, "/");
        strcat(inFile, db->names[i]);
        strcat(inFile, ".sam");

        printf("%s\n", inFile);
        outFile = (char *)calloc(strlen(inFile) + 5, sizeof(char));
        strcpy(outFile, inFile);
        strcat(outFile, ".aln");
        ReadSamFile(inFile, outFile);
    }

    if(inFile) free(inFile);
    if(outFile) free(outFile);
}

void GetNeighbors(void) {
    if(parms->verbose == 1) printf("\n***************************************\n");
    if(parms->verbose == 1) printf("\tNearest neighbor");
    if(parms->verbose == 1) printf("\n***************************************\n");
    StartNeighbor();
    if(parms->verbose == 1) printf("\n***************************************\n");
}
