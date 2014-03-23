#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <dirent.h>
#include "inputs.c"
#include "SAMParse.c"
#include "nodes.c"
#include <time.h>
#include <pthread.h>
#include <unistd.h>

/*
Authors: Lorinc Pongor, Roberto Vera, Balazs Ligeti
email: pongorlorinc@gmail.com
date: march 23, 2014

Purpose: Fast microbila classification of NGS reads
Language: C
*/

//another comment
//struct taxons: Store information about taxonomy
struct Taxons {
    int id;
    char * path;
    struct Taxons * next;
};

//struct REadChunk: store SAM information in memory
struct ReadChunk {
    char * data; //SAM line
    int datalength; //SAM length
    int accession; //accession id fom ncbi
    int seqstart;  //starting position of read
    int seqlen; //length of read
    float as;   //alignment score
    int taxon;  //taxon id (from ncbi)
    int iter;   //iterator at multithread option
    int element;
    int pos;
    char * chel;    //pointer position in *data
    int tempPos;
    int nameElement;
    struct ReadChunk * node;
    struct ReadChunk * next;
    struct ReadChunk * prev;
    struct ReadChunk * prevNode;
};

typedef struct {
    char * path;
} rint;

rint * tax = NULL;

/*multithreaded option*/
int TotalSams = 0;
int chunkSize = 40000; //number of reads for each thread
int AllTax = 100; //max taxons while multithread reading

int tempnode; //temporary storing of taxonomy

char * NeighborName = NULL;

struct Taxons * thead = NULL;
struct Taxons * tcurr = NULL;

struct ReadChunk * RMaster = NULL;
struct ReadChunk * Rhead = NULL;
struct ReadChunk * Rcurr = NULL;


/*retrieveTime: Get current time, and store in times variable*/
void retrieveTime(void) {
    time(&timer);
    tm_info = localtime(&timer);
    strftime(times, 30, "%Y-%m-%d %H:%M:%S", tm_info);
}

/*void taxonomyFree(void)
{
    int i;

    if(tax == NULL)
        return;

    for(i = 0; i < maxTax; i++)
    {
        if(tax[i].path != NULL)
            free(tax[i].path);
    }

    free(tax);
}*/


/*isVirus: Check if taxon id is a virus or other sequences, and exclude from results*/
int isVirus(int a) {
    int temptax = a;
    int l1 = 0;
    int l2 = 0;

    if (a == -1)
        return 1;

    if(FilterVirus == 0)
        return 0;

    while (n[temptax].id != n[temptax].nextId && l1 < AllTax) {
        if (temptax == 10239 || temptax == 28384)
            return 1;
        temptax = n[temptax].nextId;
        l1++;
    }

    return 0;
}

//ReadChunk * RemoveVirus(struct ReadChunk * psr)
//Remove virus reads form ReadChunk linked list
struct ReadChunk * RemoveVirus(struct ReadChunk * psr) {
    int counts = 0;
    struct ReadChunk * tempStruct = psr->next;
    struct ReadChunk * ptr = psr;

    while (ptr != NULL) {
        //printf("%s\n", psr->data);
        if (isVirus(ptr->taxon) == 1) {
            ptr->as = -1;
        }

        ptr = ptr->next;
    }

    return psr;
}

//CheckFile(char * input)
//Checks is input is an existing file
//returns 1 if exists, 0 if not
int CheckFile(char * input)
{
    if( access( input, F_OK ) != -1 ) {
        return 1;
    } else {
        return 0;
    }
}

//createTaxons(int num, char * p)
//Created Taxons struct
void createTaxons(int num, char * p) {
    struct Taxons *psr = (struct Taxons *) calloc(1, sizeof (struct Taxons));

    if (psr == NULL) {
        printf("\nNode creation failed\n");
        exit(EXIT_FAILURE);
    }

    thead = tcurr = psr;

    thead->id = num;
    thead->path = CopyString(p, strlen(p));
}

//AddTaxons(int num, char * p)
//If taxons list exists, adds new element to linked list
void AddTaxons(int num, char * p) {
    struct Taxons *psr = (struct Taxons *) calloc(1, sizeof (struct Taxons));

    if (psr == NULL) {
        printf("\nNode creation failed\n");
        exit(EXIT_FAILURE);
    }

    psr->next = NULL;
    tcurr->next = psr;

    psr->id = num;
    psr->path = CopyString(p, strlen(p));
    tcurr = psr;
}

//PrintTaxons(void)
//prints data stored in Taxons linked list
void PrintTaxons(void) {
    while (thead != NULL) {
        tcurr = thead;
        thead = tcurr->next;
        //printf("%d\t%s\n", tcurr->id, tcurr->path);

        if (tcurr->path)free(tcurr->path);
        if (tcurr)free(tcurr);
    }
}

//databaseAlignment(int num, char * infile)
//Aligns reads in 'inFile' to the 'num' element in the database
int databaseAlignment(int num, char * infile) {
    retrieveTime();
    char * partailfulldb;
    char * partialdb;
    char command[BUFSIZ]; //command for bowtie2
    char outf[BUFSIZ]; //temporary
    char Tempcat[BUFSIZ]; //stores string for command

    memset( outf, '\0', sizeof( outf ));
    memset( command, '\0', sizeof( command ));

    partailfulldb = (char *) calloc(strlen(fulldb[num]), sizeof (char));
    strncpy(partailfulldb, fulldb[num], strlen(fulldb[num]) - 6);

    partialdb = (char *) calloc(strlen(dbElements[num]), sizeof (char));
    strncpy(partialdb, dbElements[num], strlen(dbElements[num]) - 6);

    strcat(outf, outaln);
    strcat(outf, "/");
    strcat(outf, partialdb);
    strcat(outf,".sam");

    if(bowtieFolder == NULL)
        strcat(Tempcat, "bowtie2 ");

    else {
        sprintf(Tempcat, "%s", bowtieFolder);
        if(bowtieFolder[strlen(bowtieFolder) - 1] != '/')
            strcat(Tempcat, "/");

        strcat(Tempcat, "./bowtie2 ");
        //printf("Bow: %s\n", Tempcat);
        //exit(EXIT_FAILURE);
    }

    if(inputfasta == 1) //if reads are in fasta format
        strcat(Tempcat, "-f ");

    strcat(Tempcat, "-p"); //number of threads

    if (allHits == 1) {
        if (filterHost == 1 && host != NULL) {
            //printf("1\n");
            if(pairedend == 1) { //if reads are paired-end
                sprintf(command, "%s %d --no-mixed --no-discordant -I %d -X %d -a %s -1 %s -2 %s > %s/%s.sam", Tempcat, threads, minfrag, maxfrag, partailfulldb, infile, pairedFile, outaln, partialdb);
            }

            else {
                sprintf(command, "%s %d -a %s %s/%s > %s/%s.sam", Tempcat, threads, partailfulldb, outaln, infile, outaln, partialdb);
            }

            /*if(inputfasta == 0)
                sprintf(command, "bowtie2 -p %d -a %s %s/%s > %s/%s.sam", threads, partailfulldb, outaln, infile, outaln, partialdb);

            else if(inputfasta == 1)
                sprintf(command, "bowtie2 -f -p %d -a %s %s/%s > %s/%s.sam", threads, partailfulldb, outaln, infile, outaln, partialdb);

            else {
                printf("Could not align reads, unknown format?\n");
                exit(EXIT_FAILURE);
            }*/
        } else {
            //printf("2\n");
            if(pairedend == 1) {
                sprintf(command, "%s %d --no-mixed --no-discordant -I %d -X %d -a %s -1 %s -2 %s > %s/%s.sam", Tempcat, threads, minfrag, maxfrag, partailfulldb, infile, pairedFile, outaln, partialdb);
            }

            else {
                sprintf(command, "%s %d -a %s %s > %s/%s.sam", Tempcat, threads, partailfulldb, outaln, infile, outaln, partialdb);
            }
            /*if(inputfasta == 0)
                sprintf(command, "bowtie2 -p %d -a %s %s > %s/%s.sam", threads, partailfulldb, infile, outaln, partialdb);

            else if(inputfasta == 1)
                sprintf(command, "bowtie2 -f -p %d -a %s %s > %s/%s.sam", threads, partailfulldb, infile, outaln, partialdb);

            else {
                printf("Could not align reads, unknown format?\n");
                exit(EXIT_FAILURE);
            }*/
        }
    } else {
        if (filterHost == 1 && host != NULL) {
            //printf("3\n");

            if(pairedend == 1) {
                sprintf(command, "%s %d --no-mixed --no-discordant -I %d -X %d -k %d %s -1 %s -2 %s > %s/%s.sam", Tempcat, threads, minfrag, maxfrag, bowtieMaxHits, partailfulldb, infile, pairedFile, outaln, partialdb);
            }

            else {
                sprintf(command, "%s %d -k %d %s %s > %s/%s.sam", Tempcat, threads, bowtieMaxHits, partailfulldb, outaln, infile, outaln, partialdb);
            }

            /*if(pairedend == 1) {
                sprintf(command, "%s %d -I %d -X %d -k %d %s -1 %s -2 %s > %s/%s.sam", Tempcat, threads, minfrag, maxfrag, bowtieMaxHits, partailfulldb, infile, pairedFile, outaln, partialdb);
            }

            if(inputfasta == 0)
                sprintf(command, "bowtie2 -p %d -k %d %s %s > %s/%s.sam", threads, bowtieMaxHits, partailfulldb, outaln, infile, outaln, partialdb);

            else if(inputfasta == 1)
                sprintf(command, "bowtie2 -f -p %d -k %d %s %s > %s/%s.sam", threads, bowtieMaxHits, partailfulldb, outaln, infile, outaln, partialdb);

            else {
                printf("Could not align reads, unknown format?\n");
                exit(EXIT_FAILURE);
            }*/
        } else {
            //printf("4\n");
            if(pairedend == 1) {
                sprintf(command, "%s %d --no-mixed --no-discordant -I %d -X %d -k %d %s -1 %s -2 %s > %s/%s.sam", Tempcat, threads, minfrag, maxfrag, bowtieMaxHits, partailfulldb, infile, pairedFile, outaln, partialdb);
            }

            else {
                sprintf(command, "%s %d -k %d %s %s > %s/%s.sam", Tempcat, threads, bowtieMaxHits, partailfulldb, infile, outaln, partialdb);
            }

            //else if(inputfasta == 0)
              //  sprintf(command, "bowtie2 -p %d -k %d %s %s > %s/%s.sam", threads, bowtieMaxHits, partailfulldb, infile, outaln, partialdb);

            //else if(inputfasta == 1)
              //  sprintf(command, "bowtie2 -f -p %d -k %d %s %s > %s/%s.sam", threads, bowtieMaxHits, partailfulldb, infile, outaln, partialdb);

            //else {
              //  printf("Could not align reads, unknown format?\n");
               // exit(EXIT_FAILURE);
            //}
        }
    }

    printf("Command: %s\n", command);
    fflush(NULL);
    printf("Index: %s\nReads: %s\nStart time: %s\n", partialdb, infile, times);
    fflush(NULL);
    printf("Bowtie2 output:\n");
    fflush(NULL);
    system(command); //run command

    if(CheckFile(outf) == 0) {
        printf("Could not find alignment file: %s\nWas the correct format specified for reads? (fstq or fasta?)\n", outf);
        exit(EXIT_FAILURE);
    }

    printf("Getting nearest neighbors from %s database alignment\n", partailfulldb);
    fflush(NULL);
    retrieveTime();
    printf("Finished at: %s\n\n", times);
    fflush(NULL);

    if (partailfulldb) free(partailfulldb);
    if (partialdb) free(partialdb);

    return 0;
}

//hostAlignment(void)
//This function aligns reads to the specified host index
//does not suppport paired-end input reads
int hostAlignment(void) {
    char command[BUFSIZ];
    retrieveTime();
    printf("--------------------------------------\n");
    printf("Host alignment [ %s ]\n", times);
    printf("--------------------------------------\n");
    printf("Index: %s\nReads: %s\n", host, reads);
    fflush(NULL);

    if(inputfasta == 0)
        sprintf(command, "bowtie2 -p %d --un %s/%s %s %s > %s/host.sam", threads, outaln, nohost, host, reads, outaln);

    else
        sprintf(command, "bowtie2 -f -p %d --un %s/%s %s %s > %s/host.sam", threads, outaln, nohost, host, reads, outaln);

    printf("Command: %s\n\n", command);
    printf("Bowtie2 output:\n");
    system(command);
    printf("--------------------------------------\n\n");
    fflush(NULL);
    return 0;
}

//printPreDatabaseInfo(void)
//prints all the index names that are used in analysis
void printPreDatabaseInfo(void) {
    char command[BUFSIZ];
    retrieveTime();
    printf("--------------------------------------\n");
    printf("Db alignment [ %s ]\n", times);
    printf("--------------------------------------\n");
    fflush(NULL);
}

//Alignments(void)
//Calls alignment functions for each bowtie2 index specified
void Alignments(void) {
    int i = 0;

    if (filterHost == 1 && host != NULL) {
        hostAlignment();
        printPreDatabaseInfo();
        for (i = 0; i < dbe; i++) {
            databaseAlignment(i, nohost);
        }

        return; //if alignment is against host, returns
    }


    for (i = 0; i < dbe; i++) {
        databaseAlignment(i, reads);
    }
}

void FreeNames(char ** source) {
    int i = 0;

    for (i = 0; i < dbe; i++) {
        if (source[i]) free(source[i]);
    }

    if (source) free(source);
}

//ReturnName(char * input, char * dest)
//copies read name from SAM line from input to dest
char * ReturnName(char * input, char * dest) {
    int j = 0;
    int len = strlen(input) - 1;

    if (dest)
        free(dest);

    dest = (char *) calloc(len + 1, sizeof (char));

    while (input[j] != '\t' && input[j] != ' ' && input[j] != '\n' && j < len) {
        dest[j] = input[j];
        j++;
    }

    return dest;
}

//FreeReadStruct(struct ReadChunk * abc, FILE * output)
//Reads ReadStruct linked list, and prints results to file
void FreeReadStruct(struct ReadChunk * abc, FILE * output) {
    struct ReadChunk * temp = NULL;
    struct ReadChunk * Direct = NULL;
    int i;
    int printedOut = 0;
    char * tempOut = NULL;

    while (abc != NULL) {
        Direct = abc;
        abc = abc->node;

        while (Direct != NULL) {
            i = 0;
            temp = Direct;
            Direct = Direct->next;
            tempOut = ReturnName(temp->data, tempOut);

            if (temp->as > 0) //check if alignment score is ok
                fprintf(output, "%s\t%d\t%d\t%.3f\t%d\t%d\n", tempOut, temp->taxon, temp->accession, temp->as, temp->seqstart, (temp->seqlen + temp->seqstart));
            //printf( "%s\t%d\t%d\t%.3f\t%d\t%d\n", tempOut, temp->taxon, temp->accession, temp->as, temp->seqstart, temp->seqlen);

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

void FreeReadStruct2(struct ReadChunk * abc) {
    struct ReadChunk * temp = NULL;
    struct ReadChunk * Direct = NULL;
    int i;
    int printedOut = 0;
    char * tempOut = NULL;

    while (abc != NULL) {
        Direct = abc;
        abc = abc->node;

        while (Direct != NULL) {
            i = 0;
            temp = Direct;
            Direct = Direct->next;

            //printf("%s\n", temp->data);
            //if(temp->data)
            //free(temp->data);

            tempOut = ReturnName(temp->data, tempOut);
            //printf("%s\t%d\t%d\t%.3f\t%d\t%d\n", tempOut, temp->taxon, temp->accession, temp->as, temp->seqstart, temp->seqstart + temp->seqlen);

            if (temp->data != NULL)
                free(temp->data);

            if (temp != NULL)
                free(temp);
        }

    }
    if (tempOut) {
        free(tempOut);
    }
    //psr = NULL;
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

//ReturnLength(char * source, char delim, int length, int prev)
//looks for token specified by 'delim' in 'source'. 'Prev' tells function to continue
//from specified position in string.
//Returns position where 'delim' can be found in string
int ReturnLength(char * source, char delim, int length, int prev) {
    int i;
    int retCount = 0;

    for (i = prev; i < length; i++) {
        if (source[i] == delim)
            return retCount;

        retCount++;
    }

    return -1;
}

//ReturnNeighbor(int a, int b)
//looks for nearest neighbor for taxons 'a' and 'b', and returns
//nearest neighbor id (integer)
int ReturnNeighbor(int a, int b) {
    int * s1, * s2;
    int l1, l2, i, j;

    int temptax = a;

    if (a > maxTax || b > maxTax)
        return -1;

    if (b == 0 || n[b].nextId == 0)
        return a;

    s1 = (int *) calloc(AllTax * 2, sizeof (int));
    s2 = (int *) calloc(AllTax * 2, sizeof (int));

    l1 = 0;
    l2 = 0;

    while (n[temptax].id != n[temptax].nextId && l1 < AllTax) { //store taxID lineage in s1
        s1[l1] = temptax;
        temptax = n[temptax].nextId;
        l1++;
    }
    temptax = b;
    while (n[temptax].id != n[temptax].nextId && l2 < AllTax) { //store taxID lineage in s2
        s2[l2] = temptax;
        temptax = n[temptax].nextId;
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

//FreeTaxons(struct ReadChunk * psr) {
//Frees Taxon struct pointed by 'psr'
void FreeTaxons(struct ReadChunk * psr) {
    if (psr->data)
        free(psr->data);

    psr->data = NULL;

    if (psr) free(psr);
}

struct ReadChunk * CheckTaxonId(struct ReadChunk * psr) {
    struct ReadChunk * tempStruct;
    int deleteElement = 0;

    if (psr->taxon == 0 || psr->taxon > maxTax)
        deleteElement = 1;

    else {
        if (n[psr->taxon].nextId == 0)
            deleteElement = 1;
    }

    if (deleteElement == 1) {
        tempStruct = psr;
        psr = psr->next;

        if (tempStruct->data)
            free(tempStruct->data);

        if (tempStruct) free(tempStruct);
    }

    return psr;
}

//CheckIfAligned(char * input)
//checks if SAM line stored in 'input' is aligned
int CheckIfAligned(char * input) {
    if (strstr(input, "AS:i") != NULL)
        return 0;

    return 1;
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
    int asd = 0;
    psr = RemoveVirus(psr);
    while (psr != NULL) {
        if (psr->next != NULL && psr->iter > 3) {
            if (psr->next->nameElement == psr->nameElement && psr->as != -1) {
                if (psr->next->next != NULL) {
                    if (psr->next->as > NSc * psr->as)
                        psr->taxon = ReturnNeighbor(psr->taxon, psr->next->taxon);

                    psr = FreeNextBranch(psr);

                } else {
                    if (psr->next->as > NSc * psr->as)
                        psr->taxon = ReturnNeighbor(psr->taxon, psr->next->taxon);

                    psr = FreeNextBranch(psr);
                }

            } else
                psr = psr->next;
        } else
            psr = psr->next;
    }
}

//ReadChunk * UninitializeStruct(struct ReadChunk * psr)
//Ininitializes data stored in '*psr'
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

struct ReadChunk * GoThroughData(struct ReadChunk * psr) {
    struct ReadChunk * tempStruct = psr;
    int deleteElement = 0;


    if (tempStruct->taxon == -1 || tempStruct->taxon >= maxTax)
        deleteElement = 1;

    else {
        if (n[tempStruct->taxon].nextId == 0)
            deleteElement = 1;
    }


    if (deleteElement == 1) {
        //printf("\n%d (%d)", psr->taxon, maxTax);
        psr->nameElement = -1;
    }


    return psr;
}

//MultiParser(void * voidA)
//Parallel parsing of SAM reads
//ugly, but works
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
}

//InitializeNearestMulti(struct ReadChunk * psr)
//calls parallel nearest neighbor identification
void InitializeNearestMulti(struct ReadChunk * psr) {
    struct ReadChunk * tempStruct = psr;
    pthread_t thread_id[threads];
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

//ParseDataMulti(void)
//Calls multithreaded parsing
void ParseDataMulti(void) {
    struct ReadChunk * tempStruct = RMaster;
    pthread_t thread_id[threads];
    int i = 0;

    while (tempStruct != NULL) { //calls parser for each thread
        pthread_create(&thread_id[i], NULL, &MultiParser, tempStruct);
        tempStruct = tempStruct->node;
        i++;
    }

    i = 0;
    tempStruct = RMaster;

    while (tempStruct != NULL) {    //joins threads
        pthread_join(thread_id[i], NULL);
        tempStruct = tempStruct->node;
        i++;
    }
}

//ReadSam(char * infile, FILE * out)
//Read SAM file specified by 'infile'
void ReadSam(char * infile, FILE * out) {
    FILE * handler = fopen(infile, "r+");
    char line[BUFSIZ];
    int i = 0;
    int j = 0;
    int num = 0;
    int EqName = 0;
    char * currName, * prevName;

    currName = NULL;
    prevName = NULL;

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

                        if (num == threads) { //if thread number equals to number of threads
                            i = 0;
                            num = 0;
                            EqName = 0;

                            ParseDataMulti(); //parse data
                            InitializeNearestMulti(RMaster); //get neighbor
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

    if (currName) {
        free(currName);
    }
    if (prevName) {
        free(prevName);
    }

    ParseDataMulti();
    InitializeNearestMulti(RMaster);
    FreeReadStruct(RMaster, out);
    RMaster = NULL;

    fclose(handler);
}

//-dbPath /home/lorinc/bioinformatica/newtaxon/db/DbCreator/bowtie/ -taxpath /home/lorinc/bioinformatica/newtaxon/db/DbCreator/nodes.dmp -seq wGenome.fastq -o res/ -p 4

void copyFile(char *source, char *dest) {
    int childExitStatus;
    pid_t pid;
    int status;
    if (!source || !dest) {
        /* handle as you wish */
    }

    pid = fork();

    if (pid == 0) { /* child */
        execl("/bin/cp", "/bin/cp", source, dest, (char *) 0);
    } else if (pid < 0) {
        /* error - couldn't start process - you decide how to handle */
    } else {
        /* parent - wait for child - this has all error handling, you
         * could just call wait() as long as you are only expecting to
         * have one child process at a time.
         */
        pid_t ws = waitpid(pid, &childExitStatus, WNOHANG);
        if (ws == -1) { /* error - handle as you wish */
        }

        if (WIFEXITED(childExitStatus)) /* exit code in childExitStatus */ {
            status = WEXITSTATUS(childExitStatus); /* zero is normal exit */
            /* handle non-zero as you wish */
        } else if (WIFSIGNALED(childExitStatus)) /* killed */ {
        } else if (WIFSTOPPED(childExitStatus)) /* stopped */ {
        }
    }
}

//nearestNeighborMultithread(void)
//Get nearest neighbor for all alignment files in 'dbe'
void nearestNeighborMultithread(void) {
    char ** inFiles = (char **) calloc(dbe + 1, sizeof (char *));
    int i;
    char * name = NULL;
    char * sortedFile = NULL;
    char command[BUFSIZ];
    FILE * neighborResults;

    for (i = 0; i < dbe; i++) { //generate name array
        inFiles[i] = (char *) calloc(strlen(outaln) + strlen(dbElements[i]) + strlen(".sam") + 1, sizeof (char));
        strncpy(inFiles[i], outaln, strlen(outaln));
        strcat(inFiles[i], "/");
        strncat(inFiles[i], dbElements[i], strlen(dbElements[i]) - 6);
        strcat(inFiles[i], ".sam");
        //printf("%s\n", inFiles[i]);
    }

    NeighborName = (char *) calloc(strlen("/unsorted_temp_neighbor.aln") + strlen(outaln) + 1, sizeof (char));
    strncpy(NeighborName, outaln, strlen(outaln));
    name = (char *) calloc(strlen("/unsorted_temp_neighbor.aln") + strlen(outaln) + 1, sizeof (char));
    strncpy(name, outaln, strlen(outaln));
    strcat(name, "/unsorted_temp_neighbor.aln");
    neighborResults = fopen(name, "w+");

    for (i = 0; i < dbe; i++) {
        retrieveTime();
        if(CheckFile(inFiles[i]) == 1) {
            printf("Getting nearest neighbors of file: %s [ %s ]\n", inFiles[i], times);
            fflush(NULL);
            ReadSam(inFiles[i], neighborResults);
        }

        else {
            printf("Could not find SAM file named: %s\nSkipping from neighbor\n", inFiles[i]);

            if(NoAlignment == 1)
            {
                printf("Try realigning the reads to the database\n");
            }

            exit(EXIT_FAILURE);
        }
    }

    sortedFile = (char *) calloc(strlen("/sorted_temp_neighbors.aln") + strlen(outaln) + 1, sizeof (char));
    strncpy(sortedFile, outaln, strlen(outaln));
    strcat(sortedFile, "/sorted_temp_neighbors.aln");

    if (dbe > 1) {
        strcpy(command, "sort -k1,1 -k4,4r -s ");
        strcat(command, name);
        strcat(command, " > ");
        strcat(command, sortedFile);
        printf("%s\n", command);
        fflush(NULL);
        system(command);
    }

    if (name) free(name);
    FreeNames(inFiles);
    fclose(neighborResults);
    if (sortedFile) free(sortedFile); // by Roberto

    //free(inFiles);
}

//MultiNeighborParser(void * voidA)
//Multithreaded parsing of the final neighbor file (sorted_temp_neighbors.aln)
void * MultiNeighborParser(void * voidA) {
    struct ReadChunk * psr = (struct ReadChunk *) voidA;

    while (psr != NULL) {
        //printf("%s\n", psr->data);
        psr->pos = 0;
        psr->iter = 0;
        psr->datalength = strlen(psr->data);
        psr->element = 0;

        psr->taxon = -1;

        psr->as = -1;

        while (psr->pos > -1) {
            if (psr->iter == 1) {
                psr->data += psr->pos;
                psr->taxon = atoi(psr->data);
                psr->data -= psr->pos;
            } else if (psr->iter == 2) {
                psr->data += psr->pos;
                psr->accession = atoi(psr->data);
                psr->data -= psr->pos;
            } else if (psr->iter == 3) {
                psr->data += psr->pos;
                psr->as = atof(psr->data);
                psr->data -= psr->pos;
            } else if (psr->iter == 4) {
                psr->data += psr->pos;
                psr->seqstart = atoi(psr->data);
                psr->data -= psr->pos;
            } else if (psr->iter == 5) {
                psr->data += psr->pos;
                psr->seqlen = atoi(psr->data);
                psr->data -= psr->pos;
            }

            if (psr->pos != -1)
                psr->pos = GetToken(psr->data, '\t', psr->datalength, psr->pos);

            psr->iter++;
        }

        psr = GoThroughData(psr);
        psr = psr->next;
    }
}

//ParseMultiNeighbors(void)
//Calls multithreaded functions
void ParseMultiNeighbors(void) {
    struct ReadChunk * tempStruct = RMaster;
    pthread_t thread_id[threads];
    int i = 0;

    while (tempStruct != NULL) {
        pthread_create(&thread_id[i], NULL, &MultiNeighborParser, tempStruct);
        tempStruct = tempStruct->node;
        i++;
    }

    i = 0;
    tempStruct = RMaster;

    while (tempStruct != NULL) {
        pthread_join(thread_id[i], NULL);
        tempStruct = tempStruct->node;
        i++;
    }
}

void GetFinalFileNames(void) {
    if (InputName || OutputName) {
        printf("Output files iniotialized?\n");
        exit(EXIT_FAILURE);
    }


    InputName = (char *) calloc(strlen(outDir) + strlen("/alignments/sorted_temp_neighbors.aln") + 1, sizeof (char));
    OutputName = (char *) calloc(strlen(outDir) + strlen("/Results/Taxonomy.txt") + 1, sizeof (char));


    if (outDir[strlen(outDir) - 1] == '/') {
        strncpy(InputName, outDir, strlen(outDir) - 1);
        strncpy(OutputName, outDir, strlen(outDir) - 1);
    } else {
        strcpy(InputName, outDir);
        strcpy(OutputName, outDir);
    }

    strcat(InputName, "/alignments/sorted_temp_neighbors.aln");
    strcat(OutputName, "/Results/Taxonomy.txt");

    //printf("%s - %s\n", InputName, OutputName);
}

//MeganFormat(void)
//Converts final file 'Taxonomy.txt' to a megan compatible format
//prints 'read_name\ttaxon_id" to file
void MeganFormat(void)
{
    FILE * handle, * Mout;
    char line[BUFSIZ];
    int token = 0;
    char * ptr;
    char * meganin = (char *) calloc(strlen(outDir) + strlen("/Results/Taxonomy.txt") + 1, sizeof (char));
    char * meganres = (char *) calloc(strlen(outDir) + strlen("/Results/megan.txt") + 1, sizeof (char));

    if (outDir[strlen(outDir) - 1] == '/') {
        strncpy(meganin, outDir, strlen(outDir) - 1);
        strncpy(meganres, outDir, strlen(outDir) - 1);
    } else {
        strcpy(meganin, outDir);
        strcpy(meganres, outDir);
    }

    strcat(meganin, "/Results/Taxonomy.txt");
    strcat(meganres, "/Results/megan.txt");

    if(CheckFile(meganin) == 0) {
        printf("No %s found, can't convert result to megan format\n");
        free(meganin);
        free(meganres);
        return;
    }

    handle = fopen(meganin, "r+");
    Mout = fopen(meganres, "w+");

    while(fgets(line, sizeof(line), handle))
    {
        token = 0;
        ptr = strtok(line, "\t");

        while(ptr != NULL)
        {
            if(token == 0)
                fprintf(Mout, "%s\t", ptr);

            if(token == 1)
                fprintf(Mout, "%s", ptr);

            token++;
            ptr = strtok(NULL, "\t");
        }

        fprintf(Mout,"\n");
    }

    fclose(handle);
    fclose(Mout);

    free(meganin);
    free(meganres);
    return;
}

void CopyFileOneIndex(void)
{
    InputName = (char *) calloc(strlen(outDir) + strlen("/alignments/unsorted_temp_neighbor.aln") + 1, sizeof (char));
    OutputName = (char *) calloc(strlen(outDir) + strlen("/Results/Taxonomy.txt") + 1, sizeof (char));

    if (outDir[strlen(outDir) - 1] == '/') {
        strncpy(InputName, outDir, strlen(outDir) - 1);
        strncpy(OutputName, outDir, strlen(outDir) - 1);
    } else {
        strcpy(InputName, outDir);
        strcpy(OutputName, outDir);
    }

    strcat(InputName, "/alignments/unsorted_temp_neighbor.aln");
    strcat(OutputName, "/Results/Taxonomy.txt");

    //printf("%s\n", InputName);

    if( access( InputName, F_OK ) != -1 ) {
        copyFile(InputName, OutputName);
    } else {
        printf("No nearest nighbor file found\nCould not copy results\nExiting\n");
        exit(EXIT_FAILURE);
    }

    copyFile(InputName, OutputName);
}

//FinalNearestNeighbor(void)
//reads sorted file of taxonomies, copies data into memory, and calls parser + neighbor findwer
void FinalNearestNeighbor(void) {
    FILE * handle, * testout;
    char line[BUFSIZ];
    int i = 0;
    int j = 0;
    int num = 0;
    int EqName = 0;
    char * currName = NULL;
    char * prevName = NULL;

    retrieveTime();
    printf("Final nearest neighbor [ %s ]\n", times);
    fflush(NULL);

    GetFinalFileNames();

    if( access( InputName, F_OK ) != -1 ) {
    // file exists
    } else {
        printf("Could not locate '%d'\nNearest neighbor can not be done\nExiting\n");
        exit(EXIT_FAILURE);
    }

    handle = fopen(InputName, "r+");
    //handle = fopen("sorted.aln", "r+");
    testout = fopen(OutputName, "w+");

    while (fgets(line, sizeof (line), handle)) {
        currName = ReturnName(line, currName);
        //printf("%s\n", currName);

        if (RMaster == NULL) {
            AddReadStructNode(line, EqName);
        } else if (CompareStrings(currName, prevName) == 1) {
            AddReadElement(line, EqName);
        } else {
            EqName++;

            if (i >= chunkSize) {
                i = 0;
                num++;

                if (num == threads) {
                    i = 0;
                    num = 0;
                    EqName = 0;

                    ParseMultiNeighbors();
                    InitializeNearestMulti(RMaster);
                    FreeReadStruct(RMaster, testout);
                    RMaster = NULL;
                }

                AddReadStructNode(line, EqName);
            } else if (RMaster == NULL)
                AddReadStructNode(line, EqName);

            else
                AddReadElement(line, EqName);
        }

        prevName = ReturnName(line, prevName);
        i++;
    }
    ParseMultiNeighbors();
    InitializeNearestMulti(RMaster);
    FreeReadStruct(RMaster, testout);
    fclose(handle);
    fclose(testout); // by Roberto
    if (currName) {
        free(currName);
    }
    if (prevName) {
        free(prevName);
    }
}

int main(int argc, char *argv[]) {
    //FILE * handOut = fopen("ellenorzo.txt", "w+");
    CheckCommands(argv, argc);
    //readTaxonomyFile(taxPath);
    GetMaxNodes(taxPath);

    if (NoAlignment == 0)
        Alignments();

    nearestNeighborMultithread();

    if(dbe > 1) {
        FinalNearestNeighbor();
    }

    else {
        CopyFileOneIndex();
    }

    if(meganOut == 1)
        MeganFormat();

    InputFree();
    //taxonomyFree();
    if (NeighborName) {
        free(NeighborName);
    }

    FreeData();
    //fclose(handOut);
    return 0;
}
