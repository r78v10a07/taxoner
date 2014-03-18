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

//another comment

struct Taxons {
    int id;
    char * path;
    struct Taxons * next;
};

struct ReadChunk {
    char * data;
    int datalength;
    int accession;
    int seqstart;
    int seqlen;
    float as;
    int taxon;
    int iter;
    int element;
    int pos;
    char * chel;
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

int TotalSams = 0;
int chunkSize = 40000;
int AllTax = 100; //max taxons while multithread reading


int tempnode;

char * NeighborName = NULL;

struct Taxons * thead = NULL;
struct Taxons * tcurr = NULL;

struct ReadChunk * RMaster = NULL;
struct ReadChunk * Rhead = NULL;
struct ReadChunk * Rcurr = NULL;

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

int isVirus(int a) {
    int temptax = a;
    int l1 = 0;
    int l2 = 0;

    if (a == -1)
        return 1;

    //printf("%d\n", a);

    while (n[temptax].id != n[temptax].nextId && l1 < AllTax) {
        if (temptax == 10239 || temptax == 28384)
            return 1;
        //printf("%d;", s1[l1]);
        temptax = n[temptax].nextId;
        l1++;
    }

    return 0;
}

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

void PrintTaxons(void) {
    while (thead != NULL) {
        tcurr = thead;
        thead = tcurr->next;
        //printf("%d\t%s\n", tcurr->id, tcurr->path);

        if (tcurr->path)free(tcurr->path);
        if (tcurr)free(tcurr);
    }
}

/*void CopyTaxons(void)
{
    printf("Filling taxon id(s) and path(s) into memory\n");
    while(thead != NULL)
    {
        tcurr = thead;
        thead = tcurr->next;

        if(tcurr->id > 0 && tcurr->id < maxTax)
        {
            tax[tcurr->id].path = (char *)calloc(strlen(tcurr->path) + 1, sizeof(char));
            strncpy(tax[tcurr->id].path, tcurr->path, strlen(tcurr->path));
        }

        else
        {
            printf("There was a problem with taxon id: %d (path: %s)\n", tcurr->id, tcurr->path);
            printf("Maybe it is smaller than 0 or bigger than %d?\n", maxTax-1);
        }

        free(tcurr->path);
        free(tcurr);
    }
}*/

/*void Allocatetaxons(int num)
{
    int i = 0;
    tax = (rint *)calloc(num + 1, sizeof(rint));

    if(tax == NULL)
    {
        printf("Could not allocate struct\nExiting\n");
        exit(EXIT_FAILURE);
    }

    num++;

    for(i = 0; i < num; i++)
    {
        tax[i].path = NULL;
    }

    CopyTaxons();
}*/

/*void readTaxonomyFile(char * source)
{
    FILE * handle = fopen(source, "r+");
    char line[BUFSIZ];
    int token;
    char * ptr, * input;
    int taxa;
    int count = 0;

    printf("\n--------------------------------------\n");
    printf("Taxonomy info\n");
    printf("\nReading taxonomy file (%s)\n", source);

    while(fgets(line, sizeof(line), handle))
    {
        input = CopyString(line, strlen(line) - 1);
        ptr = strtok(input, "\t");
        token = 0;
        count++;
        //printf("%s\n", line);
        while(ptr != NULL)
        {
            if(token == 0)
                taxa = atoi(ptr);

            if(token == 1)
            {
                if(thead == NULL)
                    createTaxons(taxa, ptr);

                else
                    AddTaxons(taxa, ptr);

                if(taxa > maxTax)
                    maxTax = taxa;
            }

            token++;
            ptr = strtok(NULL, "\t");
        }

        free(input);
    }

    maxTax++;

    fclose(handle);
    printf("Found %d taxons, %d was the largest id number\n", count, maxTax);
    //PrintTaxons();
    Allocatetaxons(maxTax + 1);
    printf("--------------------------------------\n\n");
}*/

int databaseAlignment(int num, char * infile) {
    retrieveTime();
    char * partailfulldb;
    char * partialdb;
    char command[BUFSIZ];
    partailfulldb = (char *) calloc(strlen(fulldb[num]), sizeof (char));
    strncpy(partailfulldb, fulldb[num], strlen(fulldb[num]) - 6);

    partialdb = (char *) calloc(strlen(dbElements[num]), sizeof (char));
    strncpy(partialdb, dbElements[num], strlen(dbElements[num]) - 6);

    //printf("databaseAlignment: [%s] [%s] [%s]\n", infile, outaln, partialdb);
    if (allHits == 1) {
        if (filterHost == 1 && host != NULL) {
            //printf("1\n");
            sprintf(command, "bowtie2 -p %d -a %s %s/%s > %s/%s.sam", threads, partailfulldb, outaln, infile, outaln, partialdb);
        } else {
            //printf("2\n");
            sprintf(command, "bowtie2 -p %d -a %s %s > %s/%s.sam", threads, partailfulldb, infile, outaln, partialdb);
        }
    } else {
        if (filterHost == 1 && host != NULL) {
            //printf("3\n");
            sprintf(command, "bowtie2 -p %d -k %d %s %s > %s/%s.sam", threads, bowtieMaxHits, partailfulldb, outaln, infile, outaln, partialdb);
        } else {
            //printf("4\n");
            sprintf(command, "bowtie2 -p %d -k %d %s %s > %s/%s.sam", threads, bowtieMaxHits, partailfulldb, infile, outaln, partialdb);
        }
    }


    printf("Command: %s\n", command);
    fflush(NULL);
    //system(command);
    printf("Index: %s\nReads: %s\nStart time: %s\n", partialdb, infile, times);
    fflush(NULL);
    //printf("Aligning %s to database no. %d (%s)\n", infile, num, dbElements[num]);
    printf("Bowtie2 output:\n");
    fflush(NULL);
    system(command);


    printf("Getting nearest neighbors from %s database alignment\n", partailfulldb);
    fflush(NULL);
    retrieveTime();
    printf("Finished at: %s\n\n", times);
    fflush(NULL);

    if (partailfulldb) free(partailfulldb);
    if (partialdb) free(partialdb);

    return 0;
}

int hostAlignment(void) {
    char command[BUFSIZ];
    retrieveTime();
    printf("--------------------------------------\n");
    printf("Host alignment [ %s ]\n", times);
    printf("--------------------------------------\n");
    printf("Index: %s\nReads: %s\n", host, reads);
    fflush(NULL);
    
    sprintf(command, "bowtie2 -p %d --un %s/%s %s %s > %s/host.sam", threads, outaln, nohost, host, reads, outaln);
    printf("Command: %s\n\n", command);
    printf("Bowtie2 output:\n");
    system(command);
    printf("--------------------------------------\n\n");
    fflush(NULL);
    return 0;
}

void printPreDatabaseInfo(void) {
    char command[BUFSIZ];
    retrieveTime();
    printf("--------------------------------------\n");
    printf("Db alignment [ %s ]\n", times);
    printf("--------------------------------------\n");
    fflush(NULL);
}

void Alignments(void) {
    int i = 0;

    if (filterHost == 1 && host != NULL) {
        hostAlignment();
        printPreDatabaseInfo();
        for (i = 0; i < dbe; i++) {
            databaseAlignment(i, nohost);
        }

        return;
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

            //printf("%s\n", temp->data);
            //if(temp->data)
            //free(temp->data);

            tempOut = ReturnName(temp->data, tempOut);

            if (temp->as > 0)
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
    //psr = NULL;
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

    //printf("\nold: ");

    while (n[temptax].id != n[temptax].nextId && l1 < AllTax) {
        s1[l1] = temptax;
        //printf("%d;", s1[l1]);
        temptax = n[temptax].nextId;
        l1++;
    }


    //printf("\nnew: ");
    temptax = b;

    while (n[temptax].id != n[temptax].nextId && l2 < AllTax) {
        s2[l2] = temptax;
        //printf("%d;", s2[l2]);
        temptax = n[temptax].nextId;
        l2++;
    }

    //printf("\n\n");

    if (l1 > l2) {
        j = l1 - 1;

        for (i = l2 - 1; i > -1; i--) {
            //printf("%d - %d\n", s2[i], s1[j]);
            if (s2[i] == s1[j])
                temptax = s2[i];

            else
                i = -1;

            j--;
        }
    } else {
        j = l2 - 1;

        for (i = l1 - 1; i > -1; i--) {
            //printf("%d - %d\n", s2[j], s1[i]);
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

int CheckIfAligned(char * input) {
    if (strstr(input, "AS:i") != NULL)
        return 0;

    return 1;
}

struct ReadChunk * FreeNextBranch(struct ReadChunk * psr) {
    struct ReadChunk * tempStruct = psr->next;
    psr->next = NULL;
    psr->next = tempStruct->next;

    if (tempStruct->data) free(tempStruct->data);
    if (tempStruct) free(tempStruct);

    return psr;
}

void * ParallelNeighbor(void * voidA) {
    struct ReadChunk * psr = (struct ReadChunk *) voidA;
    int asd = 0;
    psr = RemoveVirus(psr);
    while (psr != NULL) {
        if (psr->next != NULL && psr->iter > 3) {
            if (psr->next->nameElement == psr->nameElement && psr->as != -1) {
                if (psr->next->next != NULL) {
                    //printf("%s\n%s\n (%.2f - %.2f)\n", psr->data, psr->next->data, psr->as, psr->next->as);
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

void CPUParallelCalc(struct ReadChunk * psr) {
    int i, j;
    char * tempplace, * toks;
    struct ReadChunk * asd = psr;

    while (asd != NULL) {
        psr = asd;
        asd = asd->node;

        while (psr != NULL) {
            tempplace = CopyString(psr->data, strlen(psr->data));
            toks = strtok(tempplace, "\t");

            while (toks != NULL) {
                for (i = 0; i < 100; i++) {

                }

                toks = strtok(NULL, "\t");
            }

            //if(tempplace)
            //free(tempplace);

            psr = psr->next;
        }
    }
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

void InitializeNearestMulti(struct ReadChunk * psr) {
    struct ReadChunk * tempStruct = psr;
    pthread_t thread_id[threads];
    int i = 0;

    while (tempStruct != NULL) {
        pthread_create(&thread_id[i], NULL, &ParallelNeighbor, tempStruct);
        tempStruct = tempStruct->node;
        i++;
    }

    i = 0;
    tempStruct = psr;

    while (tempStruct != NULL) {
        pthread_join(thread_id[i], NULL);
        tempStruct = tempStruct->node;
        i++;
    }
}

void ParseDataMulti(void) {
    struct ReadChunk * tempStruct = RMaster;
    pthread_t thread_id[threads];
    int i = 0;

    while (tempStruct != NULL) {
        pthread_create(&thread_id[i], NULL, &MultiParser, tempStruct);
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

void ReadSam(char * infile, FILE * out) {
    //FILE * handler = fopen("res/alignments/few.sam", "r+");
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
                currName = ReturnName(line, currName);

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

                            ParseDataMulti();
                            InitializeNearestMulti(RMaster);
                            FreeReadStruct(RMaster, out);
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

void nearestNeighborMultithread(void) {
    char ** inFiles = (char **) calloc(dbe + 1, sizeof (char *));
    int i;
    char * name = NULL;
    char * sortedFile = NULL;
    char command[BUFSIZ];
    FILE * neighborResults;

    for (i = 0; i < dbe; i++) {
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
        printf("Getting nearest neighbors of file: %s [ %s ]\n", inFiles[i], times);
        fflush(NULL);
        ReadSam(inFiles[i], neighborResults);
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

    //copyFile(name, finalname);
    //unlink(name);
    if (name) free(name);
    FreeNames(inFiles);
    fclose(neighborResults);
    if (sortedFile) free(sortedFile); // by Roberto

    //free(inFiles);
}

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

/*void InitializeNearestMulti(struct ReadChunk * psr)
{
    struct ReadChunk * tempStruct = psr;
    pthread_t thread_id[threads];
    int i = 0;

    while(tempStruct != NULL)
    {
        pthread_create(&thread_id[i], NULL, &ParallelNeighbor, tempStruct);
        tempStruct = tempStruct->node;
        i++;
    }

    i = 0;
    tempStruct = psr;

    while(tempStruct != NULL)
    {
        pthread_join(thread_id[i], NULL);
        tempStruct = tempStruct->node;
        i++;
    }
}*/

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
    FinalNearestNeighbor();
    InputFree();
    //taxonomyFree();
    if (NeighborName) {
        free(NeighborName);
    }

    FreeData();
    //fclose(handOut);
    return 0;
}
