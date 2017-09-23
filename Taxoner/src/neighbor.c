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

struct Reads {
    char * sam;
    char * name;
    int readID;
    float as;
    int taxon;
    int gi;
    int start;
    int length;
    char alnStat[200];
    struct Reads * next;
};

struct Data {
    int element;
    struct Reads * reads;
    struct Data * next;
};

struct Data * head = NULL;
struct Data * curr = NULL;

struct Reads * curr_read = NULL;

int chunkSize = 40000; //number of reads for each thread
int AllTax = 100; //max taxons while multithread reading

/*
 *  Returns read names
 */
char * ReturnName(char * input, char * dest) {
    int j = 0;
    int len = strlen(input) - 1;

    if (dest == NULL)
        free(dest);

    dest = (char *) calloc(len + 1, sizeof (char));

    while (input[j] != '\t' && input[j] != ' ' && input[j] != '\n' && j < len) {
        dest[j] = input[j];
        j++;
    }    //char * prevName = NULL;


    return dest;
}

/*
 *  Checks if SAM entry is aligned (0: true, 1: false)
 */
int CheckIfAligned(char * input) {
    if (strstr(input, "AS:i") != NULL)
        return 0;

    return 1;
}

/*
* ReturnMismatchScore
*/
int ReturnMismatchScore(char * input) {
    char * ptr = strstr(input, "XM:i");

    if(ptr == NULL) return -1;

    ptr += 5;
    return atoi(ptr);
}

/*
 * Returns alignment score
 */
int ReturnAlignmentScore(char * input) {
    char * ptr = strstr(input, "AS:i");

    if(ptr == NULL) return -1;

    ptr += 5;
    return atoi(ptr);
}

/*
 * Returns alignment score needed for taxoner
 */
float GetAlignmentScore(char * input, char * cigar, int seqlentemp, int templen) {
    char * seqlen = NULL;
    int digit = 0;
    int len = strlen(cigar);
    int i, j, k;
    char temp;
    int total = 0;
    int tempscore = ReturnAlignmentScore(input);
    float score = 0;

    if(parms->bowtie2_local == 0)
        return ((float)templen + (float)tempscore)/(float)templen;

    else {
        seqlen = (char *)calloc(11, sizeof(char));
        for(i = 0; i < len; i++) {
            if(isalpha(cigar[i])) {
                if(cigar[i] == 'M' || cigar[i] == 'I')
                    total += atoi(seqlen);

                 if(seqlen)
                    free(seqlen);

                seqlen = (char *)calloc(11, sizeof(char));
                digit = 0;

            }

            else {
                seqlen[digit] = cigar[i];
                digit++;
            }
        }

        if(seqlen)
            free(seqlen);

        score = (float)tempscore / ((float)total * 2);
        return score;
    }
}

/*
 *  Allocates head pointer to store data for each thread separately
 */
void AllocateHead(void) {
    head = (struct Data *)calloc(1, sizeof(struct Data));
    curr = head;

    head->next = NULL;
    head->reads = NULL;
    head->element = 0;
    curr_read = curr->reads;
}

/*
 *  Allocated data structs for next thread
 */
void AllocateData(void) {
    struct Data * ptr = NULL;
    if(head == NULL) {
        AllocateHead();
        return;
    }

    ptr = (struct Data *)calloc(1, sizeof(struct Data));
    ptr->next = NULL;
    ptr->reads = NULL;
    ptr->element = curr->element + 1;

    curr->next = ptr;
    curr = curr->next;
    curr_read = curr->reads;
}

/*
 *  Compares to input strings by length and chars
 */
int CompareStrings(char * s1, char * s2){
    if(strlen(s1) == strlen(s2)) {
        if(strncmp(s1, s2, strlen(s1)) == 0)
            return 1;
    }

    return 0;
}

/*
 *  Adds read info to Reads struct linked list
 */
void AddRead(char * source, char * name) {
    struct Reads * ptr = (struct Reads *)calloc(1, sizeof(struct Reads));
    ptr->next = NULL;
    ptr->sam = NULL;
    ptr->readID = 0;
    ptr->name = NULL;
    ptr->length = 0;
    ptr->start = 0;
    ptr->gi = 0;
    ptr->as = 0;
    ptr->taxon = 0;

    ptr->sam = (char *)calloc(strlen(source) + 1, sizeof(char));
    strcpy(ptr->sam, source);

    ptr->name = (char *)calloc(strlen(name) + 1, sizeof(char));
    strcpy(ptr->name, name);

    if(curr->reads == NULL) {
        curr->reads = ptr;
        curr_read = ptr;
        return;
    }

    curr_read->next = ptr;
    ptr->readID = curr_read->readID;

    if(CompareStrings(ptr->name, curr_read->name) == 0)
        ptr->readID++;

    curr_read = curr_read->next;
}

/*
 *  Frees allocated thread structs
 */
void FreeStruct(void) {
    struct Reads * ptr;
    while(head != NULL) {
        curr = head;
        head = head->next;

        curr_read = curr->reads;
        while(curr_read != NULL) {
            ptr = curr_read;
            curr_read = curr_read->next;
            if(ptr->sam) free(ptr->sam);
            if(ptr->name) free(ptr->name);
            free(ptr);
        }

        free(curr);
    }

    head = NULL;
    curr = NULL;
    curr_read = NULL;
}

/*isVirus: Check if taxon id is a virus or other sequences, and exclude from results*/
/*
 *  Checks if taxon ID is blacklisted
 */
int Blacklist(int a) {
    int temptax = a;
    int l1 = 0;

    if (a == -1)
        return 1;

    while (nodes[temptax].curr != nodes[temptax].next && l1 < AllTax) {
        if (nodes[temptax].ok == 0)
            return 1;
        temptax = nodes[temptax].next;
        l1++;
    }

    return 0;
}

/*
 *  Deletes imported read
 */
void FreeReadElement(struct Reads * ptr) {
    if(ptr->name) free(ptr->name);
    if(ptr->sam) free(ptr->sam);

    free(ptr);
    ptr = NULL;
}
/*
 * ReturnNeighbor(int a, int b)
 * looks for nearest neighbor for taxons 'a' and 'b', and returns
 * nearest neighbor id (integer)
 */
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
    while (nodes[temptax].curr != nodes[temptax].next && l2 < AllTax) { //store taxID lineage in s2
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

/*
 *  Parse stored reads in a multithreaded-safe way
 */
void ParseRead(struct Reads * psr) {
    int token = 0;
    char * ptr;
    char * semicolon;
    char *storage;
    char * cigar = NULL; //lorinc
    char * temp = (char *)calloc(strlen(psr->sam) + 1, sizeof(char));
    strcpy(temp, psr->sam);

    ptr = strtok_r( temp, "\t", &storage);
    do
    {
        if(token == 2) {
            psr->gi = atoi(ptr);

            semicolon = strstr(ptr, ";");
            if(semicolon != NULL) {
                semicolon += 1;
                psr->taxon = atoi(semicolon);
            }
        }

        if(token == 3) psr->start = atoi(ptr);
        if(token == 5) {//lorinc
            cigar = (char *)calloc(strlen(ptr) + 1, sizeof(char));//lorinc
            strncpy(cigar, ptr, strlen(ptr));//lorinc
        }//lorinc
        if(token == 9) psr->length = strlen(ptr);

        token++;
    }
    while( ptr = strtok_r( NULL, "\t", &storage));

    if(parms->alnstats == 1) {
        strncat(psr->alnStat, cigar, strlen(cigar));
        strcat(psr->alnStat, ",");
        strcat(psr->alnStat, "XM:");
        sprintf(psr->alnStat + strlen(psr->alnStat), "%d", ReturnMismatchScore(psr->sam));
        strcat(psr->alnStat, ",");
        strcat(psr->alnStat, "AS:");
        sprintf(psr->alnStat + strlen(psr->alnStat), "%d", ReturnAlignmentScore(psr->sam));
        strcat(psr->alnStat, ",");
    }

    psr->as = GetAlignmentScore(psr->sam, cigar, psr->length, psr->length);
    free(temp);
    if(cigar)//lorinc
        free(cigar);//lorinc
}

/*
 *  Starts read parsing and checks if reads have to be filtered
 */
void ParseAndFilterReads(struct Data * ptr) {
    struct Reads * psr = ptr->reads;
    struct Reads * temp = ptr->reads;

    while(psr != NULL) {
        ParseRead(psr);

        if(Blacklist(psr->taxon) == 1) {
            if(temp == ptr->reads) {
                ptr->reads = psr->next;
                temp = ptr->reads;
            }

            else {
                temp = psr->next;
            }

            free(psr);
            psr = temp;
        }

        else
            psr = psr->next;
    }
}

/*
 *  Parse data in coming from .aln file
 */
void ParseReadFinal(struct Reads * psr) {
    int token = 0;
    char * ptr;
    char * semicolon;
    char *storage;
    char * temp = (char *)calloc(strlen(psr->sam) + 1, sizeof(char));
    strcpy(temp, psr->sam);

    ptr = strtok_r( temp, "\t", &storage);
    do
    {
        if(token == 1) psr->taxon = atoi(ptr);
        if(token == 2) psr->gi = atoi(ptr);
        if(token == 3) psr->as = atof(ptr);
        if(token == 4) psr->start = atoi(ptr);
        if(token == 5) psr->length = atoi(ptr) - psr->start;
        if(token == 6) strcpy(psr->alnStat, ptr);
        token++;
    }
    while( ptr = strtok_r( NULL, "\t", &storage));

    free(temp);
}

/*
 *  Starts read parsing and checks if reads have to be filtered from .aln file
 */
void ParseAndFilterFinalReads(struct Data * ptr) {
    struct Reads * psr = ptr->reads;
    struct Reads * temp = ptr->reads;

    while(psr != NULL) {
        ParseReadFinal(psr);

        if(Blacklist(psr->taxon) == 1) {
            if(temp == ptr->reads) {
                ptr->reads = psr->next;
                temp = ptr->reads;
            }

            else {
                temp = psr->next;
            }

            free(psr);
            psr = temp;
        }

        else
            psr = psr->next;
    }
}

/*
 *  calculate nearest neighbor (called by each thread)
 */
void * NearestNeighbor(void * voidA) {
    struct Data * tData = (struct Data *) voidA;
    struct Reads * psr = tData->reads;
    struct Reads * temp = psr;

    while(psr != NULL) {
        temp = psr->next;
        if(temp == NULL) return NULL;
        if(CompareStrings(temp->name, psr->name) == 1) {
            if(temp->as > psr->as) {
                if(temp->as * parms->neighborScore < psr->as)
                    temp->taxon = ReturnNeighbor(temp->taxon, psr->taxon);

                psr->taxon = temp->taxon;
                psr->as = temp->as;
                psr->gi = temp->gi;
                psr->length = temp->length;
                psr->start = temp->start;
            }

            else {
                if(psr->as * parms->neighborScore < temp->as)
                    psr->taxon = ReturnNeighbor(psr->taxon, temp->taxon);

                temp->taxon = psr->taxon;
                temp->as = psr->as;
                temp->gi = psr->gi;
                temp->length = psr->length;
                temp->start = psr->start;
            }
        }

        psr = psr->next;
    }
    return NULL;
}

/*
 *  Print neighbor results to file
 */
void PrintResults(FILE * output) {
    struct Data * psr = head;
    struct Reads * ptr;
    int i;
    while(psr != NULL) {
        ptr = psr->reads;

       while(ptr != NULL) {
            if(ptr->next != NULL) {
                if(CompareStrings(ptr->next->name, ptr->name) == 0) {
                   fprintf(output, "%s\t%d\t%d\t%.3f\t%d\t%d", ptr->name, ptr->taxon, ptr->gi, ptr->as, ptr->start, (ptr->start + ptr->length));

                    if(parms->alnstats == 1) {
                        if(strlen(ptr->alnStat) > 1) {
                            fprintf(output, "\t");
                            for(i = 0; i < strlen(ptr->alnStat)-1; i++) {
                                if(ptr->alnStat[i] != ' ')
                                    fprintf(output, "%c", ptr->alnStat[i]);
                            }
                        }
                        //printf("%s\n", ptr->alnStat);
                        else {
                            fprintf(output, "\t-1");
                        }
                    }

                    fprintf(output, "\n");
                }
            }

            else {
                fprintf(output, "%s\t%d\t%d\t%.3f\t%d\t%d", ptr->name, ptr->taxon, ptr->gi, ptr->as, ptr->start, (ptr->start + ptr->length));

                if(parms->alnstats == 1) {
                    if(strlen(ptr->alnStat) > 1) {
                        fprintf(output, "\t");
                        for(i = 0; i < strlen(ptr->alnStat)-1; i++) {
                            if(ptr->alnStat[i] != ' ')
                                fprintf(output, "%c", ptr->alnStat[i]);
                        }
                    }
                        //printf("%s\n", ptr->alnStat);
                    else {
                        fprintf(output, "\t-1");
                    }
                }

                fprintf(output, "\n");
            }

            ptr = ptr->next;
        }

        psr = psr->next;
    }
}

/*
 *  start neighbor calculation
 */
void * ParseAndFindNeighbor(void * voidA) {
    struct Data * tData = (struct Data *) voidA;
    ParseAndFilterReads(tData);
    NearestNeighbor(tData);
}

/*
 *  call multithreaded nearest neighbor algorithm
 */
void InitializeNearestMulti(struct Data * psr) {
    struct Data * tempStruct = psr;
    pthread_t thread_id[parms->threads];
    int i = 0;

    while (tempStruct != NULL) { //calls threads
        pthread_create(&thread_id[i], NULL, &ParseAndFindNeighbor, tempStruct);
        tempStruct = tempStruct->next;
        i++;
    }

    i = 0;
    tempStruct = psr;

    while (tempStruct != NULL) { //joins threads
        pthread_join(thread_id[i], NULL);
        tempStruct = tempStruct->next;
        i++;
    }
}

/*
 *  Start nearest neighbor calculation
 */
void * ParseAndFindFinalNeighbor(void * voidA) {
    struct Data * tData = (struct Data *) voidA;
    ParseAndFilterFinalReads(tData);
    NearestNeighbor(tData);
}

/*
 *  call multithreaded nearest neighbor algorithm on sorted .aln file
 */
void InitializeFinalNearestMulti(struct Data * psr) {
    struct Data * tempStruct = psr;
    pthread_t thread_id[parms->threads];
    int i = 0;

    while (tempStruct != NULL) { //calls threads
        pthread_create(&thread_id[i], NULL, &ParseAndFindFinalNeighbor, tempStruct);
        tempStruct = tempStruct->next;
        i++;
    }

    i = 0;
    tempStruct = psr;

    while (tempStruct != NULL) { //joins threads
        pthread_join(thread_id[i], NULL);
        tempStruct = tempStruct->next;
        i++;
    }
}

/*
 *  import .aln file into memory
 */
void ReadAlnFile(char * infile, char * outfile) {
    FILE * handler = fopen(infile, "r");
    FILE * out = fopen(outfile, "w+");
    char line[BUFSIZ];
    char * currName = NULL;

    if(head != NULL) FreeStruct();
    AllocateData();

    while (fgets(line, sizeof (line), handler)) {
        currName = ReturnName(line, currName);
        AddRead(line, currName);

        if(curr_read != NULL) {
            if(curr_read->readID >= chunkSize) {
                if(curr->element == parms->threads) {
                    curr_read->readID = -1;
                    //NearestNeighbor(head);
                    InitializeFinalNearestMulti(head);
                    PrintResults(out);
                    FreeStruct();
                }

                AllocateData();
                AddRead(line, currName);
            }
        }

        if(currName) free(currName);
        currName = NULL;
    }

    if(head != NULL) {
        //NearestNeighbor(head);
        InitializeFinalNearestMulti(head);
        PrintResults(out);
    }

    FreeStruct();
    fclose(handler);
    fclose(out);
}

/*
 *  import sam file data
 */
void ReadSamFile(char * infile, FILE * out) {
    FILE * handler = fopen(infile, "r");
    char line[BUFSIZ];
    char * currName = NULL;

    if(head != NULL) FreeStruct();
    AllocateData();

    while (fgets(line, sizeof (line), handler)) {
        if (line[0] != '@' && CheckIfAligned(line) == 0) {
            currName = ReturnName(line, currName);
            AddRead(line, currName);

            if(curr_read != NULL) {
                if(curr_read->readID >= chunkSize) {
                    if(curr->element == parms->threads) {
                        curr_read->readID = -1;
                        //NearestNeighbor(head);
                        InitializeNearestMulti(head);
                        PrintResults(out);
                        FreeStruct();
                    }

                    AllocateData();
                    AddRead(line, currName);
                }
            }
            if(currName) free(currName);
            currName = NULL;
        }
    }

    if(head != NULL) {
        //NearestNeighbor(head);
        InitializeNearestMulti(head);
        PrintResults(out);
    }
    FreeStruct();
    fclose(handler);
}

void StartNeighbor(void) {
    int i = 0;
    char * inFile = NULL;
    char * outFile = NULL;
    FILE * out;
    char command[BUFSIZ];

    outFile = (char *)calloc(strlen(parms->alignOut) + strlen("unsorted.aln") + 2, sizeof(char));
    sprintf(outFile, "%s/unsorted.aln", parms->alignOut);
    out = fopen(outFile, "w+");

    for(i = 0; i < db->indexes; i++) {
        if(inFile) free(inFile);
        inFile = (char *)calloc(strlen(db->names[i]) + strlen(parms->alignOut) + 6, sizeof(char));
        sprintf(inFile, "%s/%s.sam", parms->alignOut, db->names[i]);

        if(parms->verbose == 1) printf("Analyzing: %s\n", inFile);
        ReadSamFile(inFile, out);

        if(inFile) free(inFile);
        inFile = NULL;
    }

    fclose(out);
    sprintf(command, "sort -k1,1 -k4,4r -s %s > %s/sorted.aln", outFile, parms->alignOut);
    if(parms->verbose == 1) printf("Sorting results in %s\n Cmd: %s\n", outFile, command);
    system(command);

    free(outFile);
    free(inFile);

    outFile = (char *)calloc(strlen(parms->output) + strlen("Results/Taxonomy.txt") + 2, sizeof(char));
    sprintf(outFile, "%s/Results/Taxonomy.txt", parms->output);

    inFile = (char *)calloc(strlen(parms->alignOut) + strlen("/sorted.aln") + 1, sizeof(char));
    sprintf(inFile, "%s/sorted.aln", parms->alignOut);

    if(parms->verbose == 1) printf("Getting final neighbors from: %s\n", inFile);
    ReadAlnFile(inFile, outFile);

    if(inFile) free(inFile);
    if(outFile) free(outFile);
}

/*
 *  convert results to megan format
 */
void PrintMeganFormat(void) {
    char * output;
    FILE * out;
    FILE * in;
    char * input;
    char line[BUFSIZ];
    char * ptr;
    int token = 0;

    if(parms->verbose == 1) printf("Creating megan output\n");

    input = (char *)calloc(strlen(parms->output) + strlen("/Results/Taxonomy.txt") + 1, sizeof(char));
    output = (char *)calloc(strlen(parms->output) + strlen("/Results/megan.txt") + 1, sizeof(char));

    sprintf(input, "%s/Results/Taxonomy.txt", parms->output);
    sprintf(output, "%s/Results/megan.txt", parms->output);

    out = fopen(output, "w+");
    in = fopen(input, "r+");

    while (fgets(line, sizeof (line), in)) {
        ptr = strtok(line, "\t");
        token = 0;
        while(ptr != NULL) {
            if(token == 0) fprintf(out, "%s\t", ptr);
            if(token == 1) fprintf(out, "%s\t", ptr);
            if(token == 3) fprintf(out, "%.f\n", atof(ptr) * 100);
            token++;
            ptr = strtok(NULL, "\t");
        }
    }

    free(input);
    free(output);
    fclose(in);
    fclose(out);
}

void GetNeighbors(void) {
    if(parms->verbose == 1) printf("\n***************************************\n");
    if(parms->verbose == 1) printf("\tNearest neighbor");
    if(parms->verbose == 1) printf("\n***************************************\n");
    StartNeighbor();
    if(parms->meganOutput == 1) PrintMeganFormat();
    if(parms->verbose == 1) printf("\n***************************************\n");
}
