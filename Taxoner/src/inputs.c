#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <getopt.h>
#include <time.h>
#include <stdbool.h>
#include <dirent.h>
#include <sys/stat.h>
#include <pthread.h>
#include <unistd.h>
#include "main.h"
#include <ctype.h>

extern struct TaxonNodes * nodes;
extern char * program_name;
extern struct params * parms;
extern struct database * db;

/*
 *  Free allocated data stored in params struct
 */
void FreeInputs(void) {
    int i = 0;

    if(parms != NULL) {
        if(parms->database) free(parms->database);
        if(parms->names) free(parms->names);
        if(parms->input) free(parms->input);
        if(parms->output) free(parms->output);
        if(parms->bowtieExec) free(parms->bowtieExec);
        if(parms->bwt2Indexes) free(parms->bwt2Indexes);
        if(parms->pairedReads) free(parms->pairedReads);
        if(parms->alignOut) free(parms->alignOut);
        if(parms->hostindex) free(parms->hostindex);
        free(parms);
    }

    if(db != NULL) {
        for(i = 0; i < db->indexes; i++) {
            if(db->names[i]) free(db->names[i]);
            if(db->paths[i]) free(db->paths[i]);
        }

        if(db->names) free(db->names);
        if(db->paths) free(db->paths);

        free(db);
    }

    if(nodes) free(nodes);
}

/*
 *  Prints the taxoner usage
 */
void print_usage(FILE *stream, int exit_code) {
    fprintf(stream, "\n********************************************************************************\n");
    fprintf(stream, "\nUsage: taxoner -d <folder> -s <input reads> -o <output folder> -n <nodes.dmp>\n");
    fprintf(stream, "\nTaxoner64 version 0.1.3\n");
    fprintf(stream, "\n%s options:\n\n", program_name);
    fprintf(stream, " -d <string>,\t--dbPath\tSpecifies folder path to database\n");
    fprintf(stream, " -s <string>,\t--seq\t\tInput fastq/fasta reads\n");
    fprintf(stream, " -o <string>,\t--output\tOutput folder path\n");
    fprintf(stream, " -n <string>,\t--taxpath\tSpecifies file nodes file (nodes.dmp from NCBI)\n\n");
    fprintf(stream, " Database options:\n");
    fprintf(stream, " -l,\t\t--largeGenome\tBowtie2 will use large (>4Gb) genomes\n");
    fprintf(stream, " -g,\t\t--both-genome-sizes\tBowtie2 will use large (>4Gb) and smaller genomes\n");
    fprintf(stream, " -i <string>,\t--bwt2-indexes\tSpecify comma separated names of indexes to be used at\n");
    fprintf(stream, "\t\t\talignment. Indexes present in --dbPath <folder>\n\n");
    fprintf(stream, " Read options:\n");
    fprintf(stream, " -f,\t\t--fasta\t\tInput reads are in fasta format\n");
    fprintf(stream, " -p <string>,\t--paired\tSpecify second pair of reads for paired-end mode\n");
    fprintf(stream, " -I <int>,\t--minInsert\tMinimum insert size for paired-end mode. Default: %d\n", parms->minInsert);
    fprintf(stream, " -X <int>,\t--maxInsert\tMaximum insert size for paired-end mode. Default: %d\n\n", parms->maxInsert);
    fprintf(stream, " Alignment options:\n");
    fprintf(stream, " -b <int>,\t--bt2-maxhits\tMaximum alignments to report. Default: %d\n", parms->bowtieMaxHits);
    fprintf(stream, " -a,\t\t--bt2-allhits\tReport all alignments\n");
    fprintf(stream, " -w <string>,\t--bowtie2\tSpecifies bowtie2 executable\n");
    fprintf(stream, " -c <string>,\t--host-filter\tSpecify bowtie2 index to host genome to filter with.\n");
    fprintf(stream, " -A,\t--alignStats\tAdd an extra CIGAR style column about alignments\n");
    fprintf(stream, " -y <string>,\t--bwt2-params\tInput file with all extra bowtie2 parameters.\n\n");
    fprintf(stream, " Nearest neighbor options:\n");
    fprintf(stream, " -r <float>,\t--neighbor-score\tSpecifies alignment score for nearest neighbor.\n");
    fprintf(stream, "\t\t\t\tDefault score: %.1f\n", parms->neighborScore);
    fprintf(stream, " -e,\t\t--only-neighbor\tPerform only nearest neighbor analysis without alignment\n");
    fprintf(stream, " -k <string>,\t--skip\t\tSpecify file containing taxon IDs to skip from analysis\n");
    fprintf(stream, " -m,\t\t--megan\t\tCreate megan compatible output\n");
    fprintf(stream, " -u,\t\t--virus-filter\tDon't skip viruses from analysis \n");
    fprintf(stream, "\t\t\tIf specified, IDs 10239, 28384 will be used in analysis.\n\n");

    fprintf(stream, " Performance options:\n");
    fprintf(stream, " -t <int>,\t--threads\tNumber of CPU threads to use. Default: %d\n\n", parms->threads);

    fprintf(stream, " -h,\t\t--help\t\tDisplay this usage information.\n");
    fprintf(stream, " -v,\t\t--verbose\tQuiet mode (no output in terminal\n\n");
    fprintf(stream, "********************************************************************************\n");
    fprintf(stream, "\n            Lorinc Pongor (e-mail: pongorlorinc@gmail.com)");
    fprintf(stream, "\n            Roberto Vera Alvarez (e-mail: r78v10a07@gmail.com)\n\n");
    fprintf(stream, "********************************************************************************\n");
    FreeInputs();
    exit(0);
}

/*
 *  Creates a params struct to store argument data. Returns a pointer of struct params *
 */
struct params * CreateInputs(void) {
    struct params * ptr = (struct params *)calloc(1, sizeof(struct params));

    if(ptr == NULL) {
        printf("Unable to allocate parameter structure\nExiting\n");
        FreeInputs();
        exit(EXIT_FAILURE);
    }
    ptr->verbose = 1;
    ptr->database = NULL;
    ptr->names = NULL;
    ptr->hostFilter = 1;
    ptr->input = NULL;
    ptr->output = NULL;
    ptr->bowtieAllHits = 0;
    ptr->bowtieMaxHits = 10;
    ptr->neighborOnly = 0;
    ptr->bwt2Indexes = NULL;
    ptr->fastaInput = 0;
    ptr->meganOutput = 0;
    ptr->minInsert = 0;
    ptr->maxInsert = 500;
    ptr->pairedReads = NULL;
    ptr->bowtieExec = NULL;
    ptr->threads = 2;
    ptr->largeGenome = 0;
    ptr->neighborScore = 0.99;
    ptr->alignOut = NULL;
    ptr->filterVirus = 1;
    ptr->removeTaxons = NULL;
    ptr->results = NULL;
    ptr->hostindex = NULL;
    ptr->bowtie2_params = NULL;
    ptr->bowtie2_params_file = NULL;
    ptr->bowtie2_local = 0;
    ptr->alnstats = 0;
    return ptr;
}

/*
 *  Checks if file exists (0 is false, 1 is true)
 */
int FileExists(const char *fname) {
    FILE *file;
    if (file = fopen(fname, "r")) {
        fclose(file);
        return 1;
    }
    return 0;
}

/*
 *  Checks if folder exists (0 is false, 1 is true)
 */
int DirectoryExists(char * inputTest) {
    DIR * dir = opendir(inputTest);
    if (dir)
    {
        closedir(dir);
        return 1;
    }

    return 0;
}

/*
 *  Looks for bowtie2 indexes specified in dbPath
 */
void LookForBowtieDatabases(void) {
    DIR * dir = opendir(parms->database);
    struct dirent *ent;
    int indexes = 0;
    int i = 0;

    while ((ent = readdir(dir)) != NULL) {
        //looks for bowtie2 small indexes
        if (strstr(ent->d_name, ".1.bt2") != NULL && strstr(ent->d_name, "rev.1.bt2") == NULL && parms->largeGenome == 0 && strstr(ent->d_name, ".1.bt2l") == NULL)
            indexes++;

        //looks for bowtie2 large inddexes
        else if(parms->largeGenome == 1 && strstr(ent->d_name, ".1.bt2l") != NULL && strstr(ent->d_name, "rev.1.bt2l") == NULL)
            indexes++;

        //looks for bowtie2 small and large indexes
        if(parms->largeGenome == 2) {
            if (strstr(ent->d_name, ".1.bt2") != NULL && strstr(ent->d_name, "rev.1.bt2") == NULL)
                indexes++;
        }
    }

    closedir(dir);
    db->names = (char **)calloc(indexes, sizeof(char *));
    dir = opendir(parms->database);

    while ((ent = readdir(dir)) != NULL) {

        if(parms->largeGenome == 0) {
            if (strstr(ent->d_name, ".1.bt2") != NULL && strstr(ent->d_name, "rev.1.bt2") == NULL && strstr(ent->d_name, ".1.bt2l") == NULL) {
                db->names[i] = (char *)calloc(strlen(ent->d_name) - 5, sizeof(char));
                strncpy(db->names[i], ent->d_name, strlen(ent->d_name) - 6);
                i++;
            }
        }

        else if(parms->largeGenome == 1) {
            if( strstr(ent->d_name, ".1.bt2l") != NULL && strstr(ent->d_name, ".rev.1.bt2l") == NULL) {
                db->names[i] = (char *)calloc(strlen(ent->d_name), sizeof(char));
                strncpy(db->names[i], ent->d_name, strlen(ent->d_name) - 7);
                printf("name: %s\n", db->names[i]);
                i++;
            }
        }

        else {
            if (strstr(ent->d_name, ".1.bt2") != NULL && strstr(ent->d_name, "rev.1.bt2") == NULL) {
                db->names[i] = (char *)calloc(strlen(ent->d_name) - 4, sizeof(char));
                strncpy(db->names[i], ent->d_name, strlen(ent->d_name) - 6);
                i++;
            }
        }
    }

    closedir(dir);

    db->indexes = indexes;
}

/*
 *  Parses the anme of bowtie2 indexes
 */
void ParseBowtieDatabaseNames(void) {
    char * ptr;
    int i = 0;
    int indexes = 1;
    char * temp_names = (char *)calloc(strlen(parms->bwt2Indexes) + 1, sizeof(char));
    strcpy(temp_names, parms->bwt2Indexes);
    ptr = temp_names;

    while((ptr = strchr(ptr, ',')) != NULL) {
        indexes++;
        ptr++;
    }

    ptr = strtok(temp_names, ",");
    db->names = (char **)calloc(indexes, sizeof(char *));

    while(ptr != NULL) {
        db->names[i] = (char *)calloc(strlen(ptr) + 1, sizeof(char));
        strcpy(db->names[i], ptr);
        i++;
        ptr = strtok(NULL, ",");
    }

    db->indexes = indexes;

    if(temp_names) free(temp_names);
}

/*
 *  Checks if bowtie2 index files exist
 */
void CheckBowtieIndexes(void) {
    int i = 0;

    if(db->indexes == 0) {
        printf("No indexes found in %s\n", parms->database);
    }

    db->paths = (char **)calloc(db->indexes, sizeof(char *));

    for(i = 0; i < db->indexes; i++) {
        db->paths[i] = (char *)calloc(strlen(parms->database) + strlen(db->names[i]) + 2, sizeof(char));
        sprintf(db->paths[i], "%s/%s", parms->database, db->names[i]);

        if(FileExists(db->paths[i]) == 0) {
            printf("Database file %s does not exist\n", db->paths[i]);
            FreeInputs();
            exit(EXIT_FAILURE);
        }
    }
}

/*
 *  Prints run parameters (if verbose mode specified)
 */
void PrintParameters(void) {
    int i;
    printf("\n********************************************************************************\n");
    printf("Database path: %s\n", parms->database);
    printf("Database indexes:\n");

    for(i = 0; i < db->indexes; i++)
        printf("\t%s\n", db->names[i]);

    printf("\nInput reads: %s\n", parms->input);

    if(parms->pairedReads != NULL)
        printf("Input read pairs: %s\n", parms->pairedReads);

    if(parms->pairedReads) {
        printf("\n***************************************\n");
        printf("\tPAIRED-END information\n");
        printf("Min insert size: %d\n", parms->minInsert);
        printf("Max insert size: %d\n", parms->maxInsert);
        printf("***************************************\n\n");
    }

    if(parms->fastaInput == 1)
        printf("Input is fasta format\n");

    printf("Taxon nodes:\t\t\t%s\n", parms->names);

    printf("Output folder:\t\t\t%s\n", parms->output);
    printf("Host filter:\t\t\t%d\n", parms->hostFilter);

    if(parms->bowtieAllHits == 1)
        printf("Alignment hits:\t\t\tALL\n");

    else
        printf("Alignment hits:\t\t\t%d\n", parms->bowtieMaxHits);

    if(parms->bowtie2_params) {
        printf("Using extra bowtie2 command: %s\n", parms->bowtie2_params);
    }

    printf("Megan output:\t\t\t%d\n", parms->meganOutput);
    printf("Number of threads:\t\t%d\n", parms->threads);
    printf("Nearest neighbor score:\t\t%f\n", parms->neighborScore);
    if(parms->neighborOnly == 0)
        printf("Analysis:\t\t\tFull\n");

    else
        printf("Analysis:\t\t\tnearest neighbor only\n");

    if(parms->largeGenome == 1) {
        printf("Genome size:\t\t\tlarge\n");
    }

    else if(parms->largeGenome == 0) {
        printf("Genome size:\t\t\tsmall\n");
    }

    else {
        printf("Genome size:\t\t\tsmall and large\n");
    }

    printf("\n********************************************************************************\n\n");
}

/*
 *  Checks if input files are specified and exist
 */
void CheckInputs (void) {
    int errors = 0;
    FILE * in;
    char line[BUFSIZ];
    int lcount = 0;
    char * ptr;

    if(parms->database == NULL) {
        printf("No database path specified\n");
        errors++;
    }

    if(parms->input == NULL) {
        printf("No input reads (fastq) specified\n");
        errors++;
    }

    if(parms->output == NULL) {
        printf("No output folder specified\n");
        errors++;
    }

    if(parms->names == NULL) {
        printf("No taxon nodes (eg. nodes.dmp from NCBI) specifies\n");
        errors++;
    }

    //printf("%s\n", parms->bowtie2_params_file);

    if(parms->bowtie2_params_file) {
        if(FileExists(parms->bowtie2_params_file) == 0) {
            printf("No bowtie2 extra parameter file found: %s\n", parms->bowtie2_params_file);
            errors++;
        }

        in = fopen(parms->bowtie2_params_file, "r+");
        while (fgets(line, sizeof (line), in)) {
            if(lcount == 0) {
                parms->bowtie2_params = (char *)calloc(strlen(line) + 1, sizeof(char));
                strcpy(parms->bowtie2_params, line);
                parms->bowtie2_params[strcspn(parms->bowtie2_params, "\r\n")] = 0;
                if(strstr(parms->bowtie2_params, "-local") != NULL) {
                    parms->bowtie2_local = 1;//lorinc
                }
            }

            lcount++;
        }
        fclose(in);

        //printf("%s\n", parms->bowtie2_params);
    }

    if(errors > 0) {print_usage(stdout, 0);}


    if(DirectoryExists(parms->database) == 0) {
        printf("No such database directory: %s\n", parms->database);
        errors++;
    }

    if(parms->bwt2Indexes == NULL) {
        if(parms->verbose == 1)
            printf("No database names specified, looking for bowtie2 indexes at: %s\n", parms->database);
        LookForBowtieDatabases();
    }

    else
        ParseBowtieDatabaseNames();

    CheckBowtieIndexes();

    if(parms->verbose == 1)
        PrintParameters();

    if(FileExists(parms->input) == 0) {
        printf("No input reads found: %s\n", parms->input);
        errors++;
    }

    if(parms->pairedReads) {
        if(FileExists(parms->pairedReads) == 0) {
            printf("No input read pairs found: %s\n", parms->pairedReads);
            errors++;
        }
    }

    if(errors > 0) {
        print_usage(stdout, 0);
    }

    if(DirectoryExists(parms->output) == 0)
        mkdir(parms->output, 0700);

    parms->alignOut = (char *)calloc(strlen(parms->output) + strlen("Alignments") + 2, sizeof(char));
    sprintf(parms->alignOut, "%s/Alignments", parms->output);

    parms->results = (char *)calloc(strlen(parms->output) + strlen("/Results") + 1, sizeof(char));
    sprintf(parms->results, "%s/Results", parms->output);

    if(DirectoryExists(parms->alignOut) == 0)
        mkdir(parms->alignOut, 0700);

    if(DirectoryExists(parms->results) == 0)
        mkdir(parms->results, 0700);
}
//-d /home/lorinc//bioinformatica/newtaxon/db/DbCreator/bowtie -s anthrax.fastq  -o test -v  -i 0.fasta,1.fasta -n /home/lorinc//bioinformatica/newtaxon/db/DbCreator/nodes.dmp -b 20 -t 1
/*
 *  Reads inputs from command line
 */
void readInputs(int argc, char** argv) {
    const char* const short_options = "hvmfeualAd:n:s:o:b:i:r:I:X:p:w:t:k:c:y:";
    int next_option;

    parms = CreateInputs();
    db = (struct database *)calloc(1, sizeof(struct database));
    db->maxNode = 0;

    const struct option long_options[] = {
        { "help", 0, NULL, 'h'},
        { "verbose", 0, NULL, 'v'},
        { "dbPath", 1, NULL, 'd'},
        { "no-host-filter", 0, NULL, 'c'},
        { "seq", 1, NULL, 's'},
        { "output", 1, NULL, 'o'},
        { "bt2-maxhits", 1, NULL, 'b'},
        { "bt2-allhits", 1, NULL, 'a'},
        { "only-neighbor", 0, NULL, 'e'},
        { "bwt2-indexes", 1, NULL, 'i'},
        { "neighbor-score", 1, NULL, 'r'},
        { "fasta", 0, NULL, 'f'},
        { "megan", 0, NULL, 'm'},
        { "minInsert", 1, NULL, 'I'},
        { "maxInsert", 1, NULL, 'X'},
        { "paired", 1, NULL, 'p'},
        { "bowtie2", 1, NULL, 'w'},
        { "threads", 1, NULL, 't'},
        { "largeGenome", 0, NULL, 'l'},
        { "taxpath", 1, NULL, 'n'},
        { "skip", 1, NULL, 'k'},
        { "no-virus", 0, NULL, 'u'},
        { "both-genome-sizes", 0, NULL, 'g'},
        { "bwt2-params", 1, NULL, 'y'},
        { "alignstats", 0, NULL, 'A'},
        { NULL, 0, NULL, 0} /* Required at end of array.  */
    };

    do {
        next_option = getopt_long(argc, argv, short_options, long_options, NULL);

        switch (next_option) {
            case 'h':
                print_usage(stdout, 0);

            case 'v':
                parms->verbose = 1;
                break;

            case 'd':
                parms->database = (char *) malloc(sizeof (char) * (strlen(optarg) + 1));
                strcpy(parms->database, optarg);
                break;

            case 'n':
                parms->names = (char *) malloc(sizeof (char) * (strlen(optarg) + 1));
                strcpy(parms->names, optarg);
                break;

            case 'c':
                parms->hostFilter = 0;
                parms->hostindex = (char *) malloc(sizeof (char) * (strlen(optarg) + 1));
                strcpy(parms->hostindex, optarg);
                break;

            case 's':
                parms->input = (char *) malloc(sizeof (char) * (strlen(optarg) + 1));
                strcpy(parms->input, optarg);
                break;

            case 'o':
                parms->output = (char *) malloc(sizeof (char) * (strlen(optarg) + 1));
                strcpy(parms->output, optarg);
                break;

            case 'b':
                parms->bowtieMaxHits = atoi(optarg);
                break;

            case 'a':
                parms->bowtieAllHits = 1;
                break;

            case 'e':
                parms->neighborOnly = 1;
                break;

            case 'i':
                parms->bwt2Indexes = (char *) malloc(sizeof (char) * (strlen(optarg) + 1));
                strcpy(parms->bwt2Indexes, optarg);
                break;

            case 'r':
                parms->neighborScore = atof(optarg);
                break;

            case 'f':
                parms->fastaInput = 1;
                break;

            case 'm':
                parms->meganOutput = 1;
                break;

            case 'I':
                parms->minInsert = atoi(optarg);
                break;

            case 'X':
                parms->maxInsert = atoi(optarg);
                break;

            case 't':
                parms->threads = atoi(optarg);
                break;

            case 'w':
                parms->bowtieExec = (char *) malloc(sizeof (char) * (strlen(optarg) + 1));
                strcpy(parms->bowtieExec, optarg);
                break;

            case 'p':
                parms->pairedReads = (char *) malloc(sizeof (char) * (strlen(optarg) + 1));
                strcpy(parms->pairedReads, optarg);
                break;

            case 'l':
                parms->largeGenome = 1;
                break;

            case 'k':
                parms->removeTaxons = (char *) malloc(sizeof (char) * (strlen(optarg) + 1));
                strcpy(parms->removeTaxons, optarg);
                break;

            case 'u':
                parms->filterVirus = 1;
                break;

            case 'A':
                parms->alnstats = 1;
                break;

            case 'y':
                parms->bowtie2_params_file = (char *) malloc(sizeof (char) * (strlen(optarg) + 1));
                strcpy(parms->bowtie2_params_file, optarg);
                break;
        }
    } while (next_option != -1);

    CheckInputs();
}

 int64_t timespecDiff(struct timespec *timeA_p, struct timespec *timeB_p) {
    return ((timeA_p->tv_sec * 1000000000) + timeA_p->tv_nsec) -
            ((timeB_p->tv_sec * 1000000000) + timeB_p->tv_nsec);
}
