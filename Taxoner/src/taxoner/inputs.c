char * host = NULL;
char * dbPath = NULL;
char ** dbElements = NULL;
char **fulldb = NULL;
int dbe = 0;
char * taxPath = NULL;
int threads = 1;
int filterHost = 1;
char * reads = NULL;
char nohost[] = "nonhost.fq";
char * outDir = NULL;
char * outaln = NULL;
char * outres = NULL;
int maxTax = 0;
time_t timer;
char times[30];
struct tm* tm_info = NULL;
int bowtieMaxHits = 10;
int allHits = 0;
float NSc = 0.99;

int meganOut = 0;
int inputfasta = 0;
int noUnaligned = 0;

int minfrag = 0;
int maxfrag = 500;
char * pairedFile = NULL;

int pairedend = 0;

char * specifiedDb = NULL;

int NoAlignment = 0;

int FilterVirus = 1;

char * bowtieFolder = NULL;

char * InputName = NULL;
char * OutputName = NULL;

void InputFree(void) {
    int i = 0;

    if (host)
        free(host);

    if (dbPath)
        free(dbPath);

    if (taxPath)
        free(taxPath);

    if (reads)
        free(reads);

    if (outDir)
        free(outDir);

    if (outaln)
        free(outaln);

    if (outres)
        free(outres);

    if (InputName)
        free(InputName);

    if (OutputName)
        free(OutputName);

    if (dbElements) {
        for (i = 0; i < dbe; i++) {
            if (dbElements[i]) free(dbElements[i]);
        }
        if (dbElements) free(dbElements);
    }

    if (fulldb) {
        for (i = 0; i < dbe; i++) {
            if (fulldb[i]) free(fulldb[i]);
        }
        if (fulldb) free(fulldb);
    }

    if (bowtieFolder) {
        free(bowtieFolder);
    }
}

int isInit(char * source) {
    if (source != NULL)
        return 0;

    return 1;
}

void InitializedData(void) {
    int bad = 0;

    if (isInit(dbPath) == 1) {
        bad++;
        printf("Database path not specified\n");
    }

    if (isInit(taxPath) == 1) {
        bad++;
        printf("No taxonomy path specified\n");
    }

    if (isInit(reads) == 1) {
        bad++;
        printf("No input fastq/fasta file specified\n");
    }


    if (isInit(outDir) == 1) {
        bad++;
        printf("No output directory specified\n");
    }

    if (bad > 0) {
        printf("\nFound %d initialization problem(s)\nExiting\n", bad);
        exit(EXIT_FAILURE);
    }
}

int FileExists(const char *fname) {
    FILE *file;
    if (file = fopen(fname, "r")) {
        fclose(file);
        return 1;
    }
    return 0;
}

void CreateDbElements(void) {
    int i, j, k;
    char bowtie[] = ".1.bt2";

    if (dbElements == NULL)
        return;

    fulldb = (char **) calloc(dbe, sizeof (char *));

    for (i = 0; i < dbe; i++) {
        if (dbPath[strlen(dbPath) - 1] == '/') {
            fulldb[i] = (char *) calloc(strlen(dbPath) + strlen(bowtie) + strlen(dbElements[i]) + 2, sizeof (char));
            strncpy(fulldb[i], dbPath, strlen(dbPath));
            k = 0;
            for (j = strlen(dbPath); j < strlen(dbElements[i]) + strlen(dbPath) + 1; j++) {
                fulldb[i][j] = dbElements[i][k];
                k++;
            }

            k = 0;

            for (j = strlen(dbElements[i]) + strlen(dbPath); j < strlen(dbElements[i]) + strlen(dbPath) + strlen(bowtie); j++) {
                fulldb[i][j] = bowtie[k];
                k++;
            }
        } else {
            fulldb[i] = (char *) calloc(strlen(dbPath) + strlen(bowtie) + strlen(dbElements[i]) + 2, sizeof (char));
            strncpy(fulldb[i], dbPath, strlen(dbPath));
            fulldb[i][strlen(dbPath)] = '/';
            k = 0;
            for (j = strlen(dbPath) + 1; j < strlen(dbElements[i]) + strlen(dbPath) + 2; j++) {
                fulldb[i][j] = dbElements[i][k];
                k++;
            }

            k = 0;

            for (j = strlen(dbElements[i]) + strlen(dbPath) + 1; j < 1 + strlen(dbElements[i]) + strlen(dbPath) + strlen(bowtie); j++) {
                fulldb[i][j] = bowtie[k];
                k++;
            }
        }
        //printf("CreateDbElements: [%d] [%s]\n", i, fulldb[i]);
    }
}

void ReadFolder(void) {
    DIR *dir = opendir(dbPath);
    struct dirent *ent;
    int ok = 0;
    int i = 0;
    int length = strlen(dbPath);
    int temp;

    while ((ent = readdir(dir)) != NULL) {
        if (strstr(ent->d_name, ".1.bt2") != NULL && strstr(ent->d_name, "rev.1.bt2") == NULL)
            ok++;
    }

    closedir(dir);

    dir = opendir(dbPath);
    fulldb = (char **) calloc(ok, sizeof (char *));
    dbElements = (char **) calloc(ok, sizeof (char *));
    dbe = ok;

    while ((ent = readdir(dir)) != NULL) {
        if (strstr(ent->d_name, ".1.bt2") != NULL && strstr(ent->d_name, "rev.1.bt2") == NULL) {
            fulldb[i] = (char *) calloc(strlen(ent->d_name) + strlen(dbPath) + 2, sizeof (char));
            strncpy(fulldb[i], dbPath, strlen(dbPath));
            fulldb[i][strlen(dbPath)] = '\0'; // by Roberto

            dbElements[i] = (char *) calloc(strlen(ent->d_name) + 2, sizeof (char));
            strncpy(dbElements[i], ent->d_name, strlen(ent->d_name));

            dbElements[i][strlen(ent->d_name)] = '\0'; // by Roberto
            temp = length;

            if (fulldb[i][strlen(dbPath) - 1] != '/') {
                fulldb[i][strlen(dbPath)] = '/';
                temp++;
            }

            fulldb[i] += temp;
            strncpy(fulldb[i], ent->d_name, strlen(ent->d_name));
            fulldb[i][strlen(ent->d_name)] = '\0'; // by Roberto
            fulldb[i] -= temp;

            //printf("ReadFolder: [%d] [%s]\n", i, fulldb[i]);
            i++;
        }
    }

    closedir(dir);

}

void CreateSpecifiedDb(void) {

}

int CheckDbDirectory(DIR * inputTest)
{
    if (inputTest) {
        closedir(inputTest);
        return 1;

    } else {
        printf("No such database directory: %s\n", dbPath);
        InputFree();
        exit(EXIT_FAILURE);
    }

    return 0;
}


void checkParams(void) {
    DIR * dir = opendir(dbPath);
    int bad = 0;
    int i = 0;
    int temp = 0;
    char tempBuffer[BUFSIZ];
    InitializedData();
    CreateDbElements();

    if(CheckDbDirectory(dir) == 1) {
        if (dbElements == NULL)
            ReadFolder();
    }

    dir = opendir(outDir);

    if (dir) {
        closedir(dir);
        //printf("Deirectory: %s exists\n", outDir);
    } else {
        printf("Creating output directory: %s\n", outDir);
        mkdir(outDir, 0700);

        dir = opendir(outDir);

        if (dir) {

        } else {
            printf("Could not create directory: %s\nExiting (check permissions of program)\n", outDir);
            exit(EXIT_FAILURE);
        }

        closedir(dir);
    }

    outaln = (char *) calloc(strlen(outDir) + strlen("alignments") + 2, sizeof (char));
    outres = (char *) calloc(strlen(outDir) + strlen("Results") + 2, sizeof (char));

    strncpy(outaln, outDir, strlen(outDir));
    strncpy(outres, outDir, strlen(outDir));

    if (outDir[strlen(outDir) - 1] != '/') {
        outaln[strlen(outDir)] = '/';
        outres[strlen(outDir)] = '/';
        temp++;
    }

    outaln += strlen(outDir) + temp;
    outres += strlen(outDir) + temp;

    strncpy(outaln, "alignments", strlen("alignments"));
    strncpy(outres, "Results", strlen("Results"));

    outaln -= strlen(outDir) + temp;
    outres -= strlen(outDir) + temp;

    mkdir(outaln, 0700);
    mkdir(outres, 0700);

    dir = opendir(outaln);

    if (dir) {
        closedir(dir);
    } else {
        printf("Could not create directory: %s\nChweck permissions\n", outaln);
        exit(EXIT_FAILURE);
    }

    dir = opendir(outres);

    if (dir) {
        closedir(dir);
    } else {
        printf("Could not create directory: %s\nChweck permissions\n", outres);
        exit(EXIT_FAILURE);
    }

    if (FileExists(reads) == 0) {
        bad++;
        printf("No such read file: %s\n", reads);
    }


    if (FileExists(taxPath) == 0) {
        bad++;
        printf("No such taxonomy file exists: %s\n", taxPath);
    }

    for (i = 0; i < dbe; i++) {
        if (FileExists(fulldb[i]) == 0 && NoAlignment == 0) {
            bad++;
            printf("No bowtie index found: %s\n", fulldb[i]);
        }
    }

    if (pairedend == 1) {
        if (FileExists(pairedFile) == 0) {
            bad++;
            printf("No paired-end file named %s founf\n", pairedFile);
        }
    }

    if (bowtieFolder != NULL) {
        strcpy(tempBuffer, bowtieFolder);

        if (bowtieFolder[strlen(bowtieFolder) - 1] != '/')
            strcat(tempBuffer, "/");

        strcat(tempBuffer, "./bowtie2");

        if (CheckFile(tempBuffer) == 0) {
            printf("No executable found named: %s\n", tempBuffer);
            printf("Exiting\n");
            exit(EXIT_FAILURE);
        }
    }else{
        strcat(tempBuffer, "bowtie2");
    }

    if (bad > 0) {
        printf("Found %d problem(s) with input files\nExiting\n", bad);
        exit(EXIT_FAILURE);
    }
}

void PrintInputParams(void) {
    int i;
    checkParams();

    printf("======================================\n");
    printf("Taxoner version: 2.0\n\n");
    printf("Input parameters\n");
    printf("--------------------------------------\n");
    fflush(NULL);

    if (host != NULL) {
        printf("Host genome: %s\n", host);
        printf("Filter host reads: ");

        if (filterHost == 0)
            printf("no\n");

        else
            printf("yes\n");
    }

    printf("Database index folder: %s\n", dbPath);

    if (dbElements != NULL) {
        for (i = 0; i < dbe; i++)
            printf("\tDatabase index (%d): %s\n", i, dbElements[i]);
    }

    printf("Taxonomy file: %s\n", taxPath);
    
    if (noUnaligned == 0)
        printf("Keeping unaligned sequences from final output\n");
    else
        printf("Discarding unaligned sequences from final output\n");
    
    if (inputfasta == 0)
        printf("Read format: fastq format\n");

    else
        printf("Read format: fasta format\n");

    printf("Input reads: %s\n", reads);

    if (pairedend == 1) {
        printf("Paired-end reads: %s\n", pairedFile);
        printf("Min insert size: %d\nmaximum insert size: %d\n", minfrag, maxfrag);
    }

    printf("Threads for analysis: %d\n", threads);
    printf("--------------------------------------\n");
    fflush(NULL);
}

int IsEqual(char * a, char * b) {
    if (strlen(a) == strlen(b)) {
        if (strncmp(a, b, strlen(a)) == 0)
            return 0;
    }

    return 1;
}

char * CopyString(char * source, int length) {
    char * destination = (char *) calloc(length + 1, sizeof (char));
    strncpy(destination, source, length);

    return destination;
}

void ParseDbNames(char * source) {
    int elements = 1;
    int i;
    char * ptr;

    for (i = 0; i < strlen(source); i++) {
        if (source[i] == ',')
            elements++;
    }

    dbe = elements;
    dbElements = (char **) calloc(elements + 1, sizeof (char *));

    ptr = strtok(source, ",");
    i = 0;

    while (ptr != NULL) {
        dbElements[i] = (char *) calloc(strlen(ptr) + 1, sizeof (char));
        strncpy(dbElements[i], ptr, strlen(ptr));

        i++;
        ptr = strtok(NULL, ",");
    }

    return;
}

void HelpMessage(void) {
    printf("Taxoner version 2.0\n");
    printf("\tAlign NGS reads to large databases like NCBI nr/nt/complete genomes...\n");
    printf("\n\tPlease index database with other program first\n");
    printf("\nUsage: ./taxoner -dbPath <folder_path> -taxpath <taxonomy_file> -seq <NGS_reads or #1 mate reads> -o <output directory>\n");
    printf("\nOther Commands:\n\n");
    printf(" Input read options\n");
    printf("\t-seq\t\tInput file with NGS reads or #1 mate pair for paired-end\n");
    printf("\t-paired\t\tFile with #2 mate reads [file name]\n");
    printf("\t-fasta\t\tInput reads are in fasta format (not fastq)\n\n");

    printf(" Input database options\n");
    printf("\t-dbPath\t\tFolder path where bowtie2 indexes can be found [folder path]\n");
    printf("\t-taxpath\t\tNodes.dmp file from NCBI database\n");
    printf("\t-dbNames\tNames of databases (indexes) to include instead of all [eg. bact.0.fasta,bact.1.fasta]\n\n");

    printf(" Host filtering options\n");
    printf("\t-host\t\tindex of host genome [index path and name]\n");
    printf("\t-no-host-filter\tDo not filter host genome\n\n");

    printf(" Alignment options\n");
    printf("\t-bt2-maxhits\tNumber of alignments reported by bowtie2 [def: 10]\n");
    printf("\t-bt2-allhits\tTells bowtie2 to report all alignments (very slowww)\n");
    printf("\t-no-unal\tTell bowtie2 to discard all unaligned reads\n\n");

    printf(" Paired-end alignment options\n");
    printf("\t-I\t\tMinimum fragment length for mates [def: 0]\n");
    printf("\t-D\t\tMaximum fragment length for mates [def: 500]\n\n");

    printf(" Nearest neighbor options\n");
    printf("\t-only-neighbor\tTells Taxoner to look for neighbors in the specified output folder. Only if there was a prior alignment\n");
    printf("\t-neighbor-score\tAlignment score for neighbor lookup between alignments [def: 0.99]\n\n");

    printf(" Output options\n");
    printf("\t-o\t\tOutput folder\n");
    printf("\t-megan\t\tCreates a Megan compatible output\n\n");

    printf(" Bowtie2 path option (if no bowtie2 installed on system)\n");
    printf("\t-bowtie2\tIf bowtie2 is not installed in the OS, bowtie2 can be specified in command line [path with executable]\n\n");

    printf(" Performance options\n");
    printf("\t-p\t\tnumber of threads [def: 1]\n\n");
    //printf("\t-dbNames\tNames of databases (indexes) to include instead of all\n");
    //printf("\t-host\t\tindex of host genome\n");
    //printf("\t-bt2-maxhits\tNumber of alignments reported by bowtie2 (def: 10)\n");
    //printf("\t-bt2-allhits\tTells bowtie2 to report all alignments (very slowww)\n");
    //printf("\t-only-neighbor\tTells Taxoner to look for neighbors in the specified output folder. Only if there was a prior alignment\n");
    //printf("\t-neighbor-score\tAlignment score for neighbor lookup between alignments (def: 0.99)\n");
    //printf("\t-fasta\t\tInput reads are in fasta format (not fastq)\n");
    //printf("\t-megan\t\tCreates a Megan compatible output\n");
    //printf("\t-paired\t\tFile with #2 mate reads\n");
    //printf("\t-I\t\tMinimum fragment length for mates (def: 0)\n");
    //printf("\t-D\t\tMaximum fragment length for mates (def: 500)\n");
    //printf("\t-bowtie2\tIf bowtie2 is not installed in the OS, bowtie2 can be specified in command line\n");
    //printf("\t-no-host-filter\tDo not filter host genome\n");
    //printf("\t-o\t\tOutput folder\n\n");

    exit(EXIT_SUCCESS);
}



void CheckCommands(char * source[], int num) {
    int i;
    if (num == 1)
        HelpMessage();

    for (i = 0; i < num; i++) {
        if (IsEqual(source[i], "-host") == 0 && (i + 1) < num)
            host = CopyString(source[i + 1], strlen(source[i + 1]));

        else if (IsEqual(source[i], "-dbPath") == 0 && (i + 1) < num)
            dbPath = CopyString(source[i + 1], strlen(source[i + 1]));

        else if (IsEqual(source[i], "-p") == 0 && (i + 1) < num)
            threads = atoi(source[i + 1]);

        else if (IsEqual(source[i], "-dbNames") == 0 && (i + 1) < num)
            ParseDbNames(source[i + 1]);

        else if (IsEqual(source[i], "-taxpath") == 0 && (i + 1) < num)
            taxPath = CopyString(source[i + 1], strlen(source[i + 1]));

        else if (IsEqual(source[i], "-no-host-filter") == 0)
            filterHost = 0;

        else if (IsEqual(source[i], "-seq") == 0 && (i + 1) < num)
            reads = CopyString(source[i + 1], strlen(source[i + 1]));

        else if (IsEqual(source[i], "-h") == 0)
            HelpMessage();

        else if (IsEqual(source[i], "-help") == 0)
            HelpMessage();

        else if (IsEqual(source[i], "-o") == 0 && (i + 1) < num)
            outDir = CopyString(source[i + 1], strlen(source[i + 1]));

        else if (IsEqual(source[i], "-bt2-maxhits") == 0 && (i + 1) < num)
            bowtieMaxHits = atoi(source[i + 1]);

        else if (IsEqual(source[i], "-bt2-allhits") == 0)
            allHits = 1;

        else if (IsEqual(source[i], "-only-neighbor") == 0)
            NoAlignment = 1;

        else if (IsEqual(source[i], "-bowtie2-indexes") == 0)
            specifiedDb = CopyString(source[i + 1], strlen(source[i + 1])); //printf("%s\n", source[i + 1]);

        else if (IsEqual(source[i], "-neighbor-score") == 0)
            NSc = atof(source[i + 1]);

        else if (IsEqual(source[i], "-fasta") == 0)
            inputfasta = 1;
        
        else if (IsEqual(source[i], "-no-unal") == 0)
            noUnaligned = 1;

        else if (IsEqual(source[i], "-megan") == 0)
            meganOut = 1;

        else if (IsEqual(source[i], "-I") == 0)
            minfrag = atoi(source[i + 1]);

        else if (IsEqual(source[i], "-X") == 0)
            maxfrag = atoi(source[i + 1]);

        else if (IsEqual(source[i], "-paired") == 0) {
            pairedFile = CopyString(source[i + 1], strlen(source[i + 1]));
            pairedend = 1;
        }

        else if (IsEqual(source[i], "-bowtie2") == 0)
            bowtieFolder = CopyString(source[i + 1], strlen(source[i + 1]));

        else {
            if(strncmp(source[i], "-", 1) == 0) {
                printf("Unknown command: \"%s\"\n", source[i]);
                InputFree();
                HelpMessage();
            }
        }

    }

    PrintInputParams();

    //exit(EXIT_SUCCESS);
}
