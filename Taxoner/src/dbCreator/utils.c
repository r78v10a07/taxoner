
#define _GNU_SOURCE
#define _FILE_OFFSET_BITS 64
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <stdbool.h>

#include "files.h"

extern char *program_name;

void print_usage(FILE *stream, int exit_code) {
    fprintf(stream, "\n********************************************************************************\n");
    fprintf(stream, "\nUsage: %s \n", program_name);
    fprintf(stream, "\n\n%s options:\n\n", program_name);
    fprintf(stream, "\t-h,\t--help\t\tDisplay this usage information.\n");
    fprintf(stream, "\t-n,\t--nt\t\tNT database fasta file (from NCBI)\n");
    fprintf(stream, "\t-g,\t--gi\t\tgi_taxid_nucl.dmp file (from NCBI)\n");
    fprintf(stream, "\t-o,\t--nodes\t\tnodes.dmp file (from NCBI)\n\n");
    fprintf(stream, " Optional:\n");
    fprintf(stream, "\t-s,\t--skip\t\tFile containing a list of taxons to skip\n");
    fprintf(stream, "\t-i,\t--include\tFile containing a list of taxons to include\n\n");
    fprintf(stream, "********************************************************************************\n");
    fprintf(stream, "\n            Roberto Vera Alvarez (e-mail: r78v10a07@gmail.com)\n\n");
    fprintf(stream, "\n            Lorinc Pongor (e-mail: pongorlorinc@gmail.com)\n\n");
    fprintf(stream, "********************************************************************************\n");
    exit(0);
}

void InputCheck(char *ntfile, char *gifile, char *nodesfile, char *skipfile, char *includefile) {
    if (ntfile) {
        CheckFiles(ntfile);
    } else {
        printf("No input NT database specified\n");
        print_usage(stderr, -1);
    }

    if (gifile) {
        CheckFiles(gifile);
    } else {
        printf("No input gi_taxid_nucl.dmp database specified\n");
        print_usage(stderr, -1);
    }

    if (nodesfile) {
        CheckFiles(nodesfile);
    } else {
        printf("No input nodes.dmp database specified\n");
        print_usage(stderr, -1);
    }

    if (skipfile)
        CheckFiles(skipfile);

    if (includefile)
        CheckFiles(includefile);
}

int64_t timespecDiff(struct timespec *timeA_p, struct timespec *timeB_p) {
    return ((timeA_p->tv_sec * 1000000000) + timeA_p->tv_nsec) -
            ((timeB_p->tv_sec * 1000000000) + timeB_p->tv_nsec);
}

char * CopyString(char * source) {
    char * destination = (char *) calloc(strlen(source) + 1, sizeof (char));
    strncpy(destination, source, strlen(source));
    return destination;
}

int CompareStrings(char * a, char * b) {
    if (strlen(a) == strlen(b)) {
        if (strncmp(a, b, strlen(a)) == 0)
            return 0;
    }
    return 1;
}
