/* 
 * File:   main.c
 * Author: roberto
 *
 * Created on February 4, 2014, 3:09 PM
 */
#define _GNU_SOURCE
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <getopt.h>
#include <stdbool.h>
#include <time.h>
#include "btree.h"
#include "writer.h"
#include "index.h"
#include "genbank.h"

char *program_name;

void print_usage(FILE *stream, int exit_code) {
    fprintf(stream, "\n********************************************************************************\n");
    fprintf(stream, "\nUsage: %s \n", program_name);
    fprintf(stream, "\n\n%s options:\n\n", program_name);
    fprintf(stream, "-h,   --help                        Display this usage information.\n");
    fprintf(stream, "-w,                                 Set the process to create the binary file and index\n");
    fprintf(stream, "-t,   --text                        The input text file\n");
    fprintf(stream, "-b,   --bin                         The binary file\n");
    fprintf(stream, "-i,   --index                       The index file\n");
    fprintf(stream, "-o,   --output                      The output file\n");
    fprintf(stream, "********************************************************************************\n");
    fprintf(stream, "\n            Roberto Vera Alvarez (e-mail: r78v10a07@gmail.com)\n\n");
    fprintf(stream, "********************************************************************************\n");
    exit(0);
}

/*
 *  Main function
 */
int main(int argc, char** argv) {

    int next_option, verbose, write;
    const char* const short_options = "hvwi:b:t:o:";
    char *text, *bin, *index, *output;

    program_name = argv[0];

    const struct option long_options[] = {
        { "help", 0, NULL, 'h'},
        { "verbose", 0, NULL, 'v'},
        { "write", 0, NULL, 'w'},
        { "text", 1, NULL, 't'},
        { "bin", 1, NULL, 'b'},
        { "index", 1, NULL, 'i'},
        { "output", 1, NULL, 'o'},
        { NULL, 0, NULL, 0} /* Required at end of array.  */
    };

    write = verbose = 0;
    text = bin = index = output = NULL;
    do {
        next_option = getopt_long(argc, argv, short_options, long_options, NULL);

        switch (next_option) {
            case 'h':
                print_usage(stdout, 0);

            case 'v':
                verbose = 1;
                break;

            case 'w':
                write = 1;
                break;

            case 'o':
                output = (char *) malloc(sizeof (char) * (strlen(optarg) + 1));
                strcpy(output, optarg);
                break;
                break;

            case 't':
                text = (char *) malloc(sizeof (char) * (strlen(optarg) + 1));
                strcpy(text, optarg);
                break;

            case 'b':
                bin = (char *) malloc(sizeof (char) * (strlen(optarg) + 1));
                strcpy(bin, optarg);
                break;

            case 'i':
                index = (char *) malloc(sizeof (char) * (strlen(optarg) + 1));
                strcpy(index, optarg);
                break;
        }
    } while (next_option != -1);

    if (text == NULL) {
        fprintf(stdout, "Set the input text file\n");
        print_usage(stdout, 1);
    }

    if (bin == NULL) {
        fprintf(stdout, "Set the output bin file\n");
        print_usage(stdout, 1);
    }

    if (index == NULL) {
        fprintf(stdout, "Set the index file\n");
        print_usage(stdout, 1);
    }

    if (write == 1) {
        writer(text, bin, index, verbose);
    } else {
        if (output == NULL) {
            fprintf(stdout, "Set the output file\n");
            print_usage(stdout, 1);
        }
        createBTreeIndex(text, bin, index, output, verbose);
    }
    if (text) free(text);
    if (bin) free(bin);
    if (index) free(index);
    if (output) free(output);
    return (EXIT_SUCCESS);
}

