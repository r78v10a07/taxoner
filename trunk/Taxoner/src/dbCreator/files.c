#define _GNU_SOURCE
#define _FILE_OFFSET_BITS 64
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>

#include "utils.h"

void CheckFiles(char * inptfile) {
    if (access(inptfile, F_OK) != -1) {
        return;
    } else {
        printf("No such file \"%s\", exiting\n", inptfile);
        print_usage(stderr, -1);
    }
}
