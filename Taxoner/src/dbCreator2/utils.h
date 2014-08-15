#ifndef UTILS_H
#define	UTILS_H

#ifdef	__cplusplus
extern "C" {
#endif

    extern void print_usage(FILE *stream, int exit_code);
    extern void print_parameters(void);
    extern void InputCheck(char *ntfile, char *gifile, char *nodesfile, char *skipfile, char *includefile);
    extern int64_t timespecDiff(struct timespec *timeA_p, struct timespec *timeB_p);
    extern char * CopyString(char * source);
    extern int CompareStrings(char * a, char * b);

#ifdef	__cplusplus
}
#endif

#endif	/* UTILS_H */

