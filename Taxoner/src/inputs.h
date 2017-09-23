#ifndef INPUTS_H
#define	INPUTS_H

#ifdef	__cplusplus
extern "C" {
#endif
    extern void print_usage(FILE *stream, int exit_code);
    extern void FreeInputs(void);
    extern struct params * CreateInputs(void);
    extern int FileExists(const char *fname);
    extern int DirectoryExists(char * inputTest);
    extern void LookForBowtieDatabases(void);
    extern void ParseBowtieDatabaseNames(void);
    extern void CheckBowtieIndexes(void);
    extern void PrintParameters(void);
    extern void CheckInputs (void);
    extern struct params * readInputs(int argc, char** argv);
    extern int64_t timespecDiff(struct timespec *timeA_p, struct timespec *timeB_p);

#ifdef	__cplusplus
}
#endif

#endif	/* INPUTS_H */
