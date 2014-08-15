#ifndef INCLUDEEXCLUDE_H
#define	INCLUDEEXCLUDE_H

#ifdef	__cplusplus
extern "C" {
#endif
    extern int CheckInclude(int source, int status);
    extern int CheckExclude(int source, int status);
    extern void ImportInclude(char * infile);
    extern void ImportExclude(char * infile);
#ifdef	__cplusplus
}
#endif

#endif	/* INCLUDEEXCLUDE_H */


