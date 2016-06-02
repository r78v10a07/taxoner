#ifndef IMPORTGI_H
#define	IMPORTGI_H

#ifdef	__cplusplus
extern "C" {
#endif
    extern void GetMaxGiNuclDmp(char * infile);
    extern void ParseGiNuclLine(char * input);
    extern void ImportGiNuclDmp(char * infile);
    extern void FreeGiId(void);
#ifdef	__cplusplus
}
#endif

#endif	/* IMPORTGI_H */

