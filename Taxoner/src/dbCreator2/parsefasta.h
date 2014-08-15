#ifndef PARSEFASTA_H
#define	PARSEFASTA_H

#ifdef	__cplusplus
extern "C" {
#endif
    extern void OpenFastaFile(void);
    extern void CountData(int a);
    extern int FindTaxonWithGi(int gids);
    extern int compareStrings(char * a, char * b);
    extern void ParseFastaTitle(char * input);
    extern void ReadFasta(char * inFile);
#ifdef	__cplusplus
}
#endif

#endif	/* PARSEFASTA_H */


