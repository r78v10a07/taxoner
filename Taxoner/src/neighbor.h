#ifndef NEIGHBOR_H
#define	NEIGHBOR_H

#ifdef	__cplusplus
extern "C" {
#endif
    extern char * ReturnName(char * input, char * dest);
    extern int CheckIfAligned(char * input);
    extern float GetAlignmentScore(char * input, char * cigar, int seqlentemp);
    extern int ReturnAlignmentScore(char * input);
    extern void AllocateHead(void);
    extern void AllocateData(void);
    extern int CompareStrings(char * s1, char * s2);
    extern void AddRead(char * source, char * name);
    extern void FreeStruct(void);
    extern int Blacklist(int a);
    extern void FreeReadElement(struct Reads * ptr);
    extern int ReturnNeighbor(int a, int b);
    extern void ParseRead(struct Reads * psr);
    extern void ParseAndFilterReads(struct Data * ptr);
    extern void ParseReadFinal(struct Reads * psr);
    extern void ParseAndFilterFinalReads(struct Data * ptr);
    extern void * NearestNeighbor(void * voidA);
    extern void PrintResults(FILE * output);
    extern void * ParseAndFindNeighbor(void * voidA);
    extern void InitializeNearestMulti(struct Data * psr);
    extern void * ParseAndFindFinalNeighbor(void * voidA);
    extern void InitializeFinalNearestMulti(struct Data * psr);
    extern void ReadAlnFile(char * infile, char * outfile);
    extern void ReadSamFile(char * infile, FILE * out);
    extern void StartNeighbor(void);
    extern void PrintMeganFormat(void);
    extern void GetNeighbors(void);
#ifdef	__cplusplus
}
#endif
#endif	/* NEIGHBOR_H */

