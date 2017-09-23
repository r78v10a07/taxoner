#ifndef ALIGNMENTS_H
#define	ALIGNMENTS_H

#ifdef	__cplusplus
extern "C" {
#endif
    extern char * ReturnBowtieCommand(char * InputIndex, char * indexName)
    extern char * ReturnBowtieCommandHost(char * InputIndex, char * indexName)
    extern void StartAlignment(void);
    extern void AlignReads(void);
#ifdef	__cplusplus
}
#endif
#endif	/* ALIGNMENTS_H */

