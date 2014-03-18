/* 
 * File:   define.h
 * Author: roberto
 *
 * Created on February 4, 2014, 10:04 PM
 */

#ifndef DEFINE_H
#define	DEFINE_H

#ifdef	__cplusplus
extern "C" {
#endif

    /*
        typedef struct {
            int position[2];
        } location_t;
     */
    typedef struct {
        char *proteinId;
        char *locusName;
        int pFrom;
        int pTo;
        int hits;
        int cog_number;
        char **cog;
        int protClust_number;
        char **protClust;
    } cds_t;

    typedef struct {
        int gi;
        char *locusName;
        int taxId;
        int cds_number;
        cds_t *cds;
    } genBank_t;

    extern genBank_t *allocGenBank();
    extern void freeGenBank(genBank_t *g);
    extern genBank_t *initGenBank(int gi, char *glocusName, int taxId, char *proteinId, char *locusName, int pFrom, int pTo, char *cog, char *protClust);
    extern void addCDS(cds_t **cds, int *index, char *proteinId, char *locusName, int pFrom, int pTo, char *cog, char *protClust);
    extern void printGenBank(genBank_t *g);
    extern void printGenBankBinary(FILE *fb, genBank_t *g);
    extern genBank_t *readGenBankBinary(FILE *fb);
    extern genBank_t *readGenBankBinaryOffSet(FILE *fb, long int offset);
    extern int assignGenes(node **root, genBank_t **gs, int *count, FILE *fb, off_t offset, int pFrom, int pTo);
    extern int cmpGenBank(const void *a, const void * b);
    extern int cmpCDSSum(const void *a, const void * b);

#ifdef	__cplusplus
}
#endif

#endif	/* DEFINE_H */

