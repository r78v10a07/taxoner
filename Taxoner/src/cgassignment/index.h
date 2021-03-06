/* 
 * File:   index.h
 * Author: roberto
 *
 * Created on February 6, 2014, 10:21 AM
 */

#ifndef INDEX_H
#define	INDEX_H

#ifdef	__cplusplus
extern "C" {
#endif

    typedef struct {
        int gi;
        off_t offset;
    } giOffset_t;

    extern void createBTreeIndex(char *text, char *bin, char *index, char *output, int verbose);


#ifdef	__cplusplus
}
#endif

#endif	/* INDEX_H */

