/* 
 * File:   index.h
 * Author: roberto
 *
 * Created on April 5, 2014, 11:50 AM
 */

#ifndef INDEX_H
#define	INDEX_H

#ifdef	__cplusplus
extern "C" {
#endif

    extern node *createGiIndex(char * giFileName);
    extern node *createTaxIndex(char * filename);
    extern void freeBtreeData(void *d);

#ifdef	__cplusplus
}
#endif

#endif	/* INDEX_H */

