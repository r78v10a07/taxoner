#ifndef NODES_H
#define	NODES_H

#ifdef	__cplusplus
extern "C" {
#endif
    extern void AllocateNodes(void);
    extern void GetMaxNodes(char * infile);
    extern void FillNodes(char * infile);
    extern void ImportNodes();
#ifdef	__cplusplus
}
#endif
#endif	/* NODES_H */
