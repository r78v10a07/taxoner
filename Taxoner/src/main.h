#ifndef MAIN_H
#define	MAIN_H

#ifdef	__cplusplus
extern "C" {
#endif

/*
 * Input parameter structure
 */
struct params {
    int verbose;
    char * database;
    char * names;
    int hostFilter;
    char * input;
    char * output;
    int bowtieMaxHits;
    int bowtieAllHits;
    int neighborOnly;
    char * bwt2Indexes;
    float neighborScore;
    int fastaInput;
    int meganOutput;
    int minInsert;
    int maxInsert;
    char * pairedReads;
    char * bowtieExec;
    int threads;
    int largeGenome;
    char * alignOut;
    int filterVirus;
    char * removeTaxons;
    char * results;
    char * hostindex;
    char * bowtie2_params_file;
    char * bowtie2_params;
    int bowtie2_local;//lorinc
    int alnstats;
};

struct database {
    char ** paths;
    char ** names;
    int indexes;

    int maxNode;

};

struct TaxonNodes {
    int curr;
    int next;
    int ok;
};

#ifdef	__cplusplus
}
#endif

#endif	/* MAIN_H */
