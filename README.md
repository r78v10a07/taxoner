# Taxoner

*Taxoner* publications (<font size="1">If you use *Taxoner* as part of a paper, please include the reference above</font>):
    * Pongor LS, Vera R, Ligeti B. (2014) <a href="http://www.plosone.org/article/info%3Adoi%2F10.1371%2Fjournal.pone.0103441">Fast and Sensitive Alignment of Microbial Whole Genome Sequencing Reads to Large Sequence Datasets on a Desktop PC: Application to Metagenomic Datasets and Pathogen Identification</a>. <i>PLoS ONE</i> 9(7): e103441.

*Taxoner*, a simple, parallelizable workflow that uses database partitioning in conjunction with [http://bowtie-bio.sourceforge.net/bowtie2/index.shtml Bowtie2] searches, which allows one to align millions of reads against the full [ftp://ftp.ncbi.nih.gov/blast/db/FASTA/nt.gz nt] database of [https://www.ncbi.nlm.nih.gov NCBI], typically in a few hours on a standard four threaded processor. 

The program can accept reads of all major platforms (Illumina, Ion Torrent, 454, !PacBio) and the output can be either further analyzed with sophisticated tools such as [http://ab.inf.uni-tuebingen.de/software/megan/ MEGAN], or directly evaluated by a built-in taxon assignment procedure that also allows the identification of genes and potential functions. 

The speed and accuracy of the method compares favorably to other programs and it can run on ordinary desktop or laptop computers without significant loss in speed.  

A demo server for testing porpoises is available in [http://pongor.itk.ppke.hu/taxoner/ here] 

http://taxoner.googlecode.com/svn/wiki/image/code3.jpg

http://taxoner.googlecode.com/svn/wiki/image/code7.jpg
