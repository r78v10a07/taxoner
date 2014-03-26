#!/bin/bash

if [ -e "databases" ]
then
   cd databases
   wget ftp://ftp.ncbi.nih.gov/pub/taxonomy/taxdmp.zip
   unzip -d ./ taxdmp.zip nodes.dmp names.dmp
   grep "scientific name" names.dmp | awk -F\\t\|\\t '{printf("%s\t%s\n",$1,$3)}' > names
   rm taxdmp.zip names.dmp
   wget -r --reject "index.html*" -e robots=off -np -nH --cut-dirs=1 http://pongor.itk.ppke.hu/~roberto/geneassignment/   
   wget -r --reject "index.html*" -e robots=off -np -nH --cut-dirs=2 http://pongor.itk.ppke.hu/~roberto/taxoner/bowtie2/
else
   echo "Please, go to the Taxoner folder"
fi

exit 0