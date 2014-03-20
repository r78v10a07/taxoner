#!/bin/bash

if [ -e "databases" ]
then
   cd databases
   wget ftp://ftp.ncbi.nih.gov/pub/taxonomy/taxdmp.zip
   unzip -d ./ taxdmp.zip
   cd geneassignment
   wget -np -nH -r --cut-dirs=1 http://pongor.itk.ppke.hu/~roberto/geneassigment/
else
   echo "Please, go to the Taxoner folder"
fi

exit 0