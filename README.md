# NFTOOLS
java package to localise recombination events in genetic data


*************************
**     NucFamTools     **
*************************




INTRODUCTION
------------

NucFamTools.jar is a java package to localise recombination events in pedigree genetic data. The algorithm is a modified version of the one described in Coop et al. Science (2008). Details are presented in Hussin et al. PLoS Genetics (2011).



INSTALLATION
------------

Options:

1) Download NFTOOLS.tar.gz

> mkdir NFTOOLS  
> cd NFTOOLS  
> wget https://github.com/juhussin/NFTOOLS/raw/master/NFTOOLS.tar.gz  
> tar -xvf NFTOOLS.tar.gz  

Latest version has been compiled with :

java version "1.7.0_45"  
Java(TM) SE Runtime Environment (build 1.7.0_45-b18)  
Java HotSpot(TM) 64-Bit Server VM (build 24.45-b08, mixed mode)  

To test, try :

> java -cp NucFamTools.jar Recomb

If you get an error, try using java version 1.7, or try recompiling, as explained below (Clone from github)


2) Clone from github

> git clone --recursive https://github.com/juhussin/NFTOOLS.git  
> cd NFTOOLS  

To test, try :

> java -cp NucFamTools.jar Recomb

If you get an error, try recompiling the code with your current version of java :

> cd NUCFAMTOOLS/  
> ./recompile.sh  
> cd ..  
> java -cp NucFamTools.jar Recomb  



USAGE
-----

The package is a collection of modules. The typical command line to run each module is :

> java -cp $path_to_directory/NFTOOLS/NucFamTools.jar module arguments [options]


Module PlinkToRecomb :
--------------------
This module convert .ped files in input file for Recomb module   
 usage : java -cp NucFamTools.jar PlinkToRecomb file.ped [-o outfile]

Mandatory argument (1):  
 file.ped : the 6 first columns of ped files are mandatory (see http://pngu.mgh.harvard.edu/~purcell/plink/data.shtml#ped), all delimited by tab characters. If generated with PLINK, use the --tab option in PLINK.  

Options :  
 -o outfile : name for output file. If not specified, default value is <file>\_recomb.ped  
 


Module Recomb :
-------------
This module produces reports on recombination events, errors and informative markers for each family in the dataset.  
 usage : java -cp NucFamTools.jar Recomb file.ped file.map [-o outfile -x -k N -err value(Kb) -valerr file.err]

Mandatory arguments (2):
 file.ped : an input file created by PlinkToRecomb module  
 file.map : a map file (see http://pngu.mgh.harvard.edu/~purcell/plink/data.shtml#map)  

Options :  
 -o outfile : name for output file. If not specified, default value is <file.ped>\_   
 -x : to perform the analysis on the X chromosome.  
 -k N : minimum number of informative markers allowed between recombination events. Below this number, double recombinants are ignored. Report files are generated for 0 to N. Default value is 2. For stringent analysis of crossovers, a value of 5 is suggested.  
 -err value : single markers causing a double recombinant are considered genotyping error, <value> is the max distance separating them from the nearest informative markers. Default value is to call as genotyping error all single markers causing a double recombinant.  
 -valerr file.err : a list of marker positions that will not be considered as errors (one position per line).  


Module RecombToEvents :
---------------------
This modules produces a list of events from the reports obtained with the Recomb module  
 usage : java -cp NucFamTools.jar RecombToEvents -l/-f file.recomb [-o outfile -fam excludeFamilies.txt -i individuals.txt -chrfile -printk]

Mandatory argument (1):  
 -l/-f file.recomb : a report file (-f) or a list of report files (-l) obtained from the Recomb module.

Options :  
 -o outfile : name for output file. If not specified, default value is <file.recomb>.events  
 -fam excludeFamilies.txt : a list of families to exclude (one family per line). The ID for the families are as specified in Recomb input file. Default is none.  
 -i individuals.txt : individuals to include in remaining families (one individual per line). Default is all.  
 -chrfile : takes the chromosome number from the names of the input files. The chromosome number has to be given after 'chr' and followed by '\_' (eg. *chr22\_*) in the filenames. By default, the chr ids in the output is the rank of each file in the list provided with -l option, or, for a single file, the chr id is 1 (-f).   
 -printk : prints in outfile the number of informative markers separating recombination events (before and after each one).  



Module BasicPerFamily :
---------------------
This module computes basic recombination statistics per family (informative markers, total number of maternal and paternal events, errors, ...)  
 usage : java -cp NucFamTools.jar BasicPerFamily nbfamily -l/-f file.recomb [-o outfile.res]

Mandatory arguments (2):  
 nbfamily : total number of families in input files  
 -l/-f file.recomb : a report file (-f) or a list of report files (-l) obtained from the Recomb module  

Options :  
 -o outfile : name for output file. If not specified, default value is <file.recomb>.res




RUNNING AN EXAMPLE
------------------

An example of a pipeline used to localise events in pedigree data is given in the perl script Recomb\_script.pl  
It performs the analysis on 2 chromosomes from 3 families (DATA/chr21\_3fams.ped and DATA/chr22\_3fams.ped)  

> perl Recomb\_script.pl




CONTACT INFORMATION
-------------------

For help or to report bugs/problems with the modules : please send email to ju.hussin@gmail.com
