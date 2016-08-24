#!/usr/bin/perl

system "mkdir DATA/ERRORS";
system "mkdir DATA/INFOMARKERS";

for $i (21 .. 22)
{
    #PLINK TO RECOMB
    system "java -Xmx512m -cp NucFamTools.jar PlinkToRecomb DATA/chr".$i."_3fams.ped -o DATA/chr".$i."_3fams_inrecomb.ped";
    
    #RECOMB
    system "java -Xmx512m -cp NucFamTools.jar Recomb DATA/chr".$i."_3fams_inrecomb.ped DATA/chr".$i."_3fams.map -o DATA/chr".$i."_3fams -k 5";

}

system "mv DATA/*info* DATA/INFOMARKERS";
system "mv DATA/*allErrors* DATA/ERRORS";

#Keeping only results for k=0, k=2 and k=5

for $i (1,3,4)
{
    system "rm -f DATA/chr*k".$i.".recomb";
}

for $i (0,2,5)
{
    #RECOMB TO EVENTS LIST
    system "ls DATA/chr*k".$i."* > DATA/listchr_k".$i.".txt";
    system "java -Xmx512m -cp NucFamTools.jar RecombToEvents -l DATA/listchr_k".$i.".txt -chrfile";
}

for $i (0,2,5)
{
    #RECOMB TO SUMMARY STATS
    system "java -Xmx512m -cp NucFamTools.jar BasicPerFamily 3 -l DATA/listchr_k".$i.".txt ";
}