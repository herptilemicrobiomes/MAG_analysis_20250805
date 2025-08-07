#!/usr/bin/bash -l

for a in $(ls -d /bigdata/stajichlab/shared/projects/Herptile/Metagenome/PhaseGenomics/*/*)
do
    UHM=$(basename $a)
    mkdir -p input/$UHM/bins
    for b in $(ls $a/all/microbial_genomes/*_clusters/*.fasta)
    do
	binname=$(basename $b .fasta | perl -p -e 's/bin_/bin./')
#	echo " $b -> input/$UHM/bins/$UHM.$binname.fa"
	ln -s $b input/$UHM/bins/$UHM.$binname.fa
    done
done
