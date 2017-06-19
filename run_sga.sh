#!/usr/bin/env bash

usage='runsga.sh <r1.fq> <r2.fq> [output prefix]'

if [ $# -lt 2 ]; then
	echo $usage 
	exit 1
fi

opref=${3:-wlk}

THREADS=16
KMER=13

echo -e "Running SGA:"
echo -e "[SGA] Preprocessing..."
sga preprocess -p 1 -o sga_pre $1 $2 
echo -e "[SGA] Indexing..."
sga index -t $THREADS --no-reverse sga_pre
echo -e "[SGA] Correcting..."
sga correct -t $THREADS -k $KMER -o sga_corr sga_pre 
echo -e "[SGA] Indexing again..."
sga index -t $THREADS sga_corr
echo -e "[SGA] Removing duplicates..."
sga rmdup -t $THREADS -o sga_rmd sga_corr
echo -e "[SGA] Structuring string graph..."
sga overlap -t $THREADS -m 15 sga_rmd.fa
echo -e "[SGA] Assembling contigs..."
sga assemble -m 20 -d 0 -g 0 -b 0 -l 100 -o sga_asm sga_rmd.asqg.gz
echo -e "[SGA] Finding walks..."
sga walk -d 10000 --component-walks --longest-n 10000 -o sga_wlk sga_asm-graph.asqg.gz

mv sga_wlk ${opref}.fa 
rm -rf sga_*
