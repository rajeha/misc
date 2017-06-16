#!/usr/bin/env bash
# get_reads.sh <bam> [region]
# TODO: add argument for reads header
if [ $# -lt 2 ]; then 
	echo "get_reads.sh <bam> <region> [output prefix]"
	exit
fi

opref=${3:-r}

samtools view -h $1 $2 | samtools fastq -1 r1.tmp -2 r2.tmp -

perl -pi -e 's/\@SR.+\S/$&\/1/' r1.tmp 
perl -pi -e 's/\@SR.+\S/$&\/2/' r2.tmp 

fix_pairs.pl r1.tmp r2.tmp $opref

rm -f *.tmp ${opref}_unpaired.fq

