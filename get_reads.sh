#!/usr/bin/env bash
# Because of advanced, sophisticated optional arguemnt handling,
# the user must set arg 3 to be able to set arg 4.
 
usage="get_reads.sh <bam> <region> [output prefix] [read header pattern]"

if [ $# -lt 2 ]; then 
	echo $usage; 
	exit 1
fi

opref=${3:-r}
rh=${4:-@SR} # default read header starts with @SR

samtools view -h $1 $2 | samtools fastq -1 r1.tmp -2 r2.tmp -

perl -spi -e 's/$rh.+\S/$&\/1/' -- -rh="$rh" r1.tmp 
perl -spi -e 's/$rh.+\S/$&\/2/' -- -rh="$rh" r2.tmp 

fix_pairs.pl r1.tmp r2.tmp $opref

rm -f *.tmp ${opref}_unpaired.fq

