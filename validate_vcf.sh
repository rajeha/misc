#!/usr/bin/env bash
# validate calls in VCF format
# requires:
#   yaha             (GregoryFaust/yaha)
#   fermi-lite       (lh3/fermi-lite)
#   is_*.pl scripts  (rajeha/validate-SV)
#   get_reads.sh     (rajeha/misc)
#   run_fermi.sh     (rajeha/misc)
#   fix_pairs.pl     (rajeha/fixpairs)

set -e

usage='cat svs.vcf | validate_vcf.sh <bam> <yaha_index> [read header pattern]'

if [ $# -ne 2 ]; then 
	echo $usage
	exit 1
fi

rh=${3:-@SR}

function unqdir() {
	if [ ! -d "$1" ]; then
		mkdir "$1" && cd "$1"
	else 
		newdir="${1}_"
		unqdir "$newdir"
	fi
}

bam=`realpath $1`
idx=`realpath $2`

THREADS=16

mkdir valsv.tmpdir
cd valsv.tmpdir

while read i; do
	if [[ $i =~ ^# ]]; then
		continue
	fi

	event=`echo "$i" | egrep -o 'SVTYPE=[^;]*' | sed 's/SVTYPE=//'`
		
	>&2 echo $event
	case "$event" in
	"DEL")
		sv=`echo "$i" | cut -f1,2`
		reg=`echo "$sv" | awk '{print $1, $2-1000, $2+1000}'| tr ' ' '-' | sed 's/-/:/'`
		dir=`echo "del_$sv" | sed 's/\t/_/'`
		unqdir "$dir"
		get_reads.sh $bam $reg r $rh &> .errgetreads
		run_fermi.sh r1.fq r2.fq &> .errfrm
		yaha -t $THREADS -x $idx -q frm.fa -o8 yh &> .erryaha
		val=`cat yh | is_del.pl $sv 2> .errisdel`
		if [ $val -eq 1 ]; then
			echo "$i"
		fi
		cd ../
		;;
	"INS")
		sv=`echo "$i" | cut -f1,2`
		reg=`echo "$sv" | awk '{print $1, $2-1000, $2+1000}'| tr ' ' '-' | sed 's/-/:/'`
		dir=`echo "ins_$sv" | sed 's/\t/_/'`
		unqdir "$dir"
		get_reads.sh $bam $reg r $rh &> .errgetreads
		run_fermi.sh r1.fq r2.fq &> .errfrm
		yaha -t $THREADS -x $idx -q frm.fa -o8 yh &> .erryaha
		val=`cat yh | is_ins.pl $sv 2> .errisins`
		if [ $val -eq 1 ]; then
			echo "$i"
		fi
		cd ../
		;;
	"INV")
		sv=`echo "$i" | cut -f1,2`
		reg=`echo "$sv" | awk '{print $1, $2-1000, $2+1000}'| tr ' ' '-' | sed 's/-/:/'`
		dir=`echo "inv_$sv" | sed 's/\t/_/'`
		unqdir "$dir"
		get_reads.sh $bam $reg r $rh &> .errgetreads
		run_fermi.sh r1.fq r2.fq &> .errfrm
		yaha -t $THREADS -x $idx -q frm.fa -o8 yh &> .erryaha
		val=`cat yh | is_inv.pl $sv 2> .errisinv`
		if [ $val -eq 1 ]; then
			echo "$i"
		fi
		cd ../
		;;
	"DUP" | "DUP:TANDEM")
		sv=`echo "$i" | cut -f1,2`
		reg=`echo "$sv" | awk '{print $1, $2-1000, $2+1000}'| tr ' ' '-' | sed 's/-/:/'`
		dir=`echo "dup_$sv" | sed 's/\t/_/'`
		unqdir "$dir"
		get_reads.sh $bam $reg r $rh &> .errgetreads
		run_fermi.sh r1.fq r2.fq &> .errfrm
		yaha -t $THREADS -x $idx -q frm.fa -FBS Y -o8 yh &> .erryaha
		val=`cat yh | is_dup.pl $sv 2> .errisdup`
		if [ $val -eq 1 ]; then
			echo "$i"
		fi
		cd ../
		;;
	"TRA" | "BND")	 
		sv=`echo "$i" | cut -f1,2,8 | sed 's/[^\t]*CHR2=//' | sed 's/;END=\([0-9]*\).*/\t\1/'` # delly
		#sv=`echo "$i" | cut -f1,2,5 | perl -pe 's/\S*?([A-Z]*):([0-9]*)[^\s]*/\1\t\2/'` # break-end format
		reg=`echo "$sv" | awk '{print $1, $2-1000, $2+1000}'| tr ' ' '-' | sed 's/-/:/'`
		dir=`echo "bnd_$sv" | sed 's/\t/_/g'`
		unqdir "$dir"
		get_reads.sh $bam $reg r $rh &> .errgetreads	
		run_fermi.sh r1.fq r2.fq &> .errfrm
		yaha -t $THREADS -x $idx -q frm.fa -FBS Y > yh.sam 2> .erryaha
		samtools view -bh yh.sam > yh.bam
		bedtools bamtobed -bedpe -i yh.bam > yh.bedpe 2> .errbtb
		val=`cat yh.bedpe | is_tra.pl $sv 2> .erristra`
		if [ $val -eq 1 ]; then
			echo "$i"
		fi
		cd ../
		;;
	esac
done 

cd ..
# comment out this line to keep intermediate files
rm -rf valsv.tmpdir
