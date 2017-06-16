#!/usr/bin/env bash

usage='run_ferm.sh <r1.fq> <r2.fq> [output prefix]'

if [ $# -lt 2 ]; then
	echo $usage 
	exit;
fi

opref=${3:-frm}

fml-asm $1 $2 > frm.tmp1
sed -n '1~4s/^@/>/p;2~4p' frm.tmp1 > frm.tmp2
cat frm.tmp2 | sed 's/>\([^\t]*\).*/>\1/' > ${opref}.fa 

rm -f frm.tmp1 frm.tmp2
