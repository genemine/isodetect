#!/bin/sh
t=$1
usearch9.2.64_i86linux32 -sortbylength $1 -fastaout $s1.clu.sort.fa 
usearch9.2.64_i86linux32 -cluster_fast $s1.clu.sort.fa -consout uclu.$1.fax -id 0.75 
echo date >> flog.log
