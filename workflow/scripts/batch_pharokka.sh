#!/bin/bash
for file in *.fasta; do
	tag=${file%.fasta};
	string=${tag##*_};
	echo $string	
	pharokka -i $file -l $string -o $tag -p $tag -t 8
done