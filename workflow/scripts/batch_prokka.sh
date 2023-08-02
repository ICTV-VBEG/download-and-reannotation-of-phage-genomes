#!/bin/bash
for file in *.fasta; do
	tag=${file%.fasta};
	string=${tag##*_};
	echo $string	
	prokka --hmms /bioinformatics/prokka/db/hmm/all_phrogs_ann.hmm --locustag $string --outdir $tag --prefix $tag --gcode 11 --cpus 8 --rfam $file
done