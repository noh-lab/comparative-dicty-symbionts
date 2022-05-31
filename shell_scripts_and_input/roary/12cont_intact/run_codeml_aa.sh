#!/bin/bash

list="paraburk_core/core_genes.aa"

for gene in `cat "$list"`; do
	sed '2,$s/[0-9_]//g' ../roary_to_codeml/"$gene"/"$gene".muscle.pal2nal > ../roary_to_codeml/"$gene"/"$gene".clean.pal2nal
	sed "s/INFILE/$gene.clean.pal2nal/" codeml.ctl | sed "s/OUTFILE/$gene.codeml.m0.txt/" > ../roary_to_codeml/"$gene"/codeml.ctl
	(cd ../roary_to_codeml/"$gene"; codeml)
	mv ../roary_to_codeml/"$gene"/rst1 ../roary_to_codeml/"$gene"/"$gene".codeml.rst1
	sed "s/INFILE/$gene.clean.pal2nal/" codeml.m2.h2a | sed "s/OUTFILE/$gene.codeml.m2a.txt/" > ../roary_to_codeml/"$gene"/codeml.ctl
	(cd ../roary_to_codeml/"$gene"; codeml)
	cat ../roary_to_codeml/"$gene"/rst1 >> ../roary_to_codeml/"$gene"/"$gene".codeml.rst1
	sed "s/INFILE/$gene.clean.pal2nal/" codeml.m2.h2b | sed "s/OUTFILE/$gene.codeml.m2b.txt/" > ../roary_to_codeml/"$gene"/codeml.ctl
	(cd ../roary_to_codeml/"$gene"; codeml)
	cat ../roary_to_codeml/"$gene"/rst1 >> ../roary_to_codeml/"$gene"/"$gene".codeml.rst1
	sed "s/INFILE/$gene.clean.pal2nal/" codeml.m2.h2c | sed "s/OUTFILE/$gene.codeml.m2c.txt/" > ../roary_to_codeml/"$gene"/codeml.ctl
	(cd ../roary_to_codeml/"$gene"; codeml)
	cat ../roary_to_codeml/"$gene"/rst1 >> ../roary_to_codeml/"$gene"/"$gene".codeml.rst1
done;

