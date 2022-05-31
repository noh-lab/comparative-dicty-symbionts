#!/bin/bash

# muscle alignment and pal2nal conversion

list="paraburk_core/core_genes.list"

for gene in `cat "$list"`; do
	mkdir ../roary_to_codeml/"$gene"
	muscle -in paraburk_core/"$gene".faa -out "$gene".temp.muscle
	grep ">" paraburk_core/"$gene".ffn | tr '\n' ' ' | sed 's/>//g' > "$gene".temp.order
	samtools faidx "$gene".temp.muscle
	samtools faidx "$gene".temp.muscle `cat "$gene".temp.order` > "$gene".temp.aln
	/research/snoh/data/colby/comparative_burk/pal2nal.pl "$gene".temp.aln paraburk_core/"$gene".ffn -output paml -nogap -nomismatch -codontable 11 > ../roary_to_codeml/"$gene"/"$gene".muscle.pal2nal
	rm "$gene".temp*
done;

