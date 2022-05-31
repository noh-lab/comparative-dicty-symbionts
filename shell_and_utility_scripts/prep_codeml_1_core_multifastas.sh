#!/bin/bash

# prepare for codeml

# first, get core genes 

#singularity run -i -e -H /export/groups/snoh/snoh/12cont_intact/ /usr/local/singularity/roary_latest.sif query_pan_genome -a intersection *gff

#mkdir paraburk_core

#sed 's/:\s/\t/g' pan_genome_results > paraburk_core/core_genes.table
#cut -f 1 paraburk_core/core_genes.table > paraburk_core/core_genes.list

# now get faa and fna multifastas for each core gene

list="paraburk_core/core_genes.list"
#list="paraburk_core/test.list"
table="paraburk_core/core_genes.table" # the file where you keep your string names
#list="codeml/aroA.list"
#table="codeml/aroA.table"

for gene in `cat "$list"`; do
	agri=$(grep -w "$gene" "$table" | cut -f 2)
	bonn=$(grep -w "$gene" "$table" | cut -f 3)
	hayl=$(grep -w "$gene" "$table" | cut -f 4)
	cale=$(grep -w "$gene" "$table" | cut -f 5)
	fung=$(grep -w "$gene" "$table" | cut -f 6)
	mega=$(grep -w "$gene" "$table" | cut -f 7)
	phem=$(grep -w "$gene" "$table" | cut -f 8)
	phex=$(grep -w "$gene" "$table" | cut -f 9)
	phym=$(grep -w "$gene" "$table" | cut -f 10)
	phyt=$(grep -w "$gene" "$table" | cut -f 11)
	sart=$(grep -w "$gene" "$table" | cut -f 12)
	spre=$(grep -w "$gene" "$table" | cut -f 13)
	tera=$(grep -w "$gene" "$table" | cut -f 14)
	tere=$(grep -w "$gene" "$table" | cut -f 15)
	xeno=$(grep -w "$gene" "$table" | cut -f 16)	

	echo "processing aa sequences for $gene"
	samtools faidx ../roary_to_codeml/baqs159.prokka.v3.faa "$agri" > paraburk_core/"$gene".faa
	samtools faidx ../roary_to_codeml/bbqs859.prokka.v3.faa "$bonn" >> paraburk_core/"$gene".faa
	samtools faidx ../roary_to_codeml/bhqs11.prokka.v3.faa "$hayl" >> paraburk_core/"$gene".faa
	samtools faidx ../roary_to_codeml/pcale.prokka.v3.faa "$cale" >> paraburk_core/"$gene".faa
	samtools faidx ../roary_to_codeml/pfung.prokka.v3.faa "$fung" >> paraburk_core/"$gene".faa
	samtools faidx ../roary_to_codeml/pmega.prokka.v3.faa "$mega" >> paraburk_core/"$gene".faa
	samtools faidx ../roary_to_codeml/pphem.prokka.v3.faa "$phem" >> paraburk_core/"$gene".faa
	samtools faidx ../roary_to_codeml/pphex.prokka.v3.faa "$phex" >> paraburk_core/"$gene".faa
	samtools faidx ../roary_to_codeml/pphym.prokka.v3.faa "$phym" >> paraburk_core/"$gene".faa
	samtools faidx ../roary_to_codeml/pphyt.prokka.v3.faa "$phyt" >> paraburk_core/"$gene".faa
	samtools faidx ../roary_to_codeml/psart.prokka.v3.faa "$sart" >> paraburk_core/"$gene".faa
	samtools faidx ../roary_to_codeml/pspre.prokka.v3.faa "$spre" >> paraburk_core/"$gene".faa
	samtools faidx ../roary_to_codeml/ptera.prokka.v3.faa "$tera" >> paraburk_core/"$gene".faa
	samtools faidx ../roary_to_codeml/ptere.prokka.v3.faa "$tere" >> paraburk_core/"$gene".faa
	samtools faidx ../roary_to_codeml/pxeno.prokka.v3.faa "$xeno" >> paraburk_core/"$gene".faa
	
	echo "getting nucleotide sequences for $gene"
	samtools faidx ../roary_to_codeml/baqs159.prokka.v3.ffn "$agri" > paraburk_core/"$gene".ffn
	samtools faidx ../roary_to_codeml/bbqs859.prokka.v3.ffn "$bonn" >> paraburk_core/"$gene".ffn
	samtools faidx ../roary_to_codeml/bhqs11.prokka.v3.ffn "$hayl" >> paraburk_core/"$gene".ffn
	samtools faidx ../roary_to_codeml/pcale.prokka.v3.ffn "$cale" >> paraburk_core/"$gene".ffn
	samtools faidx ../roary_to_codeml/pfung.prokka.v3.ffn "$fung" >> paraburk_core/"$gene".ffn
	samtools faidx ../roary_to_codeml/pmega.prokka.v3.ffn "$mega" >> paraburk_core/"$gene".ffn
	samtools faidx ../roary_to_codeml/pphem.prokka.v3.ffn "$phem" >> paraburk_core/"$gene".ffn
	samtools faidx ../roary_to_codeml/pphex.prokka.v3.ffn "$phex" >> paraburk_core/"$gene".ffn
	samtools faidx ../roary_to_codeml/pphym.prokka.v3.ffn "$phym" >> paraburk_core/"$gene".ffn
	samtools faidx ../roary_to_codeml/pphyt.prokka.v3.ffn "$phyt" >> paraburk_core/"$gene".ffn
	samtools faidx ../roary_to_codeml/psart.prokka.v3.ffn "$sart" >> paraburk_core/"$gene".ffn
	samtools faidx ../roary_to_codeml/pspre.prokka.v3.ffn "$spre" >> paraburk_core/"$gene".ffn
	samtools faidx ../roary_to_codeml/ptera.prokka.v3.ffn "$tera" >> paraburk_core/"$gene".ffn
	samtools faidx ../roary_to_codeml/ptere.prokka.v3.ffn "$tere" >> paraburk_core/"$gene".ffn
	samtools faidx ../roary_to_codeml/pxeno.prokka.v3.ffn "$xeno" >> paraburk_core/"$gene".ffn

done
