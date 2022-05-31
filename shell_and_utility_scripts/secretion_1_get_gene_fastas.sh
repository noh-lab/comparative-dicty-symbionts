#!/bin/bash
for gene in `echo "virB9"`; do
	agri=$(grep "$gene" baqs159.ss_components.txt | cut -f 4)
	bonn=$(grep "$gene" bbqs859.ss_components.txt | cut -f 4)
	hayl=$(grep "$gene" bhqs11.ss_components.txt | cut -f 4)
	mall=$(grep "$gene" bmall.ss_components.txt | cut -f 4)
	pseu=$(grep "$gene" bpseu.ss_components.txt | cut -f 4)
	cale=$(grep "$gene" pcale.ss_components.txt | cut -f 4)
	fung=$(grep "$gene" pfung.ss_components.txt | cut -f 4)
	mega=$(grep "$gene" pmega.ss_components.txt | cut -f 4)
	phem=$(grep "$gene" pphem.ss_components.txt | cut -f 4)
	phex=$(grep "$gene" pphex.ss_components.txt | cut -f 4)
	phym=$(grep "$gene" pphym.ss_components.txt | cut -f 4)
	phyt=$(grep "$gene" pphyt.ss_components.txt | cut -f 4)
	sart=$(grep "$gene" psart.ss_components.txt | cut -f 4)
	spre=$(grep "$gene" pspre.ss_components.txt | cut -f 4)
	tera=$(grep "$gene" ptera.ss_components.txt | cut -f 4)
	tere=$(grep "$gene" ptere.ss_components.txt | cut -f 4)
	xeno=$(grep "$gene" pxeno.ss_components.txt | cut -f 4)
	echo "getting aa sequences for $gene"
	samtools faidx ../roary/baqs159.prokka.v3.faa `echo "$agri"` > "$gene".faa
	samtools faidx ../roary/bbqs859.prokka.v3.faa `echo "$bonn"` >> "$gene".faa
	samtools faidx ../roary/bhqs11.prokka.v3.faa `echo "$hayl"` >> "$gene".faa
	samtools faidx Burkholderia_mallei_ATCC_23344_131.faa `echo "$mall"` >> "$gene".faa
	samtools faidx Burkholderia_pseudomallei_K96243_132.faa `echo "$pseu"` >> "$gene".faa
	samtools faidx ../roary/pcale.prokka.v3.faa `echo "$cale"` >> "$gene".faa
	samtools faidx ../roary/pfung.prokka.v3.faa `echo "$fung"` >> "$gene".faa
	samtools faidx ../roary/pmega.prokka.v3.faa `echo "$mega"` >> "$gene".faa
	samtools faidx ../roary/pphem.prokka.v3.faa `echo "$phem"` >> "$gene".faa
	samtools faidx ../roary/pphex.prokka.v3.faa `echo "$phex"` >> "$gene".faa
	samtools faidx ../roary/pphym.prokka.v3.faa `echo "$phym"` >> "$gene".faa
	samtools faidx ../roary/pphyt.prokka.v3.faa `echo "$phyt"` >> "$gene".faa
	samtools faidx ../roary/psart.prokka.v3.faa `echo "$sart"` >> "$gene".faa
	samtools faidx ../roary/pspre.prokka.v3.faa `echo "$spre"` >> "$gene".faa
	samtools faidx ../roary/ptera.prokka.v3.faa `echo "$tera"` >> "$gene".faa
	samtools faidx ../roary/ptere.prokka.v3.faa `echo "$tere"` >> "$gene".faa
	samtools faidx ../roary/pxeno.prokka.v3.faa `echo "$xeno"` >> "$gene".faa
done
