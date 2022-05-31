#!/bin/bash
# type 3 - sctJ, sctN, sctV
# type 4 - virb4, virB8, virB6 --> 6 and 8 may not work because not all have them (Guglielmini et al 2014, especially Fig 6)
# type 6 - tssB, tssC\|iglB, tssF

for gene in `echo "virB6"`; do
	agri=$(grep "$gene" all.manual_ss_components.txt | grep AGRI | cut -f 4)
	bonn=$(grep "$gene" all.manual_ss_components.txt | grep BONN | cut -f 4)
	hayl=$(grep "$gene" all.manual_ss_components.txt | grep HAYL | cut -f 4)
	mall=$(grep "$gene" all.manual_ss_components.txt | grep BMA | cut -f 4)
	pseu=$(grep "$gene" all.manual_ss_components.txt | grep BPS | cut -f 4)
	cale=$(grep "$gene" all.manual_ss_components.txt | grep CALE | cut -f 4)
	fung=$(grep "$gene" all.manual_ss_components.txt | grep FUNG | cut -f 4)
	mega=$(grep "$gene" all.manual_ss_components.txt | grep MEGA | cut -f 4)
	phem=$(grep "$gene" all.manual_ss_components.txt | grep PHEM | cut -f 4)
	phex=$(grep "$gene" all.manual_ss_components.txt | grep PHEX | cut -f 4)
	phym=$(grep "$gene" all.manual_ss_components.txt | grep PHYM | cut -f 4)
	phyt=$(grep "$gene" all.manual_ss_components.txt | grep PHYT | cut -f 4)
	sart=$(grep "$gene" all.manual_ss_components.txt | grep SART | cut -f 4)
	spre=$(grep "$gene" all.manual_ss_components.txt | grep SPRE | cut -f 4)
	tera=$(grep "$gene" all.manual_ss_components.txt | grep TERA | cut -f 4)
	tere=$(grep "$gene" all.manual_ss_components.txt | grep TERE | cut -f 4)
	xeno=$(grep "$gene" all.manual_ss_components.txt | grep XENO | cut -f 4)
	echo "getting aa sequences for $gene"
	samtools faidx ../roary/baqs159.prokka.v3.faa `echo "$agri"` > "$gene"_manual.faa
	samtools faidx ../roary/bbqs859.prokka.v3.faa `echo "$bonn"` >> "$gene"_manual.faa
	samtools faidx ../roary/bhqs11.prokka.v3.faa `echo "$hayl"` >> "$gene"_manual.faa
	samtools faidx Burkholderia_mallei_ATCC_23344_131.faa `echo "$mall"` >> "$gene"_manual.faa
	samtools faidx Burkholderia_pseudomallei_K96243_132.faa `echo "$pseu"` >> "$gene"_manual.faa
	samtools faidx ../roary/pcale.prokka.v3.faa `echo "$cale"` >> "$gene"_manual.faa
	samtools faidx ../roary/pfung.prokka.v3.faa `echo "$fung"` >> "$gene"_manual.faa
	samtools faidx ../roary/pmega.prokka.v3.faa `echo "$mega"` >> "$gene"_manual.faa
	samtools faidx ../roary/pphem.prokka.v3.faa `echo "$phem"` >> "$gene"_manual.faa
	samtools faidx ../roary/pphex.prokka.v3.faa `echo "$phex"` >> "$gene"_manual.faa
	samtools faidx ../roary/pphym.prokka.v3.faa `echo "$phym"` >> "$gene"_manual.faa
	samtools faidx ../roary/pphyt.prokka.v3.faa `echo "$phyt"` >> "$gene"_manual.faa
	samtools faidx ../roary/psart.prokka.v3.faa `echo "$sart"` >> "$gene"_manual.faa
	samtools faidx ../roary/pspre.prokka.v3.faa `echo "$spre"` >> "$gene"_manual.faa
	samtools faidx ../roary/ptera.prokka.v3.faa `echo "$tera"` >> "$gene"_manual.faa
	samtools faidx ../roary/ptere.prokka.v3.faa `echo "$tere"` >> "$gene"_manual.faa
	samtools faidx ../roary/pxeno.prokka.v3.faa `echo "$xeno"` >> "$gene"_manual.faa
done
