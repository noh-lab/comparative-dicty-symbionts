#!/bin/bash

# remove pseudogenes before running roary

awk '{print $1}' baqs159.prokka.v3.ffn.fai > pagri_genes.txt
awk 'NR==FNR {a[$1]; next} !($1 in a) {print}' ../pseudofinder/pagri_pseudofinder.txt pagri_genes.txt > pagri_intact.txt
grep -f pagri_intact.txt baqs159.prokka.v3.gff > 12cont_intact/baqs159.prokka.intact.gff

awk '{print $1}' bbqs859.prokka.v3.ffn.fai > pbonn_genes.txt
awk 'NR==FNR {a[$1]; next} !($1 in a) {print}' ../pseudofinder/pbonn_pseudofinder.txt pbonn_genes.txt > pbonn_intact.txt
grep -f pbonn_intact.txt bbqs859.prokka.v3.gff > 12cont_intact/bbqs859.prokka.intact.gff

awk '{print $1}' bhqs11.prokka.v3.ffn.fai > phayl_genes.txt
awk 'NR==FNR {a[$1]; next} !($1 in a) {print}' ../pseudofinder/phayl_pseudofinder.txt phayl_genes.txt > phayl_intact.txt
grep -f phayl_intact.txt bhqs11.prokka.v3.gff > 12cont_intact/bhqs11.prokka.intact.gff

awk '{print $1}' pcale.prokka.v3.ffn.fai > pcale_genes.txt
awk 'NR==FNR {a[$1]; next} !($1 in a) {print}' ../pseudofinder/pcale_pseudofinder.txt pcale_genes.txt > pcale_intact.txt
grep -f pcale_intact.txt pcale.prokka.v3.gff > 12cont_intact/pcale.prokka.intact.gff

awk '{print $1}' pfung.prokka.v3.ffn.fai > pfung_genes.txt
awk 'NR==FNR {a[$1]; next} !($1 in a) {print}' ../pseudofinder/pfung_pseudofinder.txt pfung_genes.txt > pfung_intact.txt
grep -f pfung_intact.txt pfung.prokka.v3.gff > 12cont_intact/pfung.prokka.intact.gff

awk '{print $1}' pmega.prokka.v3.ffn.fai > pmega_genes.txt
awk 'NR==FNR {a[$1]; next} !($1 in a) {print}' ../pseudofinder/pmega_pseudofinder.txt pmega_genes.txt > pmega_intact.txt
grep -f pmega_intact.txt pmega.prokka.v3.gff > 12cont_intact/pmega.prokka.intact.gff

awk '{print $1}' pphem.prokka.v3.ffn.fai > pphem_genes.txt
awk 'NR==FNR {a[$1]; next} !($1 in a) {print}' ../pseudofinder/pphem_pseudofinder.txt pphem_genes.txt > pphem_intact.txt
grep -f pphem_intact.txt pphem.prokka.v3.gff > 12cont_intact/pphem.prokka.intact.gff

awk '{print $1}' pphex.prokka.v3.ffn.fai > pphex_genes.txt
awk 'NR==FNR {a[$1]; next} !($1 in a) {print}' ../pseudofinder/pphex_pseudofinder.txt pphex_genes.txt > pphex_intact.txt
grep -f pphex_intact.txt pphex.prokka.v3.gff > 12cont_intact/pphex.prokka.intact.gff

awk '{print $1}' pphym.prokka.v3.ffn.fai > pphym_genes.txt
awk 'NR==FNR {a[$1]; next} !($1 in a) {print}' ../pseudofinder/pphym_pseudofinder.txt pphym_genes.txt > pphym_intact.txt
grep -f pphym_intact.txt pphym.prokka.v3.gff > 12cont_intact/pphym.prokka.intact.gff

awk '{print $1}' pphyt.prokka.v3.ffn.fai > pphyt_genes.txt
awk 'NR==FNR {a[$1]; next} !($1 in a) {print}' ../pseudofinder/pphyt_pseudofinder.txt pphyt_genes.txt > pphyt_intact.txt
grep -f pphyt_intact.txt pphyt.prokka.v3.gff > 12cont_intact/pphyt.prokka.intact.gff

awk '{print $1}' psart.prokka.v3.ffn.fai > psart_genes.txt
awk 'NR==FNR {a[$1]; next} !($1 in a) {print}' ../pseudofinder/psart_pseudofinder.txt psart_genes.txt > psart_intact.txt
grep -f psart_intact.txt psart.prokka.v3.gff > 12cont_intact/psart.prokka.intact.gff

awk '{print $1}' pspre.prokka.v3.ffn.fai > pspre_genes.txt
awk 'NR==FNR {a[$1]; next} !($1 in a) {print}' ../pseudofinder/pspre_pseudofinder.txt pspre_genes.txt > pspre_intact.txt
grep -f pspre_intact.txt pspre.prokka.v3.gff > 12cont_intact/pspre.prokka.intact.gff

awk '{print $1}' ptera.prokka.v3.ffn.fai > ptera_genes.txt
awk 'NR==FNR {a[$1]; next} !($1 in a) {print}' ../pseudofinder/ptera_pseudofinder.txt ptera_genes.txt > ptera_intact.txt
grep -f ptera_intact.txt ptera.prokka.v3.gff > 12cont_intact/ptera.prokka.intact.gff

awk '{print $1}' ptere.prokka.v3.ffn.fai > ptere_genes.txt
awk 'NR==FNR {a[$1]; next} !($1 in a) {print}' ../pseudofinder/ptere_pseudofinder.txt ptere_genes.txt > ptere_intact.txt
grep -f ptere_intact.txt ptere.prokka.v3.gff > 12cont_intact/ptere.prokka.intact.gff

awk '{print $1}' pxeno.prokka.v3.ffn.fai > pxeno_genes.txt
awk 'NR==FNR {a[$1]; next} !($1 in a) {print}' ../pseudofinder/pxeno_pseudofinder.txt pxeno_genes.txt > pxeno_intact.txt
grep -f pxeno_intact.txt pxeno.prokka.v3.gff > 12cont_intact/pxeno.prokka.intact.gff

rm *_genes.txt *_intact.txt
