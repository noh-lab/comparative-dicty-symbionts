#!/bin/bash

h0_list="../12cont_intact/codeml_genes.h0.txt"
sym_list="../12cont_intact/codeml_genes.sym.txt"
dic_list="../12cont_intact/codeml_genes.dic.txt"
red_list="../12cont_intact/codeml_genes.red.txt"

for gene in `cat "$h0_list"`; do
	awk '$1=="H0" {print "model",$0}' OFS="\t" "$gene"/"$gene".summary.txt | sed "s/model/$gene/" > model_summaries/"$gene".estimates.txt
done;

for gene in `cat "$sym_list"`; do
	awk '$1=="H1a" {print "model",$0}' OFS="\t" "$gene"/"$gene".summary.txt | sed "s/model/$gene/" > model_summaries/"$gene".estimates.txt
done;

for gene in `cat "$dic_list"`; do
	awk '$1=="H1b" {print "model",$0}' OFS="\t" "$gene"/"$gene".summary.txt | sed "s/model/$gene/" > model_summaries/"$gene".estimates.txt
done;

for gene in `cat "$red_list"`; do
	awk '$1=="H1c" {print "model",$0}' OFS="\t" "$gene"/"$gene".summary.txt | sed "s/model/$gene/" > model_summaries/"$gene".estimates.txt
done;

