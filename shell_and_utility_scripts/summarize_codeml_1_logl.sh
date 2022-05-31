#!/bin/bash
# this takes key parameters from the rst files and organizes them into a table
# then logLikelihoods are taken so they can be concatenated into a single file for likelihood ratio tests

# e.g. what the lines of an rst file looks like (m2c)
#	15	243	227	0.6251	0.066	0.080	0.087	0.075	0.046	0.081	0.103	0.064	0.038	0.018	0.028	0.025	0.024	0.101	0.032	0.101	0.083	0.073	0.039	0.066	0.102	0.088	0.016	0.109	0.102	0.080	0.075	0.024	
#(kappa)2.588	(omega)0.000	0.069	(lnL)-4873.801	(tree dN)0.0130	(tree dS)1.3041


#cd roary_to_codeml

list="../12cont_intact/paraburk_core/core_genes.list"
#list="../12cont_intact/paraburk_core/aroA.list"

for gene in `cat "$list"`; do
	sed 's/\s+/\t/g' "$gene"/"$gene".codeml.rst1 | sed 's/^\t//g' | awk 'NR==1{print "model","kappa","w0","w1","lnL\nH0", $33, $34, "NA", $35}; NR==2{print "H1a", $33, $34, $35, $36}; NR==3{print "H1b", $33, $34, $35, $36}; NR==4{print "H1c", $33, $34, $35, $36}' OFS="\t" - > "$gene"/"$gene".summary.txt
	cut -f 5 "$gene"/"$gene".summary.txt | ./transpose.awk - | sed "s/lnL/$gene/" > model_summaries/"$gene".summary.txt
done;
