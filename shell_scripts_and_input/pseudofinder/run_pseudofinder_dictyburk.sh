#!/bin/bash

pseudofinder.py annotate --diamond --skip_makedb -g /research/snoh/data/colby/comparative_burk/prokka/bbqs859/PROKKA_06192020.gbk -db /research/snoh/diamonddb/refseq_protein_nr.dmnd -op bbqs859

pseudofinder.py annotate --diamond --skip_makedb -g /research/snoh/data/colby/comparative_burk/prokka/bhqs11/PROKKA_06192020.gbk -db /research/snoh/diamonddb/refseq_protein_nr.dmnd -op bhqs11

pseudofinder.py annotate --diamond --skip_makedb -g /research/snoh/data/colby/comparative_burk/prokka/baqs159/PROKKA_06192020.gbk -db /research/snoh/diamonddb/refseq_protein_nr.dmnd -op baqs159


grep PAGRI baqs159_pseudos.gff | sort -k1,1 -k4,4n | awk -F "\t" '{split($9,a,";"); split(a[3],b,"="); split(a[1],c,":"); print b[2]}' OFS="\t" > pagri_pseudofinder.txt

grep PBONN bbqs859_pseudos.gff | sort -k1,1 -k4,4n | awk -F "\t" '{split($9,a,";"); split(a[3],b,"="); split(a[1],c,":"); print b[2]}' OFS="\t" > pbonn_pseudofinder.txt

grep PHAYL bhqs11_pseudos.gff | sort -k1,1 -k4,4n | awk -F "\t" '{split($9,a,";"); split(a[3],b,"="); split(a[1],c,":"); print b[2]}' OFS="\t" > phayl_pseudofinder.txt

sed 's/,/\n/g' pagri_pseudofinder.txt | grep PAGRI | sort > temp
mv temp pagri_pseudofinder.txt

sed 's/,/\n/g' pbonn_pseudofinder.txt | grep PBONN | sort > temp
mv temp pbonn_pseudofinder.txt

sed 's/,/\n/g' phayl_pseudofinder.txt | grep PHAYL | sort > temp
mv temp phayl_pseudofinder.txt

