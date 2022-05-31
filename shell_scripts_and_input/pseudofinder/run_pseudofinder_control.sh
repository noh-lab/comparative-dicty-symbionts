#!/bin/bash

pseudofinder.py annotate --diamond --skip_makedb -g /research/snoh/data/colby/comparative_burk/prokka/pfung/PROKKA_06192020.gbk -db /research/snoh/diamonddb/refseq_protein_nr.dmnd -op pfung

pseudofinder.py annotate --diamond --skip_makedb -g /research/snoh/data/colby/comparative_burk/prokka/ptere/PROKKA_06192020.gbk -db /research/snoh/diamonddb/refseq_protein_nr.dmnd -op ptere

pseudofinder.py annotate --diamond --skip_makedb -g /research/snoh/data/colby/comparative_burk/prokka/pspre/PROKKA_06192020.gbk -db /research/snoh/diamonddb/refseq_protein_nr.dmnd -op pspre

pseudofinder.py annotate --diamond --skip_makedb -g /research/snoh/data/colby/comparative_burk/prokka/pxeno/PROKKA_06192020.gbk -db /research/snoh/diamonddb/refseq_protein_nr.dmnd -op pxeno

pseudofinder.py annotate --diamond --skip_makedb -g /research/snoh/data/colby/comparative_burk/prokka/pcale/PROKKA_06192020.gbk -db /research/snoh/diamonddb/refseq_protein_nr.dmnd -op pcale

pseudofinder.py annotate --diamond --skip_makedb -g /research/snoh/data/colby/comparative_burk/prokka/pmega/PROKKA_06192020.gbk -db /research/snoh/diamonddb/refseq_protein_nr.dmnd -op pmega

pseudofinder.py annotate --diamond --skip_makedb -g /research/snoh/data/colby/comparative_burk/prokka/pphem/PROKKA_06192020.gbk -db /research/snoh/diamonddb/refseq_protein_nr.dmnd -op pphem

pseudofinder.py annotate --diamond --skip_makedb -g /research/snoh/data/colby/comparative_burk/prokka/pphex/PROKKA_06192020.gbk -db /research/snoh/diamonddb/refseq_protein_nr.dmnd -op pphex

pseudofinder.py annotate --diamond --skip_makedb -g /research/snoh/data/colby/comparative_burk/prokka/pphym/PROKKA_06192020.gbk -db /research/snoh/diamonddb/refseq_protein_nr.dmnd -op pphym

pseudofinder.py annotate --diamond --skip_makedb -g /research/snoh/data/colby/comparative_burk/prokka/pphyt/PROKKA_06192020.gbk -db /research/snoh/diamonddb/refseq_protein_nr.dmnd -op pphyt

pseudofinder.py annotate --diamond --skip_makedb -g /research/snoh/data/colby/comparative_burk/prokka/psart/PROKKA_06192020.gbk -db /research/snoh/diamonddb/refseq_protein_nr.dmnd -op psart

pseudofinder.py annotate --diamond --skip_makedb -g /research/snoh/data/colby/comparative_burk/prokka/ptera/PROKKA_06192020.gbk -db /research/snoh/diamonddb/refseq_protein_nr.dmnd -op ptera

grep PFUNG pfung_pseudos.gff | sort -k1,1 -k4,4n | awk -F "\t" '{split($9,a,";"); split(a[3],b,"="); split(a[1],c,":"); print b[2]}' OFS="\t" > pfung_pseudofinder.txt

grep PSPRE pspre_pseudos.gff | sort -k1,1 -k4,4n | awk -F "\t" '{split($9,a,";"); split(a[3],b,"="); split(a[1],c,":"); print b[2]}' OFS="\t" > pspre_pseudofinder.txt

grep PTERE ptere_pseudos.gff | sort -k1,1 -k4,4n | awk -F "\t" '{split($9,a,";"); split(a[3],b,"="); split(a[1],c,":"); print b[2]}' OFS="\t" > ptere_pseudofinder.txt

grep PXENO pxeno_pseudos.gff | sort -k1,1 -k4,4n | awk -F "\t" '{split($9,a,";"); split(a[3],b,"="); split(a[1],c,":"); print b[2]}' OFS="\t" > pxeno_pseudofinder.txt

grep PCALE pcale_pseudos.gff | sort -k1,1 -k4,4n | awk -F "\t" '{split($9,a,";"); split(a[3],b,"="); split(a[1],c,":"); print b[2]}' OFS="\t" > pcale_pseudofinder.txt

grep PMEGA pmega_pseudos.gff | sort -k1,1 -k4,4n | awk -F "\t" '{split($9,a,";"); split(a[3],b,"="); split(a[1],c,":"); print b[2]}' OFS="\t" > pmega_pseudofinder.txt

grep PPHEM pphem_pseudos.gff | sort -k1,1 -k4,4n | awk -F "\t" '{split($9,a,";"); split(a[3],b,"="); split(a[1],c,":"); print b[2]}' OFS="\t" > pphem_pseudofinder.txt

grep PPHEX pphex_pseudos.gff | sort -k1,1 -k4,4n | awk -F "\t" '{split($9,a,";"); split(a[3],b,"="); split(a[1],c,":"); print b[2]}' OFS="\t" > pphex_pseudofinder.txt

grep PPHYM pphym_pseudos.gff | sort -k1,1 -k4,4n | awk -F "\t" '{split($9,a,";"); split(a[3],b,"="); split(a[1],c,":"); print b[2]}' OFS="\t" > pphym_pseudofinder.txt

grep PPHYT pphyt_pseudos.gff | sort -k1,1 -k4,4n | awk -F "\t" '{split($9,a,";"); split(a[3],b,"="); split(a[1],c,":"); print b[2]}' OFS="\t" > pphyt_pseudofinder.txt

grep PSART psart_pseudos.gff | sort -k1,1 -k4,4n | awk -F "\t" '{split($9,a,";"); split(a[3],b,"="); split(a[1],c,":"); print b[2]}' OFS="\t" > psart_pseudofinder.txt

grep PTERA ptera_pseudos.gff | sort -k1,1 -k4,4n | awk -F "\t" '{split($9,a,";"); split(a[3],b,"="); split(a[1],c,":"); print b[2]}' OFS="\t" > ptera_pseudofinder.txt


sed 's/,/\n/g' pcale_pseudofinder.txt | grep PCALE | sort > temp
mv temp pcale_pseudofinder.txt

sed 's/,/\n/g' pfung_pseudofinder.txt | grep PFUNG | sort > temp
mv temp pfung_pseudofinder.txt

sed 's/,/\n/g' pmega_pseudofinder.txt | grep PMEGA | sort > temp
mv temp pmega_pseudofinder.txt

sed 's/,/\n/g' pphem_pseudofinder.txt | grep PPHEM | sort > temp
mv temp pphem_pseudofinder.txt

sed 's/,/\n/g' pphex_pseudofinder.txt | grep PPHEX | sort > temp
mv temp pphex_pseudofinder.txt

sed 's/,/\n/g' pphym_pseudofinder.txt | grep PPHYM | sort > temp
mv temp pphym_pseudofinder.txt

sed 's/,/\n/g' pphyt_pseudofinder.txt | grep PPHYT | sort > temp
mv temp pphyt_pseudofinder.txt

sed 's/,/\n/g' psart_pseudofinder.txt | grep PSART | sort > temp
mv temp psart_pseudofinder.txt

sed 's/,/\n/g' pspre_pseudofinder.txt | grep PSPRE | sort > temp
mv temp pspre_pseudofinder.txt

sed 's/,/\n/g' ptera_pseudofinder.txt | grep PTERA | sort > temp
mv temp ptera_pseudofinder.txt

sed 's/,/\n/g' ptere_pseudofinder.txt | grep PTERE | sort > temp
mv temp ptere_pseudofinder.txt

sed 's/,/\n/g' pxeno_pseudofinder.txt | grep PXENO | sort > temp
mv temp pxeno_pseudofinder.txt

