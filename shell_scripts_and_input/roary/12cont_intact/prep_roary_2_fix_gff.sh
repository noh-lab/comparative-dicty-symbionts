#!/bin/bash

grep '##' baqs159.prokka.v3.gff | sed '$d' > head
echo '##FASTA' >> 12cont_intact/baqs159.prokka.intact.gff
cat head 12cont_intact/baqs159.prokka.intact.gff ../prokka/baqs159/PROKKA_06192020.fsa > temp
mv temp 12cont_intact/baqs159.prokka.intact.gff

grep '##' bbqs859.prokka.v3.gff | sed '$d' > head
echo '##FASTA' >> 12cont_intact/bbqs859.prokka.intact.gff
cat head 12cont_intact/bbqs859.prokka.intact.gff ../prokka/bbqs859/PROKKA_06192020.fsa > temp
mv temp 12cont_intact/bbqs859.prokka.intact.gff

grep '##' bhqs11.prokka.v3.gff | sed '$d' > head
echo '##FASTA' >> 12cont_intact/bhqs11.prokka.intact.gff
cat head 12cont_intact/bhqs11.prokka.intact.gff ../prokka/bhqs11/PROKKA_06192020.fsa > temp
mv temp 12cont_intact/bhqs11.prokka.intact.gff


grep '##' pcale.prokka.v3.gff | sed '$d' > head
echo '##FASTA' >> 12cont_intact/pcale.prokka.intact.gff
cat head 12cont_intact/pcale.prokka.intact.gff ../prokka/pcale/PROKKA_06192020.fsa > temp
mv temp 12cont_intact/pcale.prokka.intact.gff

grep '##' pfung.prokka.v3.gff | sed '$d' > head
echo '##FASTA' >> 12cont_intact/pfung.prokka.intact.gff
cat head 12cont_intact/pfung.prokka.intact.gff ../prokka/pfung/PROKKA_06192020.fsa > temp
mv temp 12cont_intact/pfung.prokka.intact.gff

grep '##' pmega.prokka.v3.gff | sed '$d' > head
echo '##FASTA' >> 12cont_intact/pmega.prokka.intact.gff
cat head 12cont_intact/pmega.prokka.intact.gff ../prokka/pmega/PROKKA_06192020.fsa > temp
mv temp 12cont_intact/pmega.prokka.intact.gff

grep '##' pphem.prokka.v3.gff | sed '$d' > head
echo '##FASTA' >> 12cont_intact/pphem.prokka.intact.gff
cat head 12cont_intact/pphem.prokka.intact.gff ../prokka/pphem/PROKKA_06192020.fsa > temp
mv temp 12cont_intact/pphem.prokka.intact.gff

grep '##' pphex.prokka.v3.gff | sed '$d' > head
echo '##FASTA' >> 12cont_intact/pphex.prokka.intact.gff
cat head 12cont_intact/pphex.prokka.intact.gff ../prokka/pphex/PROKKA_06192020.fsa > temp
mv temp 12cont_intact/pphex.prokka.intact.gff

grep '##' pphym.prokka.v3.gff | sed '$d' > head
echo '##FASTA' >> 12cont_intact/pphym.prokka.intact.gff
cat head 12cont_intact/pphym.prokka.intact.gff ../prokka/pphym/PROKKA_06192020.fsa > temp
mv temp 12cont_intact/pphym.prokka.intact.gff

grep '##' pphyt.prokka.v3.gff | sed '$d' > head
echo '##FASTA' >> 12cont_intact/pphyt.prokka.intact.gff
cat head 12cont_intact/pphyt.prokka.intact.gff ../prokka/pphyt/PROKKA_06192020.fsa > temp
mv temp 12cont_intact/pphyt.prokka.intact.gff

grep '##' psart.prokka.v3.gff | sed '$d' > head
echo '##FASTA' >> 12cont_intact/psart.prokka.intact.gff
cat head 12cont_intact/psart.prokka.intact.gff ../prokka/psart/PROKKA_06192020.fsa > temp
mv temp 12cont_intact/psart.prokka.intact.gff

grep '##' pspre.prokka.v3.gff | sed '$d' > head
echo '##FASTA' >> 12cont_intact/pspre.prokka.intact.gff
cat head 12cont_intact/pspre.prokka.intact.gff ../prokka/pspre/PROKKA_06192020.fsa > temp
mv temp 12cont_intact/pspre.prokka.intact.gff

grep '##' ptera.prokka.v3.gff | sed '$d' > head
echo '##FASTA' >> 12cont_intact/ptera.prokka.intact.gff
cat head 12cont_intact/ptera.prokka.intact.gff ../prokka/ptera/PROKKA_06192020.fsa > temp
mv temp 12cont_intact/ptera.prokka.intact.gff

grep '##' ptere.prokka.v3.gff | sed '$d' > head
echo '##FASTA' >> 12cont_intact/ptere.prokka.intact.gff
cat head 12cont_intact/ptere.prokka.intact.gff ../prokka/ptere/PROKKA_06192020.fsa > temp
mv temp 12cont_intact/ptere.prokka.intact.gff

grep '##' pxeno.prokka.v3.gff | sed '$d' > head
echo '##FASTA' >> 12cont_intact/pxeno.prokka.intact.gff
cat head 12cont_intact/pxeno.prokka.intact.gff ../prokka/pxeno/PROKKA_06192020.fsa > temp
mv temp 12cont_intact/pxeno.prokka.intact.gff

rm head
