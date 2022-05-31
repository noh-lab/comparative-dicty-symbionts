#!/bin/bash

# run rps blast to annotate for COG
# redoing to doublecheck results because of this (when blast filters are applied): 
# https://academic.oup.com/bioinformatics/article/35/15/2697/5239655

rpsblast -query /research/snoh/data/colby/comparative_burk/prokka/baqs159/PROKKA_06192020.faa -db /research/snoh/data/colby/rpsblast/Cog -evalue 0.01 -max_target_seqs 25 -outfmt "6 qseqid sseqid evalue qcovs stitle" > rpsblast_tables/baqs159.prokka.redo.cog

rpsblast -query /research/snoh/data/colby/comparative_burk/prokka/bbqs859/PROKKA_06192020.faa -db /research/snoh/data/colby/rpsblast/Cog -evalue 0.01 -max_target_seqs 25 -outfmt "6 qseqid sseqid evalue qcovs stitle" > rpsblast_tables/bbqs859.prokka.redo.cog

rpsblast -query /research/snoh/data/colby/comparative_burk/prokka/bhqs11/PROKKA_06192020.faa -db /research/snoh/data/colby/rpsblast/Cog -evalue 0.01 -max_target_seqs 25 -outfmt "6 qseqid sseqid evalue qcovs stitle" > rpsblast_tables/bhqs11.prokka.redo.cog

# running controls

rpsblast -query /research/snoh/data/colby/comparative_burk/prokka/pfung/PROKKA_06192020.faa -db /research/snoh/data/colby/rpsblast/Cog -evalue 0.01 -max_target_seqs 25 -outfmt "6 qseqid sseqid evalue qcovs stitle" > rpsblast_tables/pfung.prokka.redo.cog

rpsblast -query /research/snoh/data/colby/comparative_burk/prokka/pmega/PROKKA_06192020.faa -db /research/snoh/data/colby/rpsblast/Cog -evalue 0.01 -max_target_seqs 25 -outfmt "6 qseqid sseqid evalue qcovs stitle" > rpsblast_tables/pmega.prokka.redo.cog

rpsblast -query /research/snoh/data/colby/comparative_burk/prokka/ptere/PROKKA_06192020.faa -db /research/snoh/data/colby/rpsblast/Cog -evalue 0.01 -max_target_seqs 25 -outfmt "6 qseqid sseqid evalue qcovs stitle" > rpsblast_tables/ptere.prokka.redo.cog

rpsblast -query /research/snoh/data/colby/comparative_burk/prokka/pxeno/PROKKA_06192020.faa -db /research/snoh/data/colby/rpsblast/Cog -evalue 0.01 -max_target_seqs 25 -outfmt "6 qseqid sseqid evalue qcovs stitle" > rpsblast_tables/pxeno.prokka.redo.cog


rpsblast -query /research/snoh/data/colby/comparative_burk/prokka/pphex/PROKKA_06192020.faa -db /research/snoh/data/colby/rpsblast/Cog -evalue 0.01 -max_target_seqs 25 -outfmt "6 qseqid sseqid evalue qcovs stitle" > rpsblast_tables/pphex.prokka.redo.cog

rpsblast -query /research/snoh/data/colby/comparative_burk/prokka/pphym/PROKKA_06192020.faa -db /research/snoh/data/colby/rpsblast/Cog -evalue 0.01 -max_target_seqs 25 -outfmt "6 qseqid sseqid evalue qcovs stitle" > rpsblast_tables/pphym.prokka.redo.cog

rpsblast -query /research/snoh/data/colby/comparative_burk/prokka/pphyt/PROKKA_06192020.faa -db /research/snoh/data/colby/rpsblast/Cog -evalue 0.01 -max_target_seqs 25 -outfmt "6 qseqid sseqid evalue qcovs stitle" > rpsblast_tables/pphyt.prokka.redo.cog

rpsblast -query /research/snoh/data/colby/comparative_burk/prokka/pspre/PROKKA_06192020.faa -db /research/snoh/data/colby/rpsblast/Cog -evalue 0.01 -max_target_seqs 25 -outfmt "6 qseqid sseqid evalue qcovs stitle" > rpsblast_tables/pspre.prokka.redo.cog


rpsblast -query /research/snoh/data/colby/comparative_burk/prokka/pcale/PROKKA_06192020.faa -db /research/snoh/data/colby/rpsblast/Cog -evalue 0.01 -max_target_seqs 25 -outfmt "6 qseqid sseqid evalue qcovs stitle" > rpsblast_tables/pcale.prokka.redo.cog

rpsblast -query /research/snoh/data/colby/comparative_burk/prokka/pphem/PROKKA_06192020.faa -db /research/snoh/data/colby/rpsblast/Cog -evalue 0.01 -max_target_seqs 25 -outfmt "6 qseqid sseqid evalue qcovs stitle" > rpsblast_tables/pphem.prokka.redo.cog

rpsblast -query /research/snoh/data/colby/comparative_burk/prokka/psart/PROKKA_06192020.faa -db /research/snoh/data/colby/rpsblast/Cog -evalue 0.01 -max_target_seqs 25 -outfmt "6 qseqid sseqid evalue qcovs stitle" > rpsblast_tables/psart.prokka.redo.cog

rpsblast -query /research/snoh/data/colby/comparative_burk/prokka/ptera/PROKKA_06192020.faa -db /research/snoh/data/colby/rpsblast/Cog -evalue 0.01 -max_target_seqs 25 -outfmt "6 qseqid sseqid evalue qcovs stitle" > rpsblast_tables/ptera.prokka.redo.cog

