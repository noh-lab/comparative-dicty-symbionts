#!/bin/bash
# parse rpsblast results
# keep top hit and filter by qcovs>=70, get COG number from stitle

awk -F'[\t,]' '!x[$1]++ && $4>=70 {print $1,$5}' OFS="\t" rpsblast_tables/baqs159.prokka.redo.cog > baqs159.cog.redo.table
awk -F'[\t,]' '!x[$1]++ && $4>=70 {print $1,$5}' OFS="\t" rpsblast_tables/bbqs859.prokka.redo.cog > bbqs859.cog.redo.table
awk -F'[\t,]' '!x[$1]++ && $4>=70 {print $1,$5}' OFS="\t" rpsblast_tables/bhqs11.prokka.redo.cog > bhqs11.cog.redo.table

awk -F'[\t,]' '!x[$1]++ && $4>=70 {print $1,$5}' OFS="\t" rpsblast_tables/pfung.prokka.redo.cog > pfung.cog.redo.table
awk -F'[\t,]' '!x[$1]++ && $4>=70 {print $1,$5}' OFS="\t" rpsblast_tables/pmega.prokka.redo.cog > pmega.cog.redo.table
awk -F'[\t,]' '!x[$1]++ && $4>=70 {print $1,$5}' OFS="\t" rpsblast_tables/ptere.prokka.redo.cog > ptere.cog.redo.table
awk -F'[\t,]' '!x[$1]++ && $4>=70 {print $1,$5}' OFS="\t" rpsblast_tables/pxeno.prokka.redo.cog > pxeno.cog.redo.table

awk -F'[\t,]' '!x[$1]++ && $4>=70 {print $1,$5}' OFS="\t" rpsblast_tables/pphex.prokka.redo.cog > pphex.cog.redo.table
awk -F'[\t,]' '!x[$1]++ && $4>=70 {print $1,$5}' OFS="\t" rpsblast_tables/pphym.prokka.redo.cog > pphym.cog.redo.table
awk -F'[\t,]' '!x[$1]++ && $4>=70 {print $1,$5}' OFS="\t" rpsblast_tables/pphyt.prokka.redo.cog > pphyt.cog.redo.table
awk -F'[\t,]' '!x[$1]++ && $4>=70 {print $1,$5}' OFS="\t" rpsblast_tables/pspre.prokka.redo.cog > pspre.cog.redo.table

awk -F'[\t,]' '!x[$1]++ && $4>=70 {print $1,$5}' OFS="\t" rpsblast_tables/pcale.prokka.redo.cog > pcale.cog.redo.table
awk -F'[\t,]' '!x[$1]++ && $4>=70 {print $1,$5}' OFS="\t" rpsblast_tables/pphem.prokka.redo.cog > pphem.cog.redo.table
awk -F'[\t,]' '!x[$1]++ && $4>=70 {print $1,$5}' OFS="\t" rpsblast_tables/psart.prokka.redo.cog > psart.cog.redo.table
awk -F'[\t,]' '!x[$1]++ && $4>=70 {print $1,$5}' OFS="\t" rpsblast_tables/ptera.prokka.redo.cog > ptera.cog.redo.table

