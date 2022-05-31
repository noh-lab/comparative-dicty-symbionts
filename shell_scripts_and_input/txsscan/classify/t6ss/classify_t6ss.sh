#!/bin/bash

cat tssB.faa tssB_secret6.faa tssB_manual.faa > tssB_all.faa
muscle -in tssB_all.faa -out tssB.muscle.faa
python ../../stable.py tssB_all.faa tssB.muscle.faa > tssB_all.stable.faa

awk '/^>/{print ">" ++i; next}{print}' tssB_all.stable.faa > tssB.temp
awk 'FNR==NR{a[">"$1]=$2;next} $1 in a{sub(/[0-9]+/,a[$1])}1' t6ss_sub.txt tssB.temp > tssB_all.namefix.faa
 
FastTree -lg < tssB_all.namefix.faa > tssB_all.lg.tre

cat tssC.faa tssC_secret6.faa tssC_manual.faa > tssC_all.faa
muscle -in tssC_all.faa -out tssC.muscle.faa
python ../../stable.py tssC_all.faa tssC.muscle.faa > tssC_all.stable.faa

awk '/^>/{print ">" ++i; next}{print}' tssC_all.stable.faa > tssC.temp
awk 'FNR==NR{a[">"$1]=$2;next} $1 in a{sub(/[0-9]+/,a[$1])}1' t6ss_sub.txt tssC.temp > tssC_all.namefix.faa
 
FastTree -lg < tssC_all.namefix.faa > tssC_all.lg.tre

cat tssF.faa tssF_secret6.faa tssF_manual.faa > tssF_all.faa
muscle -in tssF_all.faa -out tssF.muscle.faa
python ../../stable.py tssF_all.faa tssF.muscle.faa > tssF_all.stable.faa

awk '/^>/{print ">" ++i; next}{print}' tssF_all.stable.faa > tssF.temp
awk 'FNR==NR{a[">"$1]=$2;next} $1 in a{sub(/[0-9]+/,a[$1])}1' t6ss_sub_2.txt tssF.temp > tssF_all.namefix.faa
 
FastTree -lg < tssF_all.namefix.faa > tssF_all.lg.tre


cat tss*lg.tre > t6ss.lg.infile

ASTRID -i t6ss.lg.infile -o t6ss.lg.astrid

