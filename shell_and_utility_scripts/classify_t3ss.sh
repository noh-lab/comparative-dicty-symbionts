#!/bin/bash

cat sctJ.faa sctJ_from_t3enc.faa sctJ_manual.faa > sctJ_all.faa
muscle -in sctJ_all.faa -out sctJ_t3enc.muscle.faa
cat sctN.faa sctN_from_t3enc.faa sctN_manual.faa > sctN_all.faa
muscle -in sctN_all.faa -out sctN_t3enc.muscle.faa
cat sctV.faa sctV_from_t3enc.faa sctV_manual.faa > sctV_all.faa
muscle -in sctV_all.faa -out sctV_t3enc.muscle.faa


# fix order to input order and rename (to ultimately generate species tree)

python ../../stable.py sctJ_all.faa sctJ_t3enc.muscle.faa > sctJ_all.stable.faa
python ../../stable.py sctN_all.faa sctN_t3enc.muscle.faa > sctN_all.stable.faa
python ../../stable.py sctV_all.faa sctV_t3enc.muscle.faa > sctV_all.stable.faa


# substitute names of individual genes to facilitate making species tree later

awk '/^>/{print ">" ++i; next}{print}' sctJ_all.stable.faa > sctJ.temp
awk 'FNR==NR{a[">"$1]=$2;next} $1 in a{sub(/[0-9]+/,a[$1])}1' t3ss_sub.txt sctJ.temp > sctJ_all.namefix.faa

awk '/^>/{print ">" ++i; next}{print}' sctN_all.stable.faa > sctN.temp
awk 'FNR==NR{a[">"$1]=$2;next} $1 in a{sub(/[0-9]+/,a[$1])}1' t3ss_sub.txt sctN.temp > sctN_all.namefix.faa

awk '/^>/{print ">" ++i; next}{print}' sctV_all.stable.faa > sctV.temp
awk 'FNR==NR{a[">"$1]=$2;next} $1 in a{sub(/[0-9]+/,a[$1])}1' t3ss_sub.txt sctV.temp > sctV_all.namefix.faa


# make gene trees, using Le and Gascuel a.a. Sub model here

FastTree -lg < sctJ_all.namefix.faa > sctJ_all.lg.tre
FastTree -lg < sctN_all.namefix.faa > sctN_all.lg.tre
FastTree -lg < sctV_all.namefix.faa > sctV_all.lg.tre


# concatenate these gene trees and run ASTRID for species tree

cat sct*lg.tre > t3ss.lg.infile

ASTRID -i t3ss.lg.infile -o t3ss.lg.astrid

Astral -i t3ss.lg.infile -o t3ss.lg.astral
