# Comparative genomics of *Dictyostelium discoideum*-symbiont *Paraburkholderia*
Repository to support: Comparative genomics support reduced-genome *Paraburkholderia* symbionts of *Dictyostelium discoideum* amoebas are ancestrally adapted professional symbionts

Supporting data and scripts are organized as follow:
- shell_scripts_and_input: scripts and input data to recreate analyses prior to visualization
- R_scripts_and_input: scripts and input data to recreate visualizations and statistical analyses
- additional_output: data that were not directly used for visualizations or statistics; includes data used in the genome browser @ burk.colby.edu

These data were generated from genome FASTA files available via NCBI. A descriptive version of all analyses and steps can be found below.
&nbsp;

## a. *Paraburkholderia* genome selection and gene prediction
Select the genomes you want to compare. We started with the three *D. discoideum* (dicty)-associated *Paraburkholderia* genomes and selected 12 additional *Paraburkholderia* genomes for comparative analysis with the three *D. discoideum*-associated species genomes. 

| Species | Abbreviation | RefSeq file |
| -------- | -------- | -------- |
| ***dicty-associated*** |
| *P. agricolaris* | PAGRI | GCF_009455635.1_ASM945563v1_genomic.fna |
| *P. bonniea* | PBONN | GCF_009455625.1_ASM945562v1_genomic.fna |
| *P. hayleyella* |PHAYL | GCF_009455685.1_ASM945568v1_genomic.fna |
| ***core controls*** | 
| *P. fungorum* | PFUNG | GCF_000961515.1_ASM96151v1_genomic.fna |
| *P. sprentiae*| PSPRE | GCF_001865575.1_ASM186557v1_genomic.fna |
| *P. terrae* | PTERE | GCF_002902925.1_ASM290292v1_genomic.fna |
| *P. xenovorans* | PXENO | GCF_000756045.1_ASM75604v1_genomic.fna |
| ***additional plant-associated*** |
| *P. phenoliruptrix* | PPHEX | GCF_000300095.1_ASM30009v1_genomic.fna |
| *P. phymatum* | PPHYM | GCF_000020045.1_ASM2004v1_genomic.fna |
| *P. phytofirmans* | PPHYT | GCF_000020125.1_ASM2012v1_genomic.fna |
| *P. megapolitana* | PMEGA | GCF_007556815.1_ASM755681v1_genomic.fna |
| ***additional free-living*** |
| *P. caledonica* | PCALE | GCF_003330745.1_ASM333074v1_genomic.fna |
| *P. phenazinium* | PPHEM | GCF_900100735.1_IMG2651870170_genomic.fna |
| *P. sartisoli* | PSART | GCF_900107685.1_IMG2651870102_genomic.fna |
| *P. terricola* | PTERA | GCF_003330825.1_ASM333082v1_genomic.fna |

&nbsp;

Re-predict genes in all genomes with [Prokka](https://github.com/tseemann/prokka), e.g.
```
prokka --force --outdir baqs159 --locustag PAGRI --proteins Burkholderia_pseudomallei_K96243_132.gbk --gcode 11 --genus Paraburkholderia --species agricolaris --strain BaQS159 --gram neg --addgenes --rnammer GCF_009455635.1_ASM945563v1_genomic.fna
```

We used custom scripts to execute the command(s) above for all genomes:

`run_prokka_dictyburk.sh`
`run_prokka_control.sh`

&nbsp;

### Pseudogene detections
Use [Pseudofinder](https://github.com/filip-husnik/pseudofinder) to detect likely pseudogenes from each Prokka genbank file, e.g.
```
pseudofinder.py annotate --diamond --skip_makedb -g /research/snoh/data/colby/comparative_burk/prokka/baqs159/PROKKA_06192020.gbk -db /research/snoh/diamonddb/refseq_protein_nr.dmnd -op baqs159
```

Make a list of predicted genes that were likely pseudogenes from Pseudofinder output files, e.g.
```
# second column (prokka gene id) of each output file contains detected pseudogenes
grep PAGRI baqs159_pseudos.gff | sort -k1,1 -k4,4n | awk -F "\t" '{split($9,a,";"); split(a[3],b,"="); split(a[1],c,":"); print b[2]}' OFS="\t" > pagri_pseudofinder.txt

# make list with one gene id per line for future exclusion
sed 's/,/\n/g' pagri_pseudofinder.txt | grep PAGRI | sort > temp
mv temp pagri_pseudofinder.txt
```

We used custom scripts to execute the command(s) above for all genomes:

`run_pseudofinder_dictyburk.sh`
`run_pseudofinder_control.sh`

&nbsp;

### Metadata collection
Manually gather genome metadata from Prokka report files and Pseudofinder summaries. 

Manually calculate GC% directly from NCBI fna files: 
```
for FILE in GCF*.fna; do
	echo "$FILE"
	var1=$(grep -v ">" `echo "$FILE"` | tr -d -c GCgc | wc -c)
	var2=$(grep -v ">" `echo "$FILE"` | tr -d -c ATGCatgc | wc -c)
	echo "$var1, $var2"
	awk -v v1=$var1 -v v2=$var2 'BEGIN {print v1/v2}'
done
```

&nbsp;

## b. Whole genome alignment
Use [progressiveMauve](https://darlinglab.org/mauve/user-guide/progressivemauve.html) to align all finished genomes at once. 

```
progressiveMauve --output=10cont.xmfa --seed-family ../prokka/baqs159/PROKKA_06192020.fna ../prokka/bbqs859/PROKKA_06192020.fna ../prokka/bhqs11/PROKKA_06192020.fna ../prokka/pcale/PROKKA_06192020.fna ../prokka/pfung/PROKKA_06192020.fna ../prokka/pmega/PROKKA_06192020.fna ../prokka/pphex/PROKKA_06192020.fna ../prokka/pphym/PROKKA_06192020.fna ../prokka/pphyt/PROKKA_06192020.fna ../prokka/pspre/PROKKA_06192020.fna ../prokka/ptera/PROKKA_06192020.fna ../prokka/ptere/PROKKA_06192020.fna ../prokka/pxeno/PROKKA_06192020.fna &
```

Next, extract the LCBs. We did this for the 3 dicty-burk genomes and 4 core controls (the numbers are in the order you listed the genomes during the initial alignment):
```
./projectAndStrip 10cont.xmfa 10cont.4core.strip 0 1 2 4 9 11 12

./makeBadgerMatrix 10cont.4core.strip 10cont.4core.perm 10cont.4core.lcb
```

The `*.lcb` file contains the positions of any locally colinear blocks (LCB) across the genomes you selected. The `*.perm` file contains the relative orientations (forward vs. reverse) of LCB across your genomes. 

Import into `R` for visualization. We used a custom script:
`mauve.visualize_clean.R`

&nbsp;

## c. Core genome molecular evolution
### Define core genome
First remove pseudogenes from each prokka gff file, e.g.
```
awk '{print $1}' baqs159.prokka.ffn.fai > pagri_genes.txt

awk 'NR==FNR {a[$1]; next} !($1 in a) {print}' ../pseudofinder/pagri_pseudofinder.txt pagri_genes.txt > pagri_intact.txt

grep -f pagri_intact.txt baqs159.prokka.gff > baqs159.prokka.intact.gff

# add headers; ##FASTA; and fastas, e.g.
grep '##' baqs159.prokka.gff | sed '$d' > head

echo '##FASTA' >> baqs159.prokka.gff

cat head baqs159.prokka.intact.gff ../prokka/baqs159/PROKKA_06192020.fsa > temp

mv temp baqs159.prokka.intact.gff
```

We used custom scripts to execute the command(s) above for all genomes:

`prep_roary_1_pseudoremov.sh`
`prep_roary_2_fix_gff.sh`
&nbsp;

Run Roary:
```
singularity run -i -e -H /export/groups/snoh/snoh/12cont_intact/ /usr/local/singularity/roary_latest.sif roary -r -s -p 8 -i 70 *.gff
```

Find core genes:
```
singularity run -i -e -H /export/groups/snoh/snoh/12cont_intact/ /usr/local/singularity/roary_latest.sif query_pan_genome -a intersection *gff

mkdir paraburk_core

sed 's/:\s/\t/g' pan_genome_results > paraburk_core/core_genes.table
cut -f 1 paraburk_core/core_genes.table > paraburk_core/core_genes.list
```

Get faa and fna multifastas for each core gene and create alignments. We used custom scripts to execute this step:

`prep_codeml_1_core_multifastas.sh`
`prep_codeml_2_muscle_to_pal2nal.sh`
&nbsp;

Split lists of genes for pseudo-parallel execution:
```
split -l 210 paraburk_core/core_genes.list paraburk_core/core_genes.
```

&nbsp;

### Apply different models of molecular evolution
Run [CODEML](http://abacus.gene.ucl.ac.uk/software/paml.html). We used these custom scripts:

`run_codeml_a*.sh`
`summarize_codeml_1_logl.sh`
&nbsp;

These scripts execute the following steps:
- [ ] Clean up the pal2nal alignment
- [ ] Run base model with no variation in omega across tree
- [ ] Run "symbiotic" model 
- [ ] Run "dicty" model
- [ ] Run "reduced" model
- [ ] Get loglikelihood scores for model selection


Compile results for each gene:
```
cat model_summaries/*summary.txt > ../12cont_intact/codeml_likelihood.models.txt
```

In `R`, find best model for each gene and export lists of genes for each model using `codeml_parse_visualize_clean.R`.
&nbsp;

Get model parameter estimates (kappa, w0, w1, logL) for the model that best fits each gene. We used a custom script:
`summarize_codeml_2_parameters.sh`
&nbsp;

Gather these parameter estimates for further analysis:
```
cat model_summaries/*estimates.txt > ../12cont_intact/codeml_omega.diff1.txt
```

In `R`, continue analysis with `codeml.parse_visualize_clean.R`

&nbsp;

## d. Gene functional annotation

### COG (Clusters of Orthologous Groups)
You need to start with a COG position-specific scoring matrix file: `cdd.tar.gz` downloaded from ncbi here: ftp://ftp.ncbi.nih.gov/pub/mmdb/cdd/
&nbsp;

You also need the table that decodes each COG number `cognames2003-2014.tab` and the table that names each COG functional category `fun2003-2014.tab` from here:
ftp://ftp.ncbi.nih.gov/pub/COG/COG2014/data
&nbsp;

Make a database of these profiles,
```
makeprofiledb -title COG.v.20200430 -in Cog.pn -out Cog -threshold 9.82 -scale 100.0 -dbtype rps -index true
```

Run `rpsblast`, e.g.
```
rpsblast -query /research/snoh/data/colby/comparative_burk/prokka/baqs159/PROKKA_06192020.faa -db /research/snoh/data/colby/rpsblast/Cog -evalue 0.01 -max_target_seqs 25 -outfmt "6 qseqid sseqid evalue qcovs stitle" > rpsblast_tables/baqs159.prokka.cog
```

Keep the top `rpsblast` hit and filter by qcovs>=70, and get COG number from stitle, e.g.
```
awk -F'[\t,]' '!x[$1]++ && $4>=70 {print $1,$5}' OFS="\t" rpsblast_tables/baqs159.prokka.cog > baqs159.cog.table
```

For cogs with multiple functional categories, keep first one (it’s the primary one), e.g.
```
for FILE in *table; do
	KEEP="${FILE%.c*}"
	awk -F "\t" 'NR==FNR {a[$1]=$2;next} {if ($2 in a){print $1, $2, a[$2]} else {print $0}}' OFS="\t" cognames2003-2014.tab "$FILE" > temp
	awk -F "\t" '{if ( length($3)>1 ) { $3 = substr($3, 0, 1) } else { $3 = $3 }; print}' OFS="\t" temp > temp2
	awk -F "\t" 'NR==FNR {a[$1]=$2;next} {if ($3 in a){print $0, a[$3]} else {print $0}}' OFS="\t" fun2003-2014.tab temp2 > "$KEEP".cog.categorized
done

rm temp temp2
```

Remove predicted pseudogenes, e.g.
```
awk 'NR==FNR {a[$1]; next} !($1 in a){print $0}' ../pseudofinder/pagri_pseudofinder.txt baqs159.cog.categorized > baqs159.cog.intact.categorized
```

We used custom scripts to execute the command(s) above for all genomes:

`cog_1_run_rpsblast.sh`
`cog_2_parse_rpsblast.sh`
`cog_3_categorize_cogs.sh`
&nbsp;

Output files were visualized in R with this script: `roary_cog.visualize_clean.R`

&nbsp;

### KO (KEGG Orthology)
Run [BlastKOALA](https://www.kegg.jp/blastkoala/) for each genome as it’s more specific than ghostKOALA; run with prokaryote genus and eukaryote family as target database.
&nbsp;

Save results to a text file, e.g. `baqs159.blastkoala.txt`

Remove predicted pseudogenes, e.g.
```
awk 'NR==FNR {a[$1]; next} !($1 in a){print $0}' ../pseudofinder/pagri_pseudofinder.txt baqs159.blastkoala.txt > baqs159.blastkoala.intact.txt
```

Organize these files, e.g.
```
# remove empty lines
awk '$2!="" {print $0}' baqs159.blastkoala.intact.txt > baqs159.blastkoala.intact.keggmapper.txt
```

Load these files into [KEGGmapper - Reconstruct](https://www.genome.jp/kegg/mapper/reconstruct.html)
From there, you can check modules and pathways for completeness. 

&nbsp;

## e. Essential amino acid biosynthetic repertoire

Use [GapMind](https://papers.genomics.lbl.gov/cgi-bin/gapView.cgi) to detect amino acid biosynthesis pathway genes.

Pull multi fasta of intact proteins in each genome, e.g.
```
./FastaToTbl ../roary/baqs159.prokka.v3.faa | grep -vwf ../pseudofinder/pagri_pseudofinder.txt | ./TblToFasta > baqs159.prokka.intact.faa
```

Download candidates and steps from results.

&nbsp;

## f. Protein secretion system repertoire and effector prediction

### Secretion systems
Use [TXSScan via MacSyFinder](https://github.com/gem-pasteur/macsyfinder) to detect secretion system components in your genome. We used MacSyFinder on [Institute Pasteur's Galaxy instance](https://galaxy.pasteur.fr/).

There are two output files, save both:
* output (saved as `*txsscan_output.txt`) - log of genome analysis, where information on incomplete systems can be found
* report_output (saved as `*txsscan_report.txt`) - table of predicted secretion system components for complete systems 

Reformat `*txsscan_report.txt` files and reorganize secretion systems into the same order, e.g.
```
# simplify information
cut -f 1,5,7-10 baqs159.txsscan_report.txt | sed 's/UserReplicon_//g' | sed 's/T4SS_typeG/cT4SS_typeG/g' | sed 's/T4SS_typeT/cT4SS_typeT/g' | awk '{print $3,$4,$2,$1}' OFS="\t" > temp1

# reorder secretion systems
while read c2; do grep -w $c2 temp1; done < reorder_ss_components.txt > baqs159.ss_components.txt

# examine and remove any pseudogenes manually
awk 'NR==FNR {a[$1]; next} ($4 in a){print}' ../pseudofinder/pagri_pseudofinder.txt baqs159.ss_components.txt
```

Manually scan the output to search for incomplete secretion systems in the genomes (not reported as complete by txsscan), that were incorrectly identified (e.g. TssC as IglB) or incorrectly disambiguated when adjacent to another secretion system. 
&nbsp;

To do this, we first inspected software behavior using two well studied genomes, downloaded from the [Burkholderia Genome Database](burkholderia.com):
* Burkholderia mallei ATCC 23344 
* Burkholderia pseudomallei K96243 

We identified recurring errors (above) that did not agree with known secretion systems from these genomes and manually added corrected these errors. 

Compile manually corrected secretion systems into a new file: `all.manual_ss_components.txt`

Classify type 3 and type 6 secretion systems. We downloaded amino acid sequences from [T3Enc](http://www.szu-bioinf.org/T3Enc/index.html) and [SecReT6](https://bioinfo-mml.sjtu.edu.cn/SecReT6/index.php) databases and compiled them into files, e.g. `sctJ_from_t3enc.faa`.
* T3SS - sctJ, sctN, sctV
* T6SS - tssB, tssC (iglB), tssF

Get amino acid sequences of component genes for each genome, e.g.
```
for gene in `echo "sctJ"`; do
    agri=$(grep "$gene" all.manual_ss_components.txt | grep AGRI | cut -f 4)
    samtools faidx ../roary/baqs159.prokka.v3.faa `echo "$agri"` > "$gene"_manual.faa
done
```

Next, classify, e.g.
```
# align
cat sctJ.faa sctJ_from_t3enc.faa sctJ_manual.faa > sctJ_all.faa
muscle -in sctJ_all.faa -out sctJ_t3enc.muscle.faa

# fix order
python ../../stable.py sctJ_all.faa sctJ_t3enc.muscle.faa > sctJ_all.stable.faa

# substitute names (species rather than gene)
awk '/^>/{print ">" ++i; next}{print}' sctJ_all.stable.faa > sctJ.temp
awk 'FNR==NR{a[">"$1]=$2;next} $1 in a{sub(/[0-9]+/,a[$1])}1' t3ss_sub.txt sctJ.temp > sctJ_all.namefix.faa

# estimate gene trees
FastTree -lg < sctJ_all.namefix.faa > sctJ_all.lg.tre

# estimate species tree
cat sct*lg.tre > t3ss.lg.infile
ASTRID -i t3ss.lg.infile -o t3ss.lg.astrid
Astral -i t3ss.lg.infile -o t3ss.lg.astral
```

Visualize resulting trees in [iToL](itol.embl.de).

We used custom scripts to execute the command(s) above for all genomes:

`secretion_1_get_gene_fastas.sh`
`secretion_2_get_manual_fastas.sh`
`classify_t3ss.sh`
`classify_t6ss.sh`

We visually summarized detected secretion systems in `R` with custom script:

`txsscan.visualize_clean.R`

&nbsp;

### Secreted effectors
We used the following webservers:
* [BastionHub](https://bastionhub.erc.monash.edu/)
* [EffectiveELD](https://effectors.csb.univie.ac.at/)
* [VFDB](http://www.mgc.ac.cn/VFs/)

For the first two, submit predicted protein sequences and save all significant results. We used both HMMER and BastionX predictions from BastionHub. We used a lower Minimal score threshold (2) for EffectiveELD. 
&nbsp;
For EffectiveELD results, filter results for pfams of known eukaryote-like domains expected to be in virulence from the literature:
```
awk -F "\t" 'NR==FNR {a[$1]=$0; next} ($2 in a){print $1,$4,a[$2]}' OFS="\t" eukaryote_pfams.sorted.txt baqs159.ELD > baqs159.ELD.overlap
```

For VFDB, based on classification results above, we downloaded known [T3SS effectors from *Bordatella*](http://www.mgc.ac.cn/cgi-bin/VFs/genus.cgi?Genus=Bordetella), and [T6SS effectors from *Burkholderia*](http://www.mgc.ac.cn/cgi-bin/VFs/genus.cgi?Genus=Burkholderia). 

Save these as a file: `known_t36se_vfdb.fa`

Search against each dicty-burk genomes, e.g.
```
# make diamond db of each genome
diamond makedb --in ../roary/baqs159.prokka.v3.faa -d pagri

# run a search in blastp mode
diamond blastp -d pagri -q known_t36se_vfdb.fa -o pagri_t36se_matches.tsv -b8 -c1 -p 8
```

&nbsp;

## g. Horizontally transferred genetic element detection

### Genomic Islands 
Use [IslandViewer](https://www.pathogenomics.sfu.ca/islandviewer/) to detect genomic islands. You will need to upload each chromosome separately. 

&nbsp;

### IS Elements
Use [ISFinder](https://isfinder.biotoul.fr/) to detect IS elements. 

Use `Tools - BLASTN` and process the results as you would a typical blast query. 
We used the following steps:

- [ ] Query database with e-value threshold 1e-10
- [ ] Sort results by coordinate and detect hits that overlap
- [ ] Within each overlapping set of hits, mark best evalue and best bitscore hit
- [ ] Filter results by 70% alignment length
- [ ] Use `Tools - Search` to determine IS family of best hit

&nbsp;

### Horizontally tranferred genes
Use [Darkhorse2](https://github.com/spodell/Darkhorse2) to detect individual horizontally transferred genes. Darkhorse2 must be locally installed with appropriate databases. 

Create diamond database of taxonomically informative sequences `databasename_informative.dmnd`:
```
/usr/local/src/Darkhorse2-DarkHorse-2.0_rev09/install_darkhorse2.pl -c config_template -i . -u -n 12
```

Run diamond BLASTP, e.g.
```
diamond blastp -d /export/groups/snoh/snoh/diamonddb/snohdarkhorse_20210104_informative.dmnd -q baqs159.prokka.faa -a baqs159.prokka.daa -e1e-6 -t . -k500 -p12 -b8 -c1

diamond view -a baqs159.prokka.daa -f tab  -o baqs159.prokka.dh2.m8 -k500
```

Create exclude lists, e.g.
```
generate_dh_self_keywords.pl -i baqs159_taxid -c config_template #2428
```

Modified the `*_exclude_list_strain` file to contain taxids of self and sister species (e.g. *P. agricolaris* + *P. fungorum*).

Run Darkhorse2, e.g.
```
/usr/local/src/Darkhorse2-DarkHorse-2.0_rev09/bin/darkhorse2.pl -c config_template -t baqs159.prokka.dh2.m8 -e dh_keywords_2428/Paraburkholderia_agricolaris_exclude_lists/Paraburkholderia_agricolaris_exclude_list_strain -g baqs159.prokka.faa -f 0.02
```

Keep putative results with reasonable LPI scores, e.g.
```
awk -F "\t" '$6>0.4 &&$6<0.85' darkhorse/pagri_filt_0.02/*smry | sort -k1,1 | awk -F "\t" '{print $1,$4,$6,$7,$10,$11,$14,$16,$17}' OFS="\t > pagri_f0.02.darkhorse2_filtered.txt 
```

---

