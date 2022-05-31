#!/bin/bash

prokka --force --outdir pphym --locustag PPHYM --proteins Burkholderia_pseudomallei_K96243_132.gbk --gcode 11 --genus Paraburkholderia --species phymatum --strain STM815 --gram neg --addgenes --rnammer GCF_000020045.1_ASM2004v1_genomic.fna 
 
prokka --force --outdir pphyt --locustag PPHYT --proteins Burkholderia_pseudomallei_K96243_132.gbk --gcode 11 --genus Paraburkholderia --species phytofirmans --strain PsJN --gram neg --addgenes --rnammer GCF_000020125.1_ASM2012v1_genomic.fna 
 
prokka --force --outdir pphex --locustag PPHEX --proteins Burkholderia_pseudomallei_K96243_132.gbk --gcode 11 --genus Paraburkholderia --species phenoliruptrix --strain BR3459a --gram neg --addgenes --rnammer GCF_000300095.1_ASM30009v1_genomic.fna
 
prokka --force --outdir pxeno --locustag PXENO --proteins Burkholderia_pseudomallei_K96243_132.gbk --gcode 11 --genus Paraburkholderia --species xenovorans --strain LB400 --gram neg --addgenes --rnammer GCF_000756045.1_ASM75604v1_genomic.fna
 
prokka --force --outdir pfung --locustag PFUNG --proteins Burkholderia_pseudomallei_K96243_132.gbk --gcode 11 --genus Paraburkholderia --species fungorum --strain BAA-463 --gram neg --addgenes --rnammer GCF_000961515.1_ASM96151v1_genomic.fna
 
prokka --force --outdir ptere --locustag PTERE --proteins Burkholderia_pseudomallei_K96243_132.gbk --gcode 11 --genus Paraburkholderia --species terrae --strain DSM17804 --gram neg --addgenes --rnammer GCF_002902925.1_ASM290292v1_genomic.fna
 
prokka --force --outdir pcale --locustag PCALE --proteins Burkholderia_pseudomallei_K96243_132.gbk --gcode 11 --genus Paraburkholderia --species caledonica --strain PHRS4 --gram neg --addgenes --rnammer GCF_003330745.1_ASM333074v1_genomic.fna 
 
prokka --force --outdir ptera --locustag PTERA --proteins Burkholderia_pseudomallei_K96243_132.gbk --gcode 11 --genus Paraburkholderia --species terricola --strain mHS1 --gram neg --addgenes --rnammer GCF_003330825.1_ASM333082v1_genomic.fna

prokka --force --outdir pphem --locustag PPHEM --proteins Burkholderia_pseudomallei_K96243_132.gbk --gcode 11 --genus Paraburkholderia --species phenazinium --strain LMG2247 --gram neg --addgenes --rnammer GCF_900100735.1_IMG2651870170_genomic.fna

prokka --force --outdir psart --locustag PSART --proteins Burkholderia_pseudomallei_K96243_132.gbk --gcode 11 --genus Paraburkholderia --species sartisoli --strain LMG24000 --gram neg --addgenes --rnammer GCF_900107685.1_IMG2651870102_genomic.fna

prokka --force --outdir pmega --locustag PMEGA --proteins Burkholderia_pseudomallei_K96243_132.gbk --gcode 11 --genus Paraburkholderia --species megapolitana --strain LMG23650 --gram neg --addgenes --rnammer GCF_007556815.1_ASM755681v1_genomic.fna

prokka --force --outdir pspre --locustag PSPRE --proteins Burkholderia_pseudomallei_K96243_132.gbk --gcode 11 --genus Paraburkholderia --species sprentiae --strain WSM5005 --gram neg --addgenes --rnammer GCF_001865575.1_ASM186557v1_genomic.fna


#cp pphym/PROKKA_*.gff pphym.prokka.v3.gff
#cp pphyt/PROKKA_*.gff pphyt.prokka.v3.gff
#cp pphex/PROKKA_*.gff pphex.prokka.v3.gff
#cp pxeno/PROKKA_*.gff pxeno.prokka.v3.gff
#cp pfung/PROKKA_*.gff pfung.prokka.v3.gff
#cp ptere/PROKKA_*.gff ptere.prokka.v3.gff
#cp pcale/PROKKA_*.gff pcale.prokka.v3.gff
#cp ptera/PROKKA_*.gff ptera.prokka.v3.gff
#cp pphem/PROKKA_*.gff pphem.prokka.v3.gff
#cp psart/PROKKA_*.gff psart.prokka.v3.gff
#cp pmega/PROKKA_*.gff pmega.prokka.v3.gff
#cp pspre/PROKKA_*.gff pspre.prokka.v3.gff

#cp pphym/PROKKA_*.ffn pphym.prokka.v3.ffn
#cp pphyt/PROKKA_*.ffn pphyt.prokka.v3.ffn
#cp pphex/PROKKA_*.ffn pphex.prokka.v3.ffn
#cp pxeno/PROKKA_*.ffn pxeno.prokka.v3.ffn
#cp pfung/PROKKA_*.ffn pfung.prokka.v3.ffn
#cp ptere/PROKKA_*.ffn ptere.prokka.v3.ffn
#cp pcale/PROKKA_*.ffn pcale.prokka.v3.ffn
#cp ptera/PROKKA_*.ffn ptera.prokka.v3.ffn
#cp pphem/PROKKA_*.ffn pphem.prokka.v3.ffn
#cp psart/PROKKA_*.ffn psart.prokka.v3.ffn
#cp pmega/PROKKA_*.ffn pmega.prokka.v3.ffn
#cp pspre/PROKKA_*.ffn pspre.prokka.v3.ffn

#cp pphym/PROKKA_*.faa pphym.prokka.v3.faa
#cp pphyt/PROKKA_*.faa pphyt.prokka.v3.faa
#cp pphex/PROKKA_*.faa pphex.prokka.v3.faa
#cp pxeno/PROKKA_*.faa pxeno.prokka.v3.faa
#cp pfung/PROKKA_*.faa pfung.prokka.v3.faa
#cp ptere/PROKKA_*.faa ptere.prokka.v3.faa
#cp pcale/PROKKA_*.faa pcale.prokka.v3.faa
#cp ptera/PROKKA_*.faa ptera.prokka.v3.faa
#cp pphem/PROKKA_*.faa pphem.prokka.v3.faa
#cp psart/PROKKA_*.faa psart.prokka.v3.faa
#cp pmega/PROKKA_*.faa pmega.prokka.v3.faa
#cp pspre/PROKKA_*.faa pspre.prokka.v3.faa

