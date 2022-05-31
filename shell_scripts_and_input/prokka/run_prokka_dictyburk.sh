#!/bin/bash
prokka --force --outdir baqs159 --locustag PAGRI --proteins Burkholderia_pseudomallei_K96243_132.gbk --gcode 11 --genus Paraburkholderia --species agricolaris --strain BaQS159 --gram neg --addgenes --rnammer GCF_009455635.1_ASM945563v1_genomic.fna
 
prokka --force --outdir bhqs11 --locustag PHAYL --proteins Burkholderia_pseudomallei_K96243_132.gbk --gcode 11 --genus Paraburkholderia --species hayleyella --strain BhQS11 --gram neg --addgenes --rnammer GCF_009455685.1_ASM945568v1_genomic.fna
 
prokka --force --outdir bbqs859 --locustag PBONN --proteins Burkholderia_pseudomallei_K96243_132.gbk --gcode 11 --genus Paraburkholderia --species bonniea --strain BbQS859 --gram neg --addgenes --rnammer GCF_009455625.1_ASM945562v1_genomic.fna
# --rnammer now added

#cp baqs159/PROKKA_*.gff baqs159.prokka.v3.gff
#cp bhqs11/PROKKA_*.gff bhqs11.prokka.v3.gff
#cp bbqs859/PROKKA_*.gff bbqs859.prokka.v3.gff

#cp baqs159/PROKKA_*.ffn baqs159.prokka.v3.ffn
#cp bhqs11/PROKKA_*.ffn bhqs11.prokka.v3.ffn
#cp bbqs859/PROKKA_*.ffn bbqs859.prokka.v3.ffn

#cp baqs159/PROKKA_*.faa baqs159.prokka.v3.faa
#cp bhqs11/PROKKA_*.faa bhqs11.prokka.v3.faa
#cp bbqs859/PROKKA_*.faa bbqs859.prokka.v3.faa

#prokka --force --outdir baqs159 --locustag PAGRI --proteins GCF_000020045.1_ASM2004v1_genomic.gbff --proteins GCF_000020125.1_ASM2012v1_genomic.gbff --proteins GCF_000300095.1_ASM30009v1_genomic.gbff --proteins GCF_000756045.1_ASM75604v1_genomic.gbff --proteins GCF_000961515.1_ASM96151v1_genomic.gbff --proteins GCF_002902925.1_ASM290292v1_genomic.gbff --proteins GCF_003330745.1_ASM333074v1_genomic.gbff --proteins GCF_003330825.1_ASM333082v1_genomic.gbff --gcode 11 --genus Paraburkholderia --species agricolaris --strain BaQS159 --gram neg --addgenes --rnammer GCF_009455635.1_ASM945563v1_genomic.fna
#prokka --force --outdir bhqs11 --locustag PHAYL --proteins GCF_000020045.1_ASM2004v1_genomic.gbff --proteins GCF_000020125.1_ASM2012v1_genomic.gbff --proteins GCF_000300095.1_ASM30009v1_genomic.gbff --proteins GCF_000756045.1_ASM75604v1_genomic.gbff --proteins GCF_000961515.1_ASM96151v1_genomic.gbff --proteins GCF_002902925.1_ASM290292v1_genomic.gbff --proteins GCF_003330745.1_ASM333074v1_genomic.gbff --proteins GCF_003330825.1_ASM333082v1_genomic.gbff --gcode 11 --genus Paraburkholderia --species hayleyella --strain BhQS11 --gram neg --addgenes --rnammer GCF_009455685.1_ASM945568v1_genomic.fna
#prokka --force --outdir bbqs859 --locustag PBONN --proteins GCF_000020045.1_ASM2004v1_genomic.gbff --proteins GCF_000020125.1_ASM2012v1_genomic.gbff --proteins GCF_000300095.1_ASM30009v1_genomic.gbff --proteins GCF_000756045.1_ASM75604v1_genomic.gbff --proteins GCF_000961515.1_ASM96151v1_genomic.gbff --proteins GCF_002902925.1_ASM290292v1_genomic.gbff --proteins GCF_003330745.1_ASM333074v1_genomic.gbff --proteins GCF_003330825.1_ASM333082v1_genomic.gbff --gcode 11 --genus Paraburkholderia --species bonniea --strain BbQS859 --gram neg --addgenes --rnammer GCF_009455625.1_ASM945562v1_genomic.fna 

