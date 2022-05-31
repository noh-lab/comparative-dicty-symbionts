#### roary pangenome and cog summarize and visualize ####

#setwd("/personal/snoh/comparative_burk/")
library(ggplot2)
library(reshape2)
library(cluster)
library(vegan)
library(ggrepel)
library(edgeR)

# useful for formatting figures, adjust parameters as necessary
manuscript_theme = theme_classic() + theme(axis.line.x = element_line(color="black", size = .7),axis.line.y = element_line(color="black", size = .7),text=element_text(size=14))
manuscript_theme_grid = theme_classic() + theme(axis.line.x = element_line(color="black", size = .7),axis.line.y = element_line(color="black", size = .7),text=element_text(size=14),axis.text.y = element_text(size=8))
# + theme(axis.text.x = element_text(size=8))

#### processing of raw datafiles used for analyses and visualization; can skip to analyses below ####
# metadata regarding genome size, GC%, number of intact genes, number of pseudogenes
metadata <- data.frame(status = factor(c("dicty","dicty","dicty","free-living","symbiotic","symbiotic","free-living","symbiotic","symbiotic","symbiotic","free-living","symbiotic","free-living","free-living","free-living")),
                       size = c(8721420, 4098182, 4125700, 7211159, 9058983, 7627590, 8597887, 7651131, 8676562, 8214658, 5930529, 7829542, 7118039, 10062489, 9702951), 
                       GC = c(0.616337, 0.58725, 0.592403, 0.619327, 0.617501, 0.620555, 0.623499, 0.631452, 0.622924, 0.622868, 0.635233, 0.632082, 0.636686, 0.619184, 0.626334),
                       intact = c(7644, 3434, 3495, 6264, 7996, 6609, 7581, 6608, 7504, 7257, 5258, 6782, 6228, 8833, 8389),
                       pseudo = c(579, 265, 315, 631, 746, 372, 597, 547, 819, 520, 331, 857, 524, 803, 845))

rownames(metadata) <- c("pagri","pbonn","phayl","pcale", "pfung","pmega", "pphem", "pphex", "pphym", "pphyt", "psart", "pspre", "ptera", "ptere","pxeno") 

meta_table <- data.frame(size = metadata$size, 
                         GC= metadata$GC, 
                         intact = metadata$intact / (metadata$intact + metadata$pseudo) * 100, 
                         pseudo = metadata$pseudo / (metadata$intact + metadata$pseudo) * 100
)

write.table(meta_table, "all.genome_metadata.txt", append = FALSE, sep = "\t", row.names = TRUE, col.names = TRUE, quote=F)

# COG functional category names for visualization
translate <- data.frame(COG = factor(c("E", "G", "D", "N", "M", "B", "H", "Z", "V", "C", "W", "S", "R", "P", "U", "I", "X", "Y", "F", "O", "A", "L", "Q", "T", "K", "J"), levels=c("J", "A", "K", "L", "B", "D", "Y", "V", "T", "M", "N", "Z", "W", "U", "O", "X", "C", "G", "E", "F", "H", "I", "P", "Q", "R", "S")),
                        description = c("Amino acid transport and metabolism", "Carbohydrate transport and metabolism", "Cell cycle control, cell division, chromosome partitioning", "Cell motility", "Cell wall/ membrane/ envelope biogenesis", "Chromatin structure and dynamics", "Coenzyme transport and metabolism", "Cytoskeleton", "Defense mechanisms", "Energy production and conversion", "Extracellular structures", "Function unknown", "General function prediction only", "Inorganic ion transport and metabolism", "Intracellular trafficking, secretion, and vesicular transport", "Lipid transport and metabolism", "Mobilome: prophages, transposons", "Nuclear structure", "Nucleotide transport and metabolism", "Posttranslational modification, protein turnover, chaperones", "RNA processing and modification", "Replication, recombination and repair", "Secondary metabolites biosynthesis, transport and catabolism", "Signal transduction mechanisms", "Transcription", "Translation, ribosomal structure and biogenesis"),
                        stringsAsFactors=FALSE) 

cogname <- read.delim("cognames2003-2014.tab", h=T, sep="\t")

# genes categorized into cog for each genome
pagri <- read.delim("baqs159.cog.intact.categorized", h=F)
phayl <- read.delim("bhqs11.cog.intact.categorized", h=F)
pbonn <- read.delim("bbqs859.cog.intact.categorized", h=F)
pcale <- read.delim("pcale.cog.intact.categorized", h=F)
pfung <- read.delim("pfung.cog.intact.categorized", h=F)
pmega <- read.delim("pmega.cog.intact.categorized", h=F)
pphem <- read.delim("pphem.cog.intact.categorized", h=F)
pphex <- read.delim("pphex.cog.intact.categorized", h=F)
pphym <- read.delim("pphym.cog.intact.categorized", h=F)
pphyt <- read.delim("pphyt.cog.intact.categorized", h=F)
psart <- read.delim("psart.cog.intact.categorized", h=F)
pspre <- read.delim("pspre.cog.intact.categorized", h=F)
ptera <- read.delim("ptera.cog.intact.categorized", h=F)
ptere <- read.delim("ptere.cog.intact.categorized", h=F)
pxeno <- read.delim("pxeno.cog.intact.categorized", h=F)

# roary pan genome analysis results
pan <- read.delim("12cont_intact.gene_presence_absence.csv", sep=",", h=T, row.names = 1)
binary <- read.delim("12cont_intact.gene_presence_absence.Rtab", sep="\t", h=T, row.names = 1)
table(rownames(pan)==rownames(binary))

# mark genes based on presence in dicty-burk genomes
pan$dicty <- paste0(binary$baqs159.prokka.intact, binary$bhqs11.prokka.intact, binary$bbqs859.prokka.intact)
table(pan$dicty)

# mark genes based on presence in all examined paraburk genomes (core genes)
pan$core <- paste0(binary$baqs159.prokka.intact, binary$bhqs11.prokka.intact, binary$bbqs859.prokka.intact, binary$pcale.prokka.intact, binary$pfung.prokka.intact, binary$pmega.prokka.intact, binary$pphem.prokka.intact, binary$pphex.prokka.intact, binary$pphym.prokka.intact, binary$pphyt.prokka.intact, binary$psart.prokka.intact, binary$pspre.prokka.intact, binary$ptera.prokka.intact, binary$ptere.prokka.intact, binary$pxeno.prokka.intact)

dicty <- subset(pan, dicty=="111")[,c(14:16,2)] # but including core
core <- subset(pan, core=="111111111111111")[,c(14:28,2)] 

# separate dicty shared genes exclusive of core
table(rownames(dicty) %in% rownames(core))
temp <- dicty[!(rownames(dicty) %in% rownames(core)),]
dicty <- temp # 304 that are not in core; may be present in subset of other genomes

# note that some are the not-single-hit as more than one gene listed within a genome
grep("\t", core$baqs159.prokka.intact)
grep("\t", core$bbqs859.prokka.intact)
grep("\t", core$bhqs11.prokka.intact)
# for example some of these are consecutive gene ids
check <- grep("\t", core$bbqs859.prokka.intact)
core[check,"bbqs859.prokka.intact"] 
rm(check)

# separate any multi-hit core genes
core.a <- tidyr::separate_rows(core, baqs159.prokka.intact, sep="\t")[,1]
core.b <- tidyr::separate_rows(core, bbqs859.prokka.intact, sep="\t")[,2]
core.h <- tidyr::separate_rows(core, bhqs11.prokka.intact, sep="\t")[,3]
core.1 <- tidyr::separate_rows(core, pcale.prokka.intact, sep="\t")[,4]
core.f <- tidyr::separate_rows(core, pfung.prokka.intact, sep="\t")[,5]
core.2 <- tidyr::separate_rows(core, pmega.prokka.intact, sep="\t")[,6]

core.3 <- tidyr::separate_rows(core, pphem.prokka.intact, sep="\t")[,7]
core.4 <- tidyr::separate_rows(core, pphex.prokka.intact, sep="\t")[,8]
core.5 <- tidyr::separate_rows(core, pphym.prokka.intact, sep="\t")[,9]
core.6 <- tidyr::separate_rows(core, pphyt.prokka.intact, sep="\t")[,10]
core.7 <- tidyr::separate_rows(core, psart.prokka.intact, sep="\t")[,11]
core.s <- tidyr::separate_rows(core, pspre.prokka.intact, sep="\t")[,12]
core.8 <- tidyr::separate_rows(core, ptera.prokka.intact, sep="\t")[,13]

core.t <- tidyr::separate_rows(core, ptere.prokka.intact, sep="\t")[,14]
core.x <- tidyr::separate_rows(core, pxeno.prokka.intact, sep="\t")[,15]

# within each genome, mark which genes are core genes
pagri$core <- pagri$V1 %in% core.a$baqs159.prokka.intact
pbonn$core <- pbonn$V1 %in% core.b$bbqs859.prokka.intact
phayl$core <- phayl$V1 %in% core.h$bhqs11.prokka.intact
pcale$core <- pcale$V1 %in% core.1$pcale.prokka.intact
pfung$core <- pfung$V1 %in% core.f$pfung.prokka.intact
pmega$core <- pmega$V1 %in% core.2$pmega.prokka.intact

pphem$core <- pphem$V1 %in% core.3$pphem.prokka.intact
pphex$core <- pphex$V1 %in% core.4$pphex.prokka.intact
pphym$core <- pphym$V1 %in% core.5$pphym.prokka.intact
pphyt$core <- pphyt$V1 %in% core.6$pphyt.prokka.intact
psart$core <- psart$V1 %in% core.7$psart.prokka.intact
pspre$core <- pspre$V1 %in% core.s$pspre.prokka.intact
ptera$core <- ptera$V1 %in% core.8$ptera.prokka.intact

ptere$core <- ptere$V1 %in% core.t$ptere.prokka.intact
pxeno$core <- pxeno$V1 %in% core.x$pxeno.prokka.intact

# separate any multi-hit dicty-shared genes
dicty.a <- tidyr::separate_rows(dicty, baqs159.prokka.intact, sep="\t")[,1]
dicty.h <- tidyr::separate_rows(dicty, bhqs11.prokka.intact, sep="\t")[,3]
dicty.b <- tidyr::separate_rows(dicty, bbqs859.prokka.intact, sep="\t")[,2]

# within each genome, mark which genes are dicty genes
pagri$dicty <- pagri$V1 %in% dicty.a$baqs159.prokka.intact
pbonn$dicty <- pbonn$V1 %in% dicty.b$bbqs859.prokka.intact
phayl$dicty <- phayl$V1 %in% dicty.h$bhqs11.prokka.intact
pcale$dicty <- "FALSE"
pfung$dicty <- "FALSE"
pmega$dicty <- "FALSE"

pphem$dicty <- "FALSE"
pphex$dicty <- "FALSE"
pphym$dicty <- "FALSE"
pphyt$dicty <- "FALSE"
psart$dicty <- "FALSE"
pspre$dicty <- "FALSE"
ptera$dicty <- "FALSE"

ptere$dicty <- "FALSE"
pxeno$dicty <- "FALSE"


# partition annotated genes within each genome by COG category and into core, dicty but not core, and accessory (all others)
temp <- as.data.frame(table(pagri$V3))
colnames(temp) <- c("COG","all")
temp$core <- table(subset(pagri, pagri$core=="TRUE")[,3])
temp$dicty <- table(subset(pagri, pagri$dicty=="TRUE")[,3])
temp$accessory <- temp$all - temp$core - temp$dicty
cat.agri <- merge(x=translate, y=temp, by="COG", all=T, sort = T)

temp <- as.data.frame(table(pbonn$V3))
colnames(temp) <- c("COG","all")
temp$core <- table(subset(pbonn, pbonn$core=="TRUE")[,3])
temp$dicty <- table(subset(pbonn, pbonn$dicty=="TRUE")[,3])
temp$accessory <- temp$all - temp$core - temp$dicty
cat.bonn <- merge(x=translate, y=temp, by="COG", all=T, sort = T)

temp <- as.data.frame(table(phayl$V3))
colnames(temp) <- c("COG","all")
temp$core <- table(subset(phayl, phayl$core=="TRUE")[,3])
temp$dicty <- table(subset(phayl, phayl$dicty=="TRUE")[,3])
temp$accessory <- temp$all - temp$core - temp$dicty
cat.hayl <- merge(x=translate, y=temp, by="COG", all=T, sort = T)

temp <- as.data.frame(table(pcale$V3))
colnames(temp) <- c("COG","all")
temp$core <- table(subset(pcale, pcale$core=="TRUE")[,3])
temp$accessory <- temp$all - temp$core
cat.cale <- merge(x=translate, y=temp, by="COG", all=T, sort = T)

temp <- as.data.frame(table(pfung$V3))
colnames(temp) <- c("COG","all")
temp$core <- table(subset(pfung, pfung$core=="TRUE")[,3])
temp$accessory <- temp$all - temp$core
cat.fung <- merge(x=translate, y=temp, by="COG", all=T, sort = T)

temp <- as.data.frame(table(pmega$V3))
colnames(temp) <- c("COG","all")
temp$core <- table(subset(pmega, pmega$core=="TRUE")[,3])
temp$accessory <- temp$all - temp$core
cat.mega <- merge(x=translate, y=temp, by="COG", all=T, sort = T)

temp <- as.data.frame(table(pphem$V3))
colnames(temp) <- c("COG","all")
temp$core <- table(subset(pphem, pphem$core=="TRUE")[,3])
temp$accessory <- temp$all - temp$core
cat.phem <- merge(x=translate, y=temp, by="COG", all=T, sort = T)

temp <- as.data.frame(table(pphex$V3))
colnames(temp) <- c("COG","all")
temp$core <- table(subset(pphex, pphex$core=="TRUE")[,3])
temp$accessory <- temp$all - temp$core
cat.phex <- merge(x=translate, y=temp, by="COG", all=T, sort = T)

temp <- as.data.frame(table(pphym$V3))
colnames(temp) <- c("COG","all")
temp$core <- table(subset(pphym, pphym$core=="TRUE")[,3])
temp$accessory <- temp$all - temp$core
cat.phym <- merge(x=translate, y=temp, by="COG", all=T, sort = T)

temp <- as.data.frame(table(pphyt$V3))
colnames(temp) <- c("COG","all")
temp$core <- table(subset(pphyt, pphyt$core=="TRUE")[,3])
temp$accessory <- temp$all - temp$core
cat.phyt <- merge(x=translate, y=temp, by="COG", all=T, sort = T)

temp <- as.data.frame(table(psart$V3))
colnames(temp) <- c("COG","all")
temp$core <- table(subset(psart, psart$core=="TRUE")[,3])
temp$accessory <- temp$all - temp$core
cat.sart <- merge(x=translate, y=temp, by="COG", all=T, sort = T)

temp <- as.data.frame(table(pspre$V3))
colnames(temp) <- c("COG","all")
temp$core <- table(subset(pspre, pspre$core=="TRUE")[,3])
temp$accessory <- temp$all - temp$core
cat.spre <- merge(x=translate, y=temp, by="COG", all=T, sort = T)[1:26,]

temp <- as.data.frame(table(ptera$V3))
colnames(temp) <- c("COG","all")
temp$core <- table(subset(ptera, ptera$core=="TRUE")[,3])
temp$accessory <- temp$all - temp$core
cat.tera <- merge(x=translate, y=temp, by="COG", all=T, sort = T)[1:26,]

temp <- as.data.frame(table(ptere$V3))
colnames(temp) <- c("COG","all")
temp$core <- table(subset(ptere, ptere$core=="TRUE")[,3])
temp$accessory <- temp$all - temp$core
cat.tere <- merge(x=translate, y=temp, by="COG", all=T, sort = T)

temp <- as.data.frame(table(pxeno$V3))
colnames(temp) <- c("COG","all")
temp$core <- table(subset(pxeno, pxeno$core=="TRUE")[,3])
temp$accessory <- temp$all - temp$core
cat.xeno <- merge(x=translate, y=temp, by="COG", all=T, sort = T)

# tidy up
rm(core.a, core.b, core.h, core.1, core.2, core.3, core.4, core.5, core.6, core.7, core.8, core.f, core.s, core.t, core.x, dicty.a, dicty.b, dicty.h)
rm(binary, core, dicty, pan)

# save the processed (categorized) files  
write.table(cat.agri[,c(1,3:6)], "agri.cog.cog_categories.txt", append = FALSE, sep = "\t", row.names = F, col.names = TRUE, quote=F)
write.table(cat.bonn[,c(1,3:6)], "bonn.cog.cog_categories.txt", append = FALSE, sep = "\t", row.names = F, col.names = TRUE, quote=F)
write.table(cat.hayl[,c(1,3:6)], "hayl.cog.cog_categories.txt", append = FALSE, sep = "\t", row.names = F, col.names = TRUE, quote=F)
write.table(cat.cale[,c(1,3:5)], "cale.cog.cog_categories.txt", append = FALSE, sep = "\t", row.names = F, col.names = TRUE, quote=F)
write.table(cat.fung[,c(1,3:5)], "fung.cog.cog_categories.txt", append = FALSE, sep = "\t", row.names = F, col.names = TRUE, quote=F)
write.table(cat.mega[,c(1,3:5)], "mega.cog.cog_categories.txt", append = FALSE, sep = "\t", row.names = F, col.names = TRUE, quote=F)
write.table(cat.phem[,c(1,3:5)], "phem.cog.cog_categories.txt", append = FALSE, sep = "\t", row.names = F, col.names = TRUE, quote=F)
write.table(cat.phex[,c(1,3:5)], "phex.cog.cog_categories.txt", append = FALSE, sep = "\t", row.names = F, col.names = TRUE, quote=F)
write.table(cat.phym[,c(1,3:5)], "phym.cog.cog_categories.txt", append = FALSE, sep = "\t", row.names = F, col.names = TRUE, quote=F)
write.table(cat.phyt[,c(1,3:5)], "phyt.cog.cog_categories.txt", append = FALSE, sep = "\t", row.names = F, col.names = TRUE, quote=F)
write.table(cat.sart[,c(1,3:5)], "sart.cog.cog_categories.txt", append = FALSE, sep = "\t", row.names = F, col.names = TRUE, quote=F)
write.table(cat.spre[,c(1,3:5)], "spre.cog.cog_categories.txt", append = FALSE, sep = "\t", row.names = F, col.names = TRUE, quote=F)
write.table(cat.tera[,c(1,3:5)], "tera.cog.cog_categories.txt", append = FALSE, sep = "\t", row.names = F, col.names = TRUE, quote=F)
write.table(cat.tere[,c(1,3:5)], "tere.cog.cog_categories.txt", append = FALSE, sep = "\t", row.names = F, col.names = TRUE, quote=F)
write.table(cat.xeno[,c(1,3:5)], "xeno.cog.cog_categories.txt", append = FALSE, sep = "\t", row.names = F, col.names = TRUE, quote=F)


# compile cog category counts into single data frame
cat.comp <- cbind(cat.agri[,3],cat.bonn[,3],cat.hayl[,3],cat.cale[,3],cat.fung[,3],cat.mega[,3],cat.phem[,3],cat.phex[,3],cat.phym[,3],cat.phyt[,3],cat.sart[,3],cat.spre[,3],cat.tera[,3],cat.tere[,3],cat.xeno[,3])
df.cat <- as.data.frame(t(cat.comp))
rownames(df.cat) <- c("pagri","pbonn","phayl","pcale", "pfung","pmega", "pphem", "pphex", "pphym", "pphyt", "psart", "pspre", "ptera", "ptere","pxeno")
colnames(df.cat) <- cat.agri$COG

df.cat[is.na(df.cat)] <- 0
df.cat[,1:26] <- lapply(df.cat[,1:26], as.numeric)

write.table(df.cat, "all.cog.cog_categories.txt", append = FALSE, sep = "\t", row.names = TRUE, col.names = TRUE, quote=F)

# tidy up
rm(cat.agri, cat.bonn, cat.hayl, cat.cale, cat.fung, cat.mega, cat.phem, cat.phex, cat.phym, cat.phyt, cat.sart, cat.spre, cat.tera, cat.tere, cat.xeno)


# fine COGs with all control burk
# also enable comparison by individual cog numbers rather than COG categories
temp <- as.data.frame(table(pagri$V2))
colnames(temp) <- c("COG.","genes")
cog.agri <- merge(x=cogname[,c(1,3)], y=temp, by="COG.", all=T, sort = T)

temp <- as.data.frame(table(pbonn$V2))
colnames(temp) <- c("COG.","genes")
cog.bonn <- merge(x=cogname[,c(1,3)], y=temp, by="COG.", all=T, sort = T)

temp <- as.data.frame(table(phayl$V2))
colnames(temp) <- c("COG.","genes")
cog.hayl <- merge(x=cogname[,c(1,3)], y=temp, by="COG.", all=T, sort = T)

temp <- as.data.frame(table(pfung$V2))
colnames(temp) <- c("COG.","genes")
cog.fung <- merge(x=cogname[,c(1,3)], y=temp, by="COG.", all=T, sort = T)

temp <- as.data.frame(table(pmega$V2))
colnames(temp) <- c("COG.","genes")
cog.mega <- merge(x=cogname[,c(1,3)], y=temp, by="COG.", all=T, sort = T)

temp <- as.data.frame(table(ptere$V2))
colnames(temp) <- c("COG.","genes")
cog.tere <- merge(x=cogname[,c(1,3)], y=temp, by="COG.", all=T, sort = T)

temp <- as.data.frame(table(pxeno$V2))
colnames(temp) <- c("COG.","genes")
cog.xeno <- merge(x=cogname[,c(1,3)], y=temp, by="COG.", all=T, sort = T)

temp <- as.data.frame(table(pphex$V2))
colnames(temp) <- c("COG.","genes")
cog.phex <- merge(x=cogname[,c(1,3)], y=temp, by="COG.", all=T, sort = T)

temp <- as.data.frame(table(pphym$V2))
colnames(temp) <- c("COG.","genes")
cog.phym <- merge(x=cogname[,c(1,3)], y=temp, by="COG.", all=T, sort = T)

temp <- as.data.frame(table(pphyt$V2))
colnames(temp) <- c("COG.","genes")
cog.phyt <- merge(x=cogname[,c(1,3)], y=temp, by="COG.", all=T, sort = T)

temp <- as.data.frame(table(pspre$V2))
colnames(temp) <- c("COG.","genes")
cog.spre <- merge(x=cogname[,c(1,3)], y=temp, by="COG.", all=T, sort = T)

temp <- as.data.frame(table(pcale$V2))
colnames(temp) <- c("COG.","genes")
cog.cale <- merge(x=cogname[,c(1,3)], y=temp, by="COG.", all=T, sort = T)

temp <- as.data.frame(table(pphem$V2))
colnames(temp) <- c("COG.","genes")
cog.phem <- merge(x=cogname[,c(1,3)], y=temp, by="COG.", all=T, sort = T)

temp <- as.data.frame(table(psart$V2))
colnames(temp) <- c("COG.","genes")
cog.sart <- merge(x=cogname[,c(1,3)], y=temp, by="COG.", all=T, sort = T)

temp <- as.data.frame(table(ptera$V2))
colnames(temp) <- c("COG.","genes")
cog.tera <- merge(x=cogname[,c(1,3)], y=temp, by="COG.", all=T, sort = T)

cog.comp <- cbind(cog.agri[,3],cog.bonn[,3],cog.hayl[,3],cog.cale[,3],cog.fung[,3],cog.mega[,3],cog.phem[,3],cog.phex[,3],cog.phym[,3],cog.phyt[,3],cog.sart[,3],cog.spre[1:4631,3],cog.tera[1:4631,3],cog.tere[,3],cog.xeno[,3])
df.cog <- as.data.frame(t(cog.comp))
rownames(df.cog) <- c("pagri","pbonn","phayl","pcale", "pfung","pmega", "pphem", "pphex", "pphym", "pphyt", "psart", "pspre", "ptera", "ptere","pxeno")
colnames(df.cog) <- cog.agri$COG.

df.cog[is.na(df.cog)] <- 0
df.cog[,1:ncol(df.cog)] <- lapply(df.cog[,1:ncol(df.cog)], as.numeric)

# remove any columns without counts
temp <- df.cog[, colSums(df.cog != 0) > 0]
df.cog <- temp
ncol(df.cog)

write.table(df.cog, "all.cog.cog_numbers.txt", append = FALSE, sep = "\t", row.names = TRUE, col.names = TRUE, quote=F)

rm(pagri, pbonn, phayl, pcale, pfung, pmega, pphem, pphex, pphym, pphyt, psart, pspre, ptera, ptere, pxeno)
rm(cog.agri,cog.bonn,cog.hayl,cog.cale,cog.fung,cog.mega,cog.phem,cog.phex,cog.phym,cog.phyt,cog.sart,cog.spre,cog.tera,cog.tere,cog.xeno)



### broad visualization of cog categories present in each genome couched in pangenome analysis ####

cat.agri <- read.table("agri.cog.cog_categories.txt", h=T)
cat.hayl <- read.table("bonn.cog.cog_categories.txt", h=T)
cat.bonn <- read.table("hayl.cog.cog_categories.txt", h=T)
cat.fung <- read.table("fung.cog.cog_categories.txt", h=T)
cat.spre <- read.table("spre.cog.cog_categories.txt", h=T)
cat.tere <- read.table("tere.cog.cog_categories.txt", h=T)
cat.xeno <- read.table("xeno.cog.cog_categories.txt", h=T)

# to visualize cog categories, first make longform tables
lcat.agri <- subset(melt(cat.agri), variable!="all")
colnames(lcat.agri) <- c("COG","description","category","count")
lcat.agri$category <- factor(lcat.agri$category, levels=c("accessory","dicty","core"))
lcat.agri$genome <- "pagri"

lcat.hayl <- subset(melt(cat.hayl), variable!="all")
colnames(lcat.hayl) <- c("COG","description","category","count")
lcat.hayl$category <- factor(lcat.hayl$category, levels=c("accessory","dicty","core"))
lcat.hayl$genome <- "phayl"

lcat.bonn <- subset(melt(cat.bonn), variable!="all")
colnames(lcat.bonn) <- c("COG","description","category","count")
lcat.bonn$category <- factor(lcat.bonn$category, levels=c("accessory","dicty","core"))
lcat.bonn$genome <- "pbonn"

lcat.fung <- subset(melt(cat.fung), variable!="all")
colnames(lcat.fung) <- c("COG","description","category","count")
lcat.fung$category <- factor(lcat.fung$category, levels=c("accessory","core"))
lcat.fung$genome <- "pfung"

lcat.spre <- subset(melt(cat.spre), variable!="all")
colnames(lcat.spre) <- c("COG","description","category","count")
lcat.spre$category <- factor(lcat.spre$category, levels=c("accessory","core"))
lcat.spre$genome <- "pspre"

lcat.tere <- subset(melt(cat.tere), variable!="all")
colnames(lcat.tere) <- c("COG","description","category","count")
lcat.tere$category <- factor(lcat.tere$category, levels=c("accessory","core"))
lcat.tere$genome <- "ptere"

lcat.xeno <- subset(melt(cat.xeno), variable!="all")
colnames(lcat.xeno) <- c("COG","description","category","count")
lcat.xeno$category <- factor(lcat.xeno$category, levels=c("accessory","core"))
lcat.xeno$genome <- "pxeno"


# shades of "#191100"
ggplot(lcat.agri, aes(x=COG, y=count, fill=category)) + geom_bar(stat="identity") + scale_x_discrete(breaks=lcat.agri$COG, limits = rev(levels(lcat.agri$COG)))  +coord_flip() + scale_fill_manual(values=c("#d1cfcc","#8c887f","#191100")) + manuscript_theme_grid # all, dicty but not core, core , labels=lcat.agri$description,expand=c(0, 0.2)

cog1 <- ggplot(lcat.agri, aes(x=COG, y=count, fill=category)) + geom_bar(stat="identity") + scale_x_discrete(breaks=lcat.agri$COG, expand=c(0, 0.2), limits = rev(levels(lcat.agri$COG)))  +coord_flip() + scale_fill_manual(values=c("#d1cfcc","#b3cc31","#191100")) + manuscript_theme_grid + ylim(0,650) + ggtitle("P. agricolaris BaQS159")

cog2 <- ggplot(lcat.bonn, aes(x=COG, y=count, fill=category)) + geom_bar(stat="identity") + scale_x_discrete(breaks=lcat.bonn$COG, expand=c(0, 0.2), limits = rev(levels(lcat.bonn$COG)))  +coord_flip() + scale_fill_manual(values=c("#d1cfcc","#b3cc31","#191100")) + manuscript_theme_grid + ylim(0,650) + ggtitle("P. bonniea BbQS859")

cog3 <- ggplot(lcat.hayl, aes(x=COG, y=count, fill=category)) + geom_bar(stat="identity") + scale_x_discrete(breaks=lcat.hayl$COG, expand=c(0, 0.2), limits = rev(levels(lcat.hayl$COG)))  +coord_flip() + scale_fill_manual(values=c("#d1cfcc","#b3cc31","#191100")) + manuscript_theme_grid + ylim(0,650) + ggtitle("P. hayleyella BhQS11")

gridExtra::grid.arrange(cog1, cog2, cog3, ncol=1)


cog4 <- ggplot(lcat.fung, aes(x=COG, y=count, fill=category)) + geom_bar(stat="identity") + scale_x_discrete(breaks=lcat.fung$COG, expand=c(0, 0.2), limits = rev(levels(lcat.fung$COG)))  +coord_flip() + scale_fill_manual(values=c("#d1cfcc","#191100")) + manuscript_theme_grid + ggtitle("P. fungorum ATCC BAA-463")

cog5 <- ggplot(lcat.spre, aes(x=COG, y=count, fill=category)) + geom_bar(stat="identity") + scale_x_discrete(breaks=lcat.spre$COG, expand=c(0, 0.2), limits = rev(levels(lcat.spre$COG)))  +coord_flip() + scale_fill_manual(values=c("#d1cfcc","#191100")) + manuscript_theme_grid + ggtitle("P. sprentiae WSM5005")

cog6 <- ggplot(lcat.tere, aes(x=COG, y=count, fill=category)) + geom_bar(stat="identity") + scale_x_discrete(breaks=lcat.tere$COG, expand=c(0, 0.2), limits = rev(levels(lcat.tere$COG)))  +coord_flip() + scale_fill_manual(values=c("#d1cfcc","#191100")) + manuscript_theme_grid + ggtitle("P. terrae DSM17804")

cog7 <- ggplot(lcat.xeno, aes(x=COG, y=count, fill=category)) + geom_bar(stat="identity") + scale_x_discrete(breaks=lcat.xeno$COG, expand=c(0, 0.2), limits = rev(levels(lcat.xeno$COG)))  +coord_flip() + scale_fill_manual(values=c("#d1cfcc","#191100")) + manuscript_theme_grid + ggtitle("P. xenovorans LB400")

gridExtra::grid.arrange(cog1, cog4, cog5, cog6, cog7, ncol=2)


cairo_ps(filename = paste("cog_cat.dicty_burk",format(Sys.time(),"%Y%m%d"),"eps",sep="."), height=10)
gridExtra::grid.arrange(cog1, cog2, cog3, ncol=1)
dev.off()

cairo_ps(filename = paste("cog_cat.agri_control",format(Sys.time(),"%Y%m%d"),"eps",sep="."), height=10)
gridExtra::grid.arrange(cog1, cog4, cog5, cog6, cog7, ncol=2)
dev.off()


#### agglomerative clustering using agnes ####

# read in data
meta_table <- read.table("all.genome_metadata.txt", row.names=1, sep="\t", h=T)
df.cat <- read.table("all.cog.cog_categories.txt", row.names=1, sep="\t", h=T)


# broad COG comparison with all control burk
sc.cat <- wisconsin(df.cat)
hc.cat <- agnes(sc.cat, method = "ward")

# agglomerative coefficient
hc.cat$ac # wisconsin 0.8232733

pltree(hc.cat, cex = 1, hang = -1, main = "Dendrogram of agnes")
clust.cat <- cutree(hc.cat, k = 2)

# other than bonn/hayl reduced genomes, there is no clear difference based on ecology (symbiotic or not)

# add cluster result to metadata table for further exploration
clust.cat <- gsub("2", "reduced", gsub("1", "control", clust.cat))
meta_table$cluster <- as.factor(clust.cat)


#### NMDS ####
# basic description here - https://rpubs.com/CPEL/NMDS
# useful guide - http://dx.doi.org/10.1111/1574-6941.12437
# normal setup for ecology NMDS is sites in rows and species in columns; you typically want to know which species are different across relevant sites; this translates to genomes in rows in COG function in columns, where we want to know which COGs are different across relevant genomes
# useful guide - http://dx.doi.org/10.1111/1574-6941.12437

nms.cat <- metaMDS(df.cat, distance="bray", k=3) # consistent with hclust above
nms.cat # stress_k2 ~ 0.06, good fit; stress_k3 ~ 0.03, better fit
stressplot(nms.cat)

plot(nms.cat, type="t")
ord <- ordiellipse(nms.cat, clust.cat ,display = "sites", kind ="sd", conf = 0.95, label = F)

# extract coordinates for re-plotting
df.nms <- data.frame(nms.cat$points)
df.nms$type <- c("agri","reduced","reduced",rep("control",12))

df.plot <- data.frame(nms.cat$species)

ggplot(df.nms, aes(x=MDS1, y=MDS2)) + geom_label_repel(aes(color=type, label=row.names(df.nms), fontface="italic"), size=5) + geom_point(aes(color=type), size=3) + geom_text_repel(data=df.plot, aes(x=MDS1, y=MDS2, label=row.names(df.plot)), color="darkgrey", size=4) + scale_color_manual(values=c("#47d3d3","#7ee0e0","#d34747"))  + manuscript_theme + theme(legend.position="none") 

cairo_ps(filename = paste("cog_nmds.all",format(Sys.time(),"%Y%m%d"),"eps",sep="."))
ggplot(df.nms, aes(x=MDS1, y=MDS2)) + geom_label_repel(aes(color=type, label=row.names(df.nms), fontface="italic"), size=5) + geom_point(aes(color=type), size=3) + geom_text_repel(data=df.plot, aes(x=MDS1, y=MDS2, label=row.names(df.plot)), color="darkgrey", size=4) + scale_color_manual(values=c("#47d3d3","#7ee0e0","#d34747")) + xlim(-0.35, 0.15) + ylim(-0.1, 0.05) + manuscript_theme + theme(legend.position="none")
dev.off()

# based on clustering or NMDS results, show any relationships with metadata

df.m<-NULL
for(i in 1:4){
  tmp <- data.frame(meta_table[,i], clust.cat, colnames(meta_table)[i])
  if(is.null(df.m)){df.m <- tmp} else {df.m <- rbind(df.m, tmp)}
}
colnames(df.m)<-c("Value","cluster","Statistic")

ggplot(df.m,aes(cluster,Value,fill=cluster)) + ylab("normalized count") + geom_boxplot() + geom_jitter() + theme_bw() + facet_wrap( ~ Statistic , scales="free", ncol=2) + theme(axis.text.x = element_blank(), axis.ticks.x = element_blank()) + scale_fill_manual(values=c("#47d3d3","#d34747"))


#### Exact test to examine which COG categories are different ####
## implemented in edgeR to deal with overdispersed count data (e.g. presence or absence)
# https://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1003531

abund_table <- df.cat[colSums(df.cat)>15]
abund_trans <- t(abund_table)

exact.cat <- DGEList(counts=abund_trans, group=clust.cat)
exact.cat <- calcNormFactors(exact.cat)
normList <- cpm(exact.cat, normalized.lib.sizes=TRUE)
plotMDS(normList)
exact.cat <- estimateDisp(exact.cat)
exact.cat$common.dispersion
plotBCV(exact.cat)
test.cat <- exactTest(exact.cat)
topTags(test.cat, p.value=0.05)
res.cat <- topTags(test.cat, n=nrow(test.cat$table))$table

exact.order <- rownames(topTags(test.cat))
temp <- exact.order[!(exact.order %in% c("R","S"))]
exact.order <- temp

# normList (COGcat x genome) to normList[i,], groups, rownames(normList) (longform)
df.e<-NULL
for(i in exact.order){
  tmp <- data.frame(normList[i,], clust.cat, rep(i,15))
  if(is.null(df.e)){df.e <- tmp} else {df.e <- rbind(df.e, tmp)}
}
colnames(df.e)<-c("normCount","cluster","COG")
temp <- df.e$COG
df.e$description <- gsub("^P$", subset(translate, COG=="P")[,2], gsub("^F$", subset(translate, COG=="F")[,2], gsub("^D$", subset(translate, COG=="D")[,2], gsub("^L$", subset(translate, COG=="L")[,2], gsub("^U$", subset(translate, COG=="U")[,2], gsub("^J$", subset(translate, COG=="J")[,2], gsub("^G$", subset(translate, COG=="G")[,2], gsub("^N$", subset(translate, COG=="N")[,2], gsub("^K$", subset(translate, COG=="K")[,2], temp)))))))))

df.e$title <- paste(df.e$COG, df.e$description, sep=": ")
df.e$title <- factor(df.e$title, levels = unique(df.e$title))

# string wrap function with adjustable character target width
swrap <- function(string, nwrap=20) {
  paste(strwrap(string, width=nwrap), collapse="\n")
}
swrap <- Vectorize(swrap)

# create line breaks in COG description
df.e$title <- swrap(df.e$title, 25)

ggplot(df.e, aes(cluster,normCount,fill=cluster)) + ylab("normalized count") + geom_boxplot() + geom_jitter() + theme_bw() + facet_wrap( ~ title , scales="free", ncol=3) + theme( axis.text.x = element_blank(), axis.text.y = element_blank(), axis.ticks = element_blank()) + scale_fill_manual(values=c("#47d3d3","#d34747"))


## comparison with raw counts
df.f<-NULL
for(i in exact.order){
  tmp <- data.frame(abund_trans[i,], clust.cat, rep(i,15))
  if(is.null(df.f)){df.f <- tmp} else {df.f <- rbind(df.f, tmp)}
}
colnames(df.f)<-c("rawCount","cluster","COG")

temp <- df.f$COG
df.f$description <- gsub("^P$", subset(translate, COG=="P")[,2], gsub("^F$", subset(translate, COG=="F")[,2], gsub("^D$", subset(translate, COG=="D")[,2], gsub("^L$", subset(translate, COG=="L")[,2], gsub("^U$", subset(translate, COG=="U")[,2], gsub("^J$", subset(translate, COG=="J")[,2], gsub("^G$", subset(translate, COG=="G")[,2], gsub("^N$", subset(translate, COG=="N")[,2], gsub("^K$", subset(translate, COG=="K")[,2], temp)))))))))
df.f$title <- paste(df.f$COG, df.f$description, sep=": ")
df.f$title <- factor(df.f$title, levels = unique(df.f$title))

df.f$title <- swrap(df.f$title, 25)

ggplot(df.f,aes(cluster,rawCount,fill=cluster)) + ylab("raw count") + geom_boxplot() + geom_jitter() + theme_bw() + facet_wrap( ~ title , scales="free", ncol=3) + theme( axis.text.x = element_blank()) + scale_fill_manual(values=c("#47d3d3","#d34747"))


cairo_ps(filename = paste("cog_cat.exact",format(Sys.time(),"%Y%m%d"),"eps",sep="."))
ggplot(df.e, aes(cluster,normCount,fill=cluster)) + ylab("normalized count") + geom_boxplot() + geom_jitter() + theme_bw() + facet_wrap( ~ title , scales="free", ncol=3) + theme( axis.text.x = element_blank(), axis.text.y = element_blank(), axis.ticks = element_blank()) + scale_fill_manual(values=c("#47d3d3","#d34747"))
dev.off()

cairo_ps(filename = paste("cog_cat.raw",format(Sys.time(),"%Y%m%d"),"eps",sep="."))
ggplot(df.f,aes(cluster,rawCount,fill=cluster)) + ylab("raw count") + geom_boxplot() + geom_jitter() + theme_bw() + facet_wrap( ~ title , scales="free", ncol=3) + theme( axis.text.x = element_blank()) + scale_fill_manual(values=c("#47d3d3","#d34747"))
dev.off()



#### within each COG category of interest, which COG numbers are different ####

# read in fine cog details for each genome
df.cog <- read.table("all.cog.cog_numbers.txt", sep="\t", h=T)

# "K" transcription - less than expected in reduced genomes
k.cogs <- as.vector(subset(cogname, func=="K")[,1])
temp <- df.cog[,colnames(df.cog) %in% k.cogs]
k.cogs <- temp[colSums(temp)>15]
k.cogs <- k.cogs[,order(apply(k.cogs, 2, var), decreasing=T)]
k.list <- colnames(k.cogs)

k.plots = list()
for (i in 1:ncol(k.cogs)){ 
  ordi <- ordisurf(nms.cat, k.cogs[,i], plot=F, bs="ds")
  ordi.grid <- ordi$grid
  ordi.k <- expand.grid(x=ordi.grid$x, y=ordi.grid$y)
  ordi.k$z <- as.vector(ordi.grid$z)
  ordi.k.na <- data.frame(na.omit(ordi.k))
  
  title <- gsub(",",",\n",subset(cogname, COG.==k.list[i])[,3])
  
  p <- ggplot() + stat_contour(data = ordi.k.na, aes(x = x, y = y, z = z, colour = ..level..)) + geom_point(data=df.nms,aes(MDS1,MDS2,fill=type),pch=21,size=3) + theme_bw() + theme(legend.position="none") + ggtitle(title) + theme(plot.title = element_text(size = 10), axis.title.x=element_blank(), axis.title.y=element_blank())+ scale_fill_manual(values=c("#47d3d3","#c7f1f1","#d34747"))
  k.plots[[i]] = p
}

k.names = list()
for (i in 1:length(k.list)){
  k.names[i] = as.character(subset(cogname, COG.==k.list[i])[,3])
}
k.names 

grep("regulator", k.names)

gridExtra::grid.arrange(k.plots[[1]], k.plots[[2]], k.plots[[3]], k.plots[[4]], ncol=2)
gridExtra::grid.arrange(k.plots[[5]], k.plots[[6]], k.plots[[7]], k.plots[[8]], ncol=2)
gridExtra::grid.arrange(k.plots[[9]], k.plots[[11]], k.plots[[12]], k.plots[[13]], ncol=2)
gridExtra::grid.arrange(k.plots[[14]], k.plots[[15]], k.plots[[16]], k.plots[[17]], ncol=2)
gridExtra::grid.arrange(k.plots[[18]], k.plots[[19]], k.plots[[22]], ncol=2)
gridExtra::grid.arrange(k.plots[[24]], k.plots[[26]], k.plots[[30]], ncol=2)


# "N" Cell motility - more than expected in reduced genomes
n.cogs <- as.vector(subset(cogname, func=="N")[,1])
temp <- df.cog[,colnames(df.cog) %in% n.cogs]
n.cogs <- temp[colSums(temp)>15]
n.cogs <- n.cogs[,order(apply(n.cogs, 2, var), decreasing=T)]
n.list <- colnames(n.cogs)

n.plots = list()
for (i in 1:ncol(n.cogs)){ 
  ordi <- ordisurf(nms.cat, n.cogs[,i], plot=F, bs="ds")
  ordi.grid <- ordi$grid
  ordi.n <- expand.grid(x=ordi.grid$x, y=ordi.grid$y)
  ordi.n$z <- as.vector(ordi.grid$z)
  ordi.n.na <- data.frame(na.omit(ordi.n))
  
  title <- gsub("hook-associated","hook-associated\n", gsub(",",",\n",subset(cogname, COG.==n.list[i])[,3]))
  
  p <- ggplot() + stat_contour(data = ordi.n.na, aes(x = x, y = y, z = z, colour = ..level..)) + geom_point(data=df.nms,aes(MDS1,MDS2,fill=type),pch=21,size=3) + theme_bw() + theme(legend.position="none") + ggtitle(title) + theme(plot.title = element_text(size = 10), axis.title.x=element_blank(), axis.title.y=element_blank())+ scale_fill_manual(values=c("#47d3d3","#c7f1f1","#d34747"))
  n.plots[[i]] = p
}

n.names = list()
for (i in 1:length(n.list)){
  n.names[i] = as.character(subset(cogname, COG.==n.list[i])[,3])
}
n.names 
length(grep("flag", n.names, ignore.case=T)) # 21/24 are associated with flagella

gridExtra::grid.arrange(n.plots[[3]], n.plots[[5]], n.plots[[9]], n.plots[[10]], ncol=2)
gridExtra::grid.arrange(n.plots[[11]], n.plots[[12]], n.plots[[13]], n.plots[[15]], ncol=2)
gridExtra::grid.arrange(n.plots[[16]], n.plots[[17]], n.plots[[18]], n.plots[[19]], ncol=2)
gridExtra::grid.arrange(n.plots[[20]], n.plots[[21]], n.plots[[22]], n.plots[[23]], ncol=2)
gridExtra::grid.arrange(n.plots[[24]], ncol=2)

# 5, 9, 11, 12, 16, 24 flagellar biosynthesis all UP
# 3, 10, 13, 15 flagellar hook all UP
# 18-23 flagellar basal body all UP

gridExtra::grid.arrange(n.plots[[3]], n.plots[[10]], n.plots[[9]], n.plots[[18]], ncol=2)


# "G" Carbohydrate transport and metabolism - less than expected in reduced genomes
# e.g. ABC-type sugar transport

g.cogs <- as.vector(subset(cogname, func=="G")[,1])
temp <- df.cog[,colnames(df.cog) %in% g.cogs]
g.cogs <- temp[colSums(temp)>15]
g.cogs <- g.cogs[,order(apply(g.cogs, 2, var), decreasing=T)]
g.list <- colnames(g.cogs)

g.plots = list()
for (i in 1:ncol(g.cogs)){ 
  ordi <- ordisurf(nms.cat, g.cogs[,i], plot=F, bs="ds")
  ordi.grid <- ordi$grid
  ordi.g <- expand.grid(x=ordi.grid$x, y=ordi.grid$y)
  ordi.g$z <- as.vector(ordi.grid$z)
  ordi.g.na <- data.frame(na.omit(ordi.g))
  
  title <- gsub(",",",\n",subset(cogname, COG.==g.list[i])[,3])
  
  p <- ggplot() + stat_contour(data = ordi.g.na, aes(x = x, y = y, z = z, colour = ..level..)) + geom_point(data=df.nms,aes(MDS1,MDS2,fill=type),pch=21,size=3) + theme_bw() + theme(legend.position="none") + ggtitle(title) + theme(plot.title = element_text(size = 10), axis.title.x=element_blank(), axis.title.y=element_blank())+ scale_fill_manual(values=c("#47d3d3","#c7f1f1","#d34747"))
  g.plots[[i]] = p
}

g.names = list()
for (i in 1:length(g.list)){
  g.names[i] = as.character(subset(cogname, COG.==g.list[i])[,3])
}
g.names 
grep("sugar|glucose|fructose", g.names, ignore.case=T)
grep("ABC", g.names, ignore.case=T)


# fewer across all of these
gridExtra::grid.arrange(g.plots[[2]], g.plots[[3]], g.plots[[4]], g.plots[[5]], ncol=2)
gridExtra::grid.arrange(g.plots[[6]], g.plots[[7]], g.plots[[8]], g.plots[[9]], ncol=2)
gridExtra::grid.arrange(g.plots[[13]], g.plots[[14]], g.plots[[15]], g.plots[[21]], ncol=2)
gridExtra::grid.arrange(g.plots[[22]], g.plots[[41]], g.plots[[47]], ncol=2)


# "P" inorganic ion transport - less than expected in reduced genomes
p.cogs <- as.vector(subset(cogname, func=="P")[,1])
temp <- df.cog[,colnames(df.cog) %in% p.cogs]
p.cogs <- temp[colSums(temp)>15]
p.cogs <- p.cogs[,order(apply(p.cogs, 2, var), decreasing=T)]
p.list <- colnames(p.cogs)

p.plots = list()
for (i in 1:ncol(p.cogs)){ 
  ordi <- ordisurf(nms.cat, p.cogs[,i], plot=F, bs="ds")
  ordi.grid <- ordi$grid
  ordi.p <- expand.grid(x=ordi.grid$x, y=ordi.grid$y)
  ordi.p$z <- as.vector(ordi.grid$z)
  ordi.p.na <- data.frame(na.omit(ordi.p))
  
  title <- gsub(",",",\n",gsub("the","the\n",subset(cogname, COG.==p.list[i])[,3]))
  #title <- gsub(",",",\n",subset(cogname, COG.==p.list[i])[,3])
  
  p <- ggplot() + stat_contour(data = ordi.p.na, aes(x = x, y = y, z = z, colour = ..level..)) + geom_point(data=df.nms,aes(MDS1,MDS2,fill=type),pch=21,size=3) + theme_bw() + theme(legend.position="none") + ggtitle(title) + theme(plot.title = element_text(size = 10), axis.title.x=element_blank(), axis.title.y=element_blank())+ scale_fill_manual(values=c("#47d3d3","#c7f1f1","#d34747"))
  p.plots[[i]] = p
}

p.names = list()
for (i in 1:length(p.list)){
  p.names[i] = as.character(subset(cogname, COG.==p.list[i])[,3])
}
p.names 

grep("metal", p.names, ignore.case=T)
grep("Fe|Zn|Mn", p.names, ignore.case=F)


# all fewer
gridExtra::grid.arrange(p.plots[[10]], p.plots[[12]], p.plots[[24]], p.plots[[27]], ncol=2)
gridExtra::grid.arrange(p.plots[[45]], p.plots[[48]], p.plots[[64]], p.plots[[65]], ncol=2)


# representatives of each category
gridExtra::grid.arrange(n.plots[[3]], k.plots[[1]], g.plots[[5]], p.plots[[12]], ncol=2)

cairo_ps(filename = paste("cog_num.representative",format(Sys.time(),"%Y%m%d"),"eps",sep="."))
gridExtra::grid.arrange(n.plots[[3]], k.plots[[1]], g.plots[[5]], p.plots[[12]], ncol=2)
dev.off()


