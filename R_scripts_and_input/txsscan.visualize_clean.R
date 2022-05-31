#### txsscan secretion system summarize and visualize ####

library(ape)
library(phytools)
library(viridis)

## phylogeny 
# most comprehensive tree we have for paraburkholderia based on whole genomes
name.tree <- read.tree("../burk_naming/aaf_bpc_k43.tre")
plot(name.tree)

# prune this tree to only leave taxa used here
species <- c("B.agricolaris_BaQS159", "B.bonniea_BbQS859", "B.hayleyella_BhQS11", "P.caledonica_NBRC102488", "P.fungorum_ATCC.BAA-463", "P.megapolitana_LMG23650","P.phenazium_LMG2247" ,"P.phenoliruptrix_BR3459a", "P.phymatum_STM815", "P.phytofirmans_PsJn", "P.sartisoli_LMG24000", "P.sprentiae_WSM5005", "P.terrae_NBRC100964" ,"P.terricola_LMG20594", "P.xenovorans_LB400")
pruned.tree<-drop.tip(name.tree, name.tree$tip.label[-match(species, name.tree$tip.label)])
plot(pruned.tree)
nodelabels()

# write tree
write.tree(pruned.tree,file="comp_pruned.tre")
rm(name.tree)


## secretion system summary 

# reorganized and double checked files
df.txs <- read.table("txsscan_summary.txt", sep="\t", h=T)
rownames(df.txs) <- c("pagri","pbonn","phayl","pcale","pfung","pmega","pphem","pphex","pphym","pphyt","psart","pspre","ptera","ptere","pxeno")
df.txs$status <- factor(c("dicty","dicty","dicty","free-living","symbiotic","symbiotic","free-living","symbiotic","symbiotic","symbiotic","free-living","symbiotic","free-living","free-living","free-living"))

str(df.txs) # check the column types
summary(df.txs) # check max systems

# read in tree
pruned.tree <- read.newick(file="txsscan_pruned.tre")
pruned.tree$tip.label <- c("phayl","pbonn","psart","pcale","ptera","pphyt","pxeno","pagri","pfung","pphem","pspre","pphex","pphym","ptere","pmega")
temp <- ladderize(pruned.tree)
pruned.tree <- temp

# visualize 
ss.mat <- as.matrix(df.txs[,2:6])

# dot tree
par(mar=c(0,0,0,0))
dotTree(pruned.tree, ss.mat, labels=T, length=5, colors="black", mar=c(2,0,8,2))

# heat map 
vir_palette <- viridis(n = 5) # need at least one more than the max number you want to show
cols <- setNames(vir_palette, seq(0,4))

par(mar=c(5.1,4.1,4.1,2.1)) ## reset margins to default
phylo.heatmap(pruned.tree, ss.mat, standardize=F, split=c(1,2), col=vir_palette, no.margin=T, grid=T, mar=c(8,4,10,2))



cairo_ps(filename = paste("sec_sys.all",format(Sys.time(),"%Y%m%d"),"eps",sep="."))
par(mar=c(0,0,0,0))
dotTree(pruned.tree, ss.mat, labels=T, length=5, colors="black", mar=c(2,0,8,2))
dev.off()

