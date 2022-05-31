#### 2 steps analysis for codeml model selection and parameter visualization ####
## to be used in conjunction with shell scripts summarize_codeml_1_logl.sh and summarize_codeml_2_parameters.sh

library(ggplot2)
library(reshape2)

## step 1 - after running codeml models, import logLikelihood scores and compare to select best fitting model

# hypotheses tested in codeml
# H0: all rates same across tree
# H1a: Different rates of evolution when symbiotic # symbiotic
# H1b: Different rates of evolution in dicty # dicty
# H1c: Different rates of evolution when genomes reduced # reduced

# read in data and rename columns
codeml <- read.table("codeml_likelihood.models.txt") 
colnames(codeml) <- c("gene","h0", "symbiotic", "dicty", "reduced")

# calculate AIC from logLikelihood scores
codeml$aic_h0 <- -2*codeml$h0
codeml$aic_sym <- 4-2*codeml$symbiotic
codeml$aic_dic <- 4-2*codeml$dicty
codeml$aic_red <- 4-2*codeml$reduced

# find minimum AIC value for each gene
codeml$which <- apply(codeml[,6:9],1,which.min)
codeml$aic_min <- apply(codeml[,6:9],1,min)
codeml$diff <- abs(codeml$aic_h0 - codeml$aic_min)
codeml$result <- c("h0")
codeml$result[codeml$which == 2] <- c("symbiotic") 
codeml$result[codeml$which == 3] <- c("dicty") 
codeml$result[codeml$which == 4] <- c("reduced") 

table(codeml$result)
#dicty        h0   reduced symbiotic 
#157       898       424       194  

table(subset(codeml, diff>0 & diff<1)[13])
summary(subset(codeml, diff!=0)[12])

# correct best model for ones that are at least 1 AIC different from null model h0
codeml$result[codeml$diff < 1] <- c("h0")
table(codeml$result)
#dicty        h0   reduced symbiotic 
#137      1001       372       163 

# write out the names of genes for each model to retrieve appropriate parameter estimates
write.table(codeml$gene[codeml$result == "h0"], "codeml_genes.h0.txt", quote=F, col.names = F, row.names = F)
write.table(codeml$gene[codeml$result == "symbiotic"], "codeml_genes.sym.txt", quote=F, col.names = F, row.names = F)
write.table(codeml$gene[codeml$result == "dicty"], "codem_genesl.dic.txt", quote=F, col.names = F, row.names = F)
write.table(codeml$gene[codeml$result == "reduced"], "codeml_genes.red.txt", quote=F, col.names = F, row.names = F)


## step 2 - after retrieving parameter estimates for best model, summarize and visualize

# read in data
omega <- read.table("codeml_omega.diff1.txt")
colnames(omega) <- c("gene","model", "kappa", "w0", "w1","lnL")
omega$w1 <- ifelse(is.na(omega$w1), omega$w0, omega$w1)

# inspect omega estimates between models
ggplot(omega, aes(x=w0, y=w1, colour=model)) + geom_point() 

# compare omegas for each model
omega.long <- melt(omega[,c(1:2,4:5)], id.vars = 1:2, variable.name = "class", value.name = "estimate")

ggplot(omega.long, aes(x=class, y=estimate, group=class, fill=class)) + geom_boxplot() + facet_grid(.~model)
alt.lab <- c("Null", "Symbiotic", "Dicty", "Reduced")
names(alt.lab) <- c("H0", "H1a", "H1b", "H1c")

omega.sub <- subset(omega.long, !(model=="H0" & class=="w1"))
ggplot(omega.sub, aes(x=class, y=estimate, group=class, fill=class)) + geom_boxplot() + facet_grid(.~model, labeller=labeller(model=alt.lab)) + scale_fill_manual(values=c("#47d3d3","#d34747")) + theme_bw()

# compare background and alternate omega estimates within each model
wilcox.test(estimate ~ class, data = subset(omega.sub, model=="H1a"), paired = TRUE)
wilcox.test(estimate ~ class, data = subset(omega.sub, model=="H1b"), paired = TRUE)
wilcox.test(estimate ~ class, data = subset(omega.sub, model=="H1c"), paired = TRUE)


cairo_ps(filename = paste("codeml.models",format(Sys.time(),"%Y%m%d"),"eps",sep="."), width=7.5)
ggplot(omega.sub, aes(x=class, y=estimate, group=class, fill=class)) + geom_boxplot() + facet_grid(.~model, labeller=labeller(model=alt.lab)) + scale_fill_manual(values=c("#47d3d3","#d34747")) + theme_bw()
dev.off()