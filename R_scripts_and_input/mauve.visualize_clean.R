#### progressiveMauve whole genome alignment visualization ####

library(ggraph) #https://www.data-imaginist.com/2017/ggraph-introduction-layouts/
library(igraph)

## comparison to the 4 core control genomes 
# genomes are in the order: pagri, pbonn, phayl, pfung, pspre, ptere, pxeno

# read in file of lcb locations
core.lcb <- read.table("mauve_10cont.4core.lcb", h=T)
head(core.lcb)

# get lengths of lcb
core.lcb$len0 <- abs(core.lcb$seq0_leftend - core.lcb$seq0_rightend) 
core.lcb$len1 <- abs(core.lcb$seq1_leftend - core.lcb$seq1_rightend) 
core.lcb$len2 <- abs(core.lcb$seq2_leftend - core.lcb$seq2_rightend) 
core.lcb$len3 <- abs(core.lcb$seq3_leftend - core.lcb$seq3_rightend) 
core.lcb$len4 <- abs(core.lcb$seq4_leftend - core.lcb$seq4_rightend) 
core.lcb$len5 <- abs(core.lcb$seq5_leftend - core.lcb$seq5_rightend) 
core.lcb$len6 <- abs(core.lcb$seq6_leftend - core.lcb$seq6_rightend) 

nrow(core.lcb) # 153 lcb to begin with, with lcb lengths ranging from 9-283k
summary(core.lcb)

# filter lcb for those larger than median size so plot is easier to read
quantile(core.lcb$len0)
core.lcb.filt <- subset(core.lcb, len0 >=13112)

nrow(core.lcb.filt) #for reference, 94 at 3252.5, 63 at 9559

keep <- row.names(core.lcb.filt)

# read in file of lcb orientations
temp <- read.delim("mauve_10cont.4core.perm", sep=",", h=F)
lcb <- as.data.frame(t(temp[,-1]))
colnames(lcb) <- c("seq0", "seq1", "seq2", "seq3", "seq4", "seq5", "seq6")

# make vectors of lcb orientation
seq0 <- as.vector(lcb$seq0)
seq1 <- as.vector(lcb$seq1)
seq2 <- as.vector(lcb$seq2)
seq3 <- as.vector(lcb$seq3)
seq4 <- as.vector(lcb$seq4)
seq5 <- as.vector(lcb$seq5)
seq6 <- as.vector(lcb$seq6)

# keep lcb orientation over median size
agri <- seq0[abs(as.numeric(seq0)) %in% keep]
bonn <- seq1[abs(as.numeric(seq1)) %in% keep]
hayl <- seq2[abs(as.numeric(seq2)) %in% keep]
fung <- seq3[abs(as.numeric(seq3)) %in% keep]
spre <- seq4[abs(as.numeric(seq4)) %in% keep]
tere <- seq5[abs(as.numeric(seq5)) %in% keep]
xeno <- seq6[abs(as.numeric(seq6)) %in% keep]


# define utility functions dataframe for visualization
df.hive <- core.lcb.filt[,c(1,3,5,7,9,11,13,15:21)]
colnames(df.hive) <- c("agri.left","bonn.left","hayl.left","fung.left","spre.left","tere.left","xeno.left","agri.len","bonn.len","hayl.len","fung.len","spre.len","tere.len","xeno.len")

foo <- function(x) {
  y = paste0(substitute(x),".left")
  abs(df.hive[[y]])
  #return(y)
}
bar <- function(x) {
  y = paste0(substitute(x),".len")
  df.hive[[y]]
  #return(y)
}

# clean up environment
rm(temp,seq0,seq1,seq2,seq3,seq4,seq5,seq6)
rm(core.lcb, core.lcb.filt, lcb)


## start of visualization ##

# generic dataframe that describes how nodes are linked to each other
linkfile <- data.frame(id1 = as.integer(c(seq(1,length(keep)),seq(2*length(keep)+1,3*length(keep)),seq(length(keep)+1,2*length(keep)))),
                       id2 = as.integer(c(seq(length(keep)+1,2*length(keep)),seq(1,length(keep)),seq(2*length(keep)+1,3*length(keep))))
)

# 1. fung(3) vs. agri(0) vs. xeno(6) 

# distinguish lcb in the same direction between each pair (indexing here)
diraf <- as.character(sort(abs(intersect(agri, fung))))
dirax <- as.character(sort(abs(intersect(agri, xeno))))
dirfx <- as.character(sort(abs(intersect(fung, xeno))))

# define igraph object
hive.ag.1 <- graph_from_data_frame(linkfile)
# adding alphabet labels to control order, clockwise from 0/12 o'clock
V(hive.ag.1)$axis <- c(rep("A_agri",length(keep)), rep("C_xeno",length(keep)), rep("B_fung",length(keep)))
V(hive.ag.1)$edge <- c(ifelse(keep %in% diraf, "af.f","af.r"), ifelse(keep %in% dirax, "ax.f","ax.r"), ifelse(keep %in% dirfx, "fx.f","fx.r"))
V(hive.ag.1)$radius <- c(foo(agri), foo(xeno), foo(fung))
V(hive.ag.1)$width <- c(bar(agri), bar(xeno), bar(fung))

# can adjust colors, these are from paletton and colorgorical
# "#47d3d3","#d34747"
hivegraph.agri.1 <- 
ggraph(hive.ag.1, 'hive', axis = V(hive.ag.1)$axis, sort.by = V(hive.ag.1)$radius, normalize=F, use.numeric=T) + geom_edge_hive(aes(colour = V(hive.ag.1)$edge, alpha=0.1, width= V(hive.ag.1)$width)) + geom_axis_hive(size = 2,colour=c("white","white","#47d3d3")) + scale_edge_color_manual(values=rep(c("#58b5e1", "#8e2283"),3)) + theme_void() + coord_fixed() + scale_edge_alpha_continuous(range = (0.3)) + theme(legend.position="none")


# spre(4) vs. agri(0) vs. tere(5) 
diras <- as.character(sort(abs(intersect(agri, spre))))
dirat <- as.character(sort(abs(intersect(agri, tere))))
dirst <- as.character(sort(abs(intersect(spre, tere))))

hive.ag.2 <- graph_from_data_frame(linkfile)
V(hive.ag.2)$axis <- c(rep("A_agri",length(keep)), rep("C_tere",length(keep)), rep("B_spre",length(keep)))
V(hive.ag.2)$edge <- c(ifelse(keep %in% diras, "as.f","as.r"), ifelse(keep %in% dirat, "at.f","at.r"), ifelse(keep %in% dirst, "st.f","st.r"))
V(hive.ag.2)$radius <- c(foo(agri), foo(tere), foo(spre))
V(hive.ag.2)$width <- c(bar(agri), bar(tere), bar(spre))

hivegraph.agri.2 <- 
ggraph(hive.ag.2, 'hive', axis = V(hive.ag.2)$axis, sort.by = V(hive.ag.2)$radius, normalize=F, use.numeric=T) + geom_edge_hive(aes(colour = V(hive.ag.2)$edge, alpha=0.1, width= V(hive.ag.2)$width)) + geom_axis_hive(size = 2,colour=c("white","white","#47d3d3")) + scale_edge_color_manual(values=rep(c("#58b5e1", "#8e2283"),3)) + theme_void() + coord_fixed() + scale_edge_alpha_continuous(range = (0.3)) + theme(legend.position="none")


# fung(3) vs. bonn(1) vs. xeno(6) 
dirbf <- as.character(sort(abs(intersect(bonn, fung))))
dirbx <- as.character(sort(abs(intersect(bonn, xeno))))
dirfx <- as.character(sort(abs(intersect(fung, xeno))))

hive.bo.1 <- graph_from_data_frame(linkfile)
V(hive.bo.1)$axis <- c(rep("A_bonn",length(keep)), rep("C_xeno",length(keep)), rep("B_fung",length(keep)))
V(hive.bo.1)$edge <- c(ifelse(keep %in% dirbf, "bf.f","bf.r"), ifelse(keep %in% dirbx, "bx.f","bx.r"), ifelse(keep %in% dirfx, "fx.f","fx.r"))
V(hive.bo.1)$radius <- c(foo(bonn), foo(xeno), foo(fung))
V(hive.bo.1)$width <- c(bar(bonn), bar(xeno), bar(fung))

# "#d34747", "#2b7272", "#99c542"
hivegraph.bonn.1 <- 
ggraph(hive.bo.1, 'hive', axis = V(hive.bo.1)$axis, sort.by = V(hive.bo.1)$radius, normalize=F, use.numeric=T) + geom_edge_hive(aes(colour = V(hive.bo.1)$edge, alpha=0.1, width= V(hive.bo.1)$width)) + geom_axis_hive(size = 2,colour=c("white","white","#d34747")) + scale_edge_color_manual(values=rep(c("#58b5e1", "#8e2283"),3)) + theme_void() + coord_fixed() + scale_edge_alpha_continuous(range = (0.3)) + theme(legend.position="none")


# spre(4) vs. bonn(1) vs. tere(5) 
dirbs <- as.character(sort(abs(intersect(bonn, spre))))
dirbt <- as.character(sort(abs(intersect(bonn, tere))))
dirst <- as.character(sort(abs(intersect(spre, tere))))

hive.bo.2 <- graph_from_data_frame(linkfile)
V(hive.bo.2)$axis <- c(rep("A_bonn",length(keep)), rep("C_tere",length(keep)), rep("B_spre",length(keep)))
V(hive.bo.2)$edge <- c(ifelse(keep %in% dirbs, "bs.f","bs.r"), ifelse(keep %in% dirbt, "bt.f","bt.r"), ifelse(keep %in% dirst, "st.f","st.r"))
V(hive.bo.2)$radius <- c(foo(bonn), foo(tere), foo(spre))
V(hive.bo.2)$width <- c(bar(bonn), bar(tere), bar(spre))

# "#d34747", "#2b7272", "#99c542"
hivegraph.bonn.2 <- 
ggraph(hive.bo.2, 'hive', axis = V(hive.bo.2)$axis, sort.by = V(hive.bo.2)$radius, normalize=F, use.numeric=T) + geom_edge_hive(aes(colour = V(hive.bo.2)$edge, alpha=0.1, width= V(hive.bo.2)$width)) + geom_axis_hive(size = 2,colour=c("white","white","#d34747")) + scale_edge_color_manual(values=rep(c("#58b5e1", "#8e2283"),3)) + theme_void() + coord_fixed() + scale_edge_alpha_continuous(range = (0.3)) + theme(legend.position="none")


# fung(3) vs. hayl(2) vs. xeno(6) 
dirhf <- as.character(sort(abs(intersect(hayl, fung))))
dirhx <- as.character(sort(abs(intersect(hayl, xeno))))
dirfx <- as.character(sort(abs(intersect(fung, xeno))))

hive.ha.1 <- graph_from_data_frame(linkfile)
V(hive.ha.1)$axis <- c(rep("A_hayl",length(keep)), rep("C_xeno",length(keep)), rep("B_fung",length(keep)))
V(hive.ha.1)$edge <- c(ifelse(keep %in% dirhf, "hf.f","hf.r"), ifelse(keep %in% dirhx, "hx.f","hx.r"), ifelse(keep %in% dirfx, "fx.f","fx.r"))
V(hive.ha.1)$radius <- c(foo(hayl), foo(xeno), foo(fung))
V(hive.ha.1)$width <- c(bar(hayl), bar(xeno), bar(fung))

# "#d34747", "#2b7272", "#99c542"
hivegraph.hayl.1 <- 
ggraph(hive.ha.1, 'hive', axis = V(hive.ha.1)$axis, sort.by = V(hive.ha.1)$radius, normalize=F, use.numeric=T) + geom_edge_hive(aes(colour = V(hive.ha.1)$edge, alpha=0.1, width= V(hive.ha.1)$width)) + geom_axis_hive(size = 2,colour=c("white","white","#d34747")) + scale_edge_color_manual(values=rep(c("#58b5e1", "#8e2283"),3)) + theme_void() + coord_fixed() + scale_edge_alpha_continuous(range = (0.3)) + theme(legend.position="none")


# spre(4) vs. hayl(2) vs. tere(5) 
dirhs <- as.character(sort(abs(intersect(hayl, spre))))
dirht <- as.character(sort(abs(intersect(hayl, tere))))
dirst <- as.character(sort(abs(intersect(spre, tere))))

hive.ha.2 <- graph_from_data_frame(linkfile)
V(hive.ha.2)$axis <- c(rep("A_hayl",length(keep)), rep("C_tere",length(keep)), rep("B_spre",length(keep)))
V(hive.ha.2)$edge <- c(ifelse(keep %in% dirbs, "hs.f","hs.r"), ifelse(keep %in% dirbt, "ht.f","ht.r"), ifelse(keep %in% dirst, "st.f","st.r"))
V(hive.ha.2)$radius <- c(foo(hayl), foo(tere), foo(spre))
V(hive.ha.2)$width <- c(bar(hayl), bar(tere), bar(spre))

# "#d34747", "#2b7272", "#99c542"
hivegraph.hayl.2 <- 
ggraph(hive.ha.2, 'hive', axis = V(hive.ha.2)$axis, sort.by = V(hive.ha.2)$radius, normalize=F, use.numeric=T) + geom_edge_hive(aes(colour = V(hive.ha.2)$edge, alpha=0.1, width= V(hive.ha.2)$width)) + geom_axis_hive(size = 2, colour=c("white","white","#d34747"), alpha=1) + scale_edge_color_manual(values=rep(c("#58b5e1", "#8e2283"),3)) + theme_void() + coord_fixed() + scale_edge_alpha_continuous(range = (0.3)) + theme(legend.position="none")


# agri vs. bonn vs. hayl
dirab <- as.character(sort(abs(intersect(agri, bonn))))
dirah <- as.character(sort(abs(intersect(agri, hayl))))
dirbh <- as.character(sort(abs(intersect(hayl, bonn))))

hive.dicty.1 <- graph_from_data_frame(linkfile)
V(hive.dicty.1)$axis <- c(rep("A_agri",length(keep)), rep("C_hayl",length(keep)), rep("B_bonn",length(keep)))
V(hive.dicty.1)$edge <- c(ifelse(keep %in% dirab, "ab.f","ab.r"), ifelse(keep %in% dirah, "ah.f","ah.r"), ifelse(keep %in% dirbh, "bh.f","bh.r"))
V(hive.dicty.1)$radius <- c(foo(agri), foo(hayl), foo(bonn))
V(hive.dicty.1)$width <- c(bar(agri), bar(hayl), bar(bonn))

hivegraph.dicty.1 <- ggraph(hive.dicty.1, 'hive', axis = V(hive.dicty.1)$axis, sort.by = V(hive.dicty.1)$radius, normalize=F, use.numeric=T) + geom_edge_hive(aes(colour = V(hive.dicty.1)$edge, alpha=0.1, width= V(hive.dicty.1)$width)) + geom_axis_hive(size = 2,colour=c("#d34747","#d34747","#47d3d3")) + theme_void() + coord_fixed() + scale_edge_alpha_continuous(range = (0.3)) + scale_edge_color_manual(values=rep(c("#58b5e1", "#8e2283"),3)) + theme(legend.position="none")


## dicty only lcb comparison 

dicty.lcb <- read.table("mauve_10cont.dicty.lcb", h=T)
head(dicty.lcb)

dicty.lcb$len0 <- abs(dicty.lcb$seq0_leftend - dicty.lcb$seq0_rightend) 
dicty.lcb$len1 <- abs(dicty.lcb$seq1_leftend - dicty.lcb$seq1_rightend) 
dicty.lcb$len2 <- abs(dicty.lcb$seq2_leftend - dicty.lcb$seq2_rightend) 

nrow(dicty.lcb) # 180 lcb to begin with
summary(dicty.lcb$len0)

# no filtering in this version
keep.d <- row.names(dicty.lcb)

temp <- read.delim("mauve_10cont.dicty.perm", sep=",", h=F)
lcb.d <- as.data.frame(t(temp[,-1]))
colnames(lcb.d) <- c("seq0", "seq1", "seq2")
seq0 <- as.vector(lcb.d$seq0)
seq1 <- as.vector(lcb.d$seq1)
seq2 <- as.vector(lcb.d$seq2)

agri.d <- seq0[abs(as.numeric(seq0)) %in% keep.d]
bonn.d <- seq1[abs(as.numeric(seq1)) %in% keep.d]
hayl.d <- seq2[abs(as.numeric(seq2)) %in% keep.d]

# define dataframe for visualization
df.hive.d <- dicty.lcb[,c(1,3,5,7:9)]
colnames(df.hive.d) <- c("agri.left","bonn.left","hayl.left","agri.len","bonn.len","hayl.len")
df.hive.d$agri.1 <- ifelse(abs(df.hive.d$agri.left) <= 4816966, 1, 2)
df.hive.d$bonn.1 <- ifelse(abs(df.hive.d$bonn.left) <= 3175376, 1, 2)
df.hive.d$hayl.1 <- ifelse(abs(df.hive.d$hayl.left) <= 3295139, 1, 2)

foo.d <- function(x) {
  y = paste0(substitute(x),".left")
  abs(df.hive.d[[y]])
  #return(y)
}
bar.d <- function(x) {
  y = paste0(substitute(x),".len")
  df.hive.d[[y]]
  #return(y)
}
rm(temp,seq0,seq1,seq2)
rm(dicty.lcb, lcb.d)

# these are the lcb in the same direction between each pair, indexing here so we can distinguish between these and inverted ones
dirab.d <- as.character(sort(abs(intersect(agri.d, bonn.d))))
dirah.d <- as.character(sort(abs(intersect(agri.d, hayl.d))))
dirbh.d <- as.character(sort(abs(intersect(hayl.d, bonn.d))))

linkfile.d <- data.frame(id1 = as.integer(c(seq(1,length(keep.d)),seq(2*length(keep.d)+1,3*length(keep.d)),seq(length(keep.d)+1,2*length(keep.d)))),
                            id2 = as.integer(c(seq(length(keep.d)+1,2*length(keep.d)),seq(1,length(keep.d)),seq(2*length(keep.d)+1,3*length(keep.d))))
)

hive.dicty.2 <- graph_from_data_frame(linkfile.d)
V(hive.dicty.2)$axis <- c(rep("A_agri",length(keep.d)), rep("C_hayl",length(keep.d)), rep("B_bonn",length(keep.d)))
V(hive.dicty.2)$edge <- c(ifelse(keep.d %in% dirab.d, "ab.f","ab.r"), ifelse(keep.d %in% dirah.d, "ah.f","ah.r"), ifelse(keep.d %in% dirbh.d, "bh.f","bh.r"))
V(hive.dicty.2)$radius <- c(foo.d(agri), foo.d(hayl), foo.d(bonn))
V(hive.dicty.2)$width <- c(bar.d(agri), bar.d(hayl), bar.d(bonn))

synt.abh.d <- c(paste0(df.hive.d$agri.1, df.hive.d$bonn.1), paste0(df.hive.d$hayl.1, df.hive.d$agri.1), paste0(df.hive.d$bonn.1, df.hive.d$hayl.1))
V(hive.dicty.2)$chr <- synt.abh.d


hivegraph.dicty.2 <- 
ggraph(hive.dicty.2, 'hive', axis = V(hive.dicty.2)$axis, sort.by = V(hive.dicty.2)$radius, normalize=F, use.numeric=T) + geom_edge_hive(aes(colour = V(hive.dicty.2)$edge, alpha=0.1, width= V(hive.dicty.2)$width)) + geom_axis_hive(size = 2,colour=c("#d34747","#d34747","#47d3d3")) + theme_void() + coord_fixed() + scale_edge_alpha_continuous(range = (0.3)) + scale_edge_color_manual(values=rep(c("#58b5e1", "#8e2283"),3)) + theme(legend.position="none")

ggraph(hive.dicty.2, 'hive', axis = V(hive.dicty.2)$axis, sort.by = V(hive.dicty.2)$radius, normalize=F, use.numeric=T) + geom_edge_hive(aes(colour = V(hive.dicty.2)$edge, alpha=0.1, width= V(hive.dicty.2)$width)) + geom_axis_hive(size = 2,colour=c("#d34747","#d34747","#d34747","#d34747","#d34747","#d34747","#d34747","#d34747","#47d3d3","#47d3d3","#47d3d3","#47d3d3")) + theme_void() + coord_fixed() + scale_edge_alpha_continuous(range = (0.3)) + scale_edge_color_manual(values=rep(c("#58b5e1", "#8e2283"),3)) + facet_edges(~V(hive.dicty.2)$chr) + theme(legend.position="none")


## save figures in a format that can be edited in illustrator 
# open in illustrator to fix labels, save as eps, than re-export from preview for best quality 
cairo_ps(filename = paste("hive.agri_control.1",format(Sys.time(),"%Y%m%d"),"eps",sep="."), width=7.5)
hivegraph.agri.1
dev.off()

cairo_ps(filename = paste("hive.agri_control.2",format(Sys.time(),"%Y%m%d"),"eps",sep="."), width=5.6)
hivegraph.agri.2
dev.off()

cairo_ps(filename = paste("hive.bonn_control.1",format(Sys.time(),"%Y%m%d"),"eps",sep="."), width=7.5)
hivegraph.bonn.1
dev.off()

cairo_ps(filename = paste("hive.bonn_control.2",format(Sys.time(),"%Y%m%d"),"eps",sep="."), width=5.6)
hivegraph.bonn.2
dev.off()

cairo_ps(filename = paste("hive.hayl_control.1",format(Sys.time(),"%Y%m%d"),"eps",sep="."), width=7.5)
hivegraph.hayl.1
dev.off()

cairo_ps(filename = paste("hive.hayl_control.2",format(Sys.time(),"%Y%m%d"),"eps",sep="."), width=5.6)
hivegraph.hayl.2
dev.off()

cairo_ps(filename = paste("hive.dicty_burk",format(Sys.time(),"%Y%m%d"),"eps",sep="."), height=8)
hivegraph.dicty.1
dev.off()

cairo_ps(filename = paste("hive.dicty_burk.2",format(Sys.time(),"%Y%m%d"),"eps",sep="."), height=8)
hivegraph.dicty.2
dev.off()
