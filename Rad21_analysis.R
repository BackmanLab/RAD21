rm(list = ls())
gc()

# Notes
# Nearest neighbor analysis not used in paper analysis directly
# Bed file names from Encode changed to more intuitive names
# Proseq not used in final version of paper

#########################################################################################
# Paper functions
#########################################################################################

## nearest neighbor function
getNN <- function(dir.bed, bed, interval, cond, anchors){
  
  ##
  ##
  ##
  ##
  
  require(plyranges)
  require(dplyr)
  
  ##-------------------------------------------------------## BED file
  
  bed <- as.data.frame(read.table(file.path(dir.bed,bed), sep= "\t", header= T))
  colnames(bed) <- c("chr", "start", "end", "na", "score", "na", "signalValue", "pval", "qval", "peak")
  
  bed <- bed  %>% transform( seqnames= paste0(chr), start = start , end = end)  %>% as_granges()
  interval <- interval %>% mutate(size = abs(V2-V6))
  
  ##-------------------------------------------------------## Loop anchors overlap
  if ( anchors == T ){
    
    loops.5p <- interval %>% transform( seqnames= paste0(V1), start = V2, end = V3)  %>% as_granges()
    loops.3p <- interval  %>% transform( seqnames= paste0(V1), start = V5, end = V6)  %>% as_granges()
    
    nn.3p <- data.frame(plyranges::join_nearest(loops.3p, bed,  distance = T))
    nn.5p <- data.frame(plyranges::join_nearest(loops.5p, bed,  distance = T))
    
    df.out <- rbind(nn.3p, nn.5p)
    #df.out <- df.out  %>% dplyr::mutate(dups = paste0(V1, "_", V22, "_",V23)) %>% dplyr::distinct(dups, .keep_all = T) ## Filter dups
    df.out$cond <- rep(cond, nrow(df.out))
    df.out$cond <- gsub(".bed", "", df.out$cond)
    df.out$cond <- gsub("_", " ", df.out$cond)
    
    return(df.out)
    
  } else {
    ##-------------------------------------------------------## Loop centroid overlap
    loops <- interval %>% transform( seqnames= paste0(V1), start = V2, end = V6)  %>% as_granges()
    
    df.out <- data.frame(plyranges::join_nearest(loops, bed,  distance = T))
    #df.out <- df.out  %>% dplyr::mutate(dups = paste0(V1, "_", V22, "_",V23)) %>% dplyr::distinct(dups, .keep_all = T) ## Filter dups
    df.out$cond <- rep(cond, nrow(df.out))
    df.out$cond <- gsub(".bed", "", df.out$cond)
    df.out$cond <- gsub("_", " ", df.out$cond)
    
    return(df.out)
    
  }
  
}

## overlap Analysis functions 
getOverlap <- function(dir.bed, bed, interval, cond, anchors){
  
  ##
  ##
  ##
  ##
  
  require(plyranges)
  require(dplyr)
  
  ##-------------------------------------------------------## BED file
  
  bed <- as.data.frame(read.table(file.path(dir.bed,bed), sep= "\t", header= T))
  
  colnames(bed) <- c("chr", "start", "end", "na", "score", "na", "signalValue", "pval", "qval", "peak") ## bed file
  #colnames(bed) <-c("chr", "start", "end", "GENEID", "TXNAME", "TXSTRAND", "signal") ## Pro seq
  
  bed <- bed  %>% transform( seqnames= paste0(chr), start = start , end = end)  %>% as_granges()
  interval <- interval %>% mutate(size = abs(V2-V6))
  
  ##-------------------------------------------------------## Loop anchors overlap
  if ( anchors == T ){
    
    loops.5p <- interval %>% transform( seqnames= paste0(V1), start = V2, end = V3)  %>% as_granges()
    loops.3p <- interval  %>% transform( seqnames= paste0(V1), start = V5, end = V6)  %>% as_granges()
    
    nn.3p <- data.frame(plyranges::join_overlap_inner(loops.3p, bed))
    nn.5p <- data.frame(plyranges::join_overlap_inner(loops.5p, bed))
    
    df.out <- rbind(nn.3p, nn.5p)
    
    df.out<- df.out  %>% dplyr::mutate(coord = paste0(V1, "_", V2, "_",V6, "_", cond))
    #df.out$coord <- gsub(".bed", "", df.out$coord) ## For bed files
    df.out$coord <- gsub(".tsv", "", df.out$coord) ## For PRO seq files
    
    return(df.out)
    
  } else {
    ##-------------------------------------------------------## Loop centroid overlap
    loops <- interval %>% transform( seqnames= paste0(V1), start = V2, end = V6)  %>% as_granges()
    
    df.out <- data.frame(join_overlap_inner(loops, bed))
    
    df.out<- df.out  %>% dplyr::mutate(coord = paste0(V1, "_", V2, "_",V6, "_", cond))
    df.out$coord <- gsub(".bed", "", df.out$coord) ## For bed files
    #df.out$coord <- gsub(".tsv", "", df.out$coord) ## For PRO seq files
    
    return(df.out)
    
  }
  
}

#########################################################################################
# Directories and packages
#########################################################################################

library(tidyr)
library(ggplot2)
library(dplyr)
library(forcats)
library(ggridges)
library(ggdist)
library(gghalves)
library(colorspace)
library(plyr)
library(stats)
library(viridis)

## Change to your root directory
root.dir <- "/Volumes/external hd/IBiS/Backman_Lab/projects/Rad21_paper/resubmission"

## dirs
bed.dir <- file.path(root.dir,"rad21_chip")
loop.dir <- file.path(root.dir,"micro_c")

#########################################################################################
# Nearest Neighbor Analysis of ChIP marks for rad21
#########################################################################################

## Change to your root directory
root.dir <- "/Volumes/external hd/IBiS/Backman_Lab/projects/Rad21_paper/resubmission"

## dirs
bed.dir <- file.path(root.dir,"rad21_chip")
loop.dir <- file.path(root.dir,"micro_c")

loops <- read.table(file.path(loop.dir,"ctrl_loops.bedpe"), sep= "\t", header= F)
tads <- read.table(file.path(loop.dir,"ctrl_tads.bedpe"), sep= "\t", header= F)

## Loop through tads and loops to get nearest distance and signal
df.lps <- data.frame()
conds <- c("h3k9me3_ctrl.bed", "h3k9me3_6hrsAux.bed", "h3k4me3_ctrl.bed", "h3k4me3_6hrsAux.bed","h3k27ac_ctrl.bed", "h3k27ac_6hrsAux.bed", "h3k27me3_6hrsAux.bed", "h3k27me3_ctrl.bed")
for (i in 1:length(conds)){
  
  cond <- conds[i]
  print(cond)
  
  x <- getNN(dir.bed = bed.dir, bed = cond,interval=loops, cond = paste0(cond, "_loops"),anchors = F)
  df.lps <- rbind(df.lps, x)
}

df.tads <- data.frame()
for (i in 1:length(conds)){
  
  cond <- conds[i]
  print(cond)
  
  y <- getNN(dir.bed = bed.dir, bed = cond,interval=tads, cond = paste0(cond, "_tads"),anchors = F)
  df.tads <- rbind(df.tads, y)
}

## Prepare data for plotting
## equalize data frames
n.lps <- data.frame(df.lps[,c(1:4, 39:49)])
n.tds <- data.frame(df.tads[,c(1:4, 22:32)])

## Wranlge data
n.df <- rbind(n.lps,n.tds)
n.df$cond <- paste0(n.df$cond," ", ntile(n.df$size, as.numeric(7)))
n.df <- tidyr::separate_wider_delim(n.df,cond,delim = " ",names = c("mark", "status","top", "bin"), cols_remove=F)

## Get mean bin size for loops
size.lab <- n.df %>% dplyr::group_by(bin) %>% dplyr::reframe(mean.size = mean(size)) 

## filter on criteria
filt.lst <- c("h3k9me3","h3k4me3","h3k27ac", "h3k27me3")
plot <- n.df[n.df$mark == filt.lst[1],]

## Plotting
ggplot(plot, aes(x = signalValue, y = cond, fill = stat(x))) +
  geom_density_ridges_gradient(scale = 3, rel_min_height = 0.01) +
  scale_fill_viridis_c(name = "Signal", option = "C") +
  labs(title = 'h3k9me3')

ggplot(plot, aes(x = signalValue, y=cond, fill = factor(stat(quantile)))) +
  stat_density_ridges(
    geom = "density_ridges_gradient", calc_ecdf = TRUE,
    quantiles = 4, quantile_lines = TRUE
  ) + scale_x_continuous(limits = c(0,8), breaks = seq(1,8,1))+
  scale_fill_viridis_d(name = "Quartiles") + theme_classic()


#########################################################################################
# Overlap Analysis of ChIP marks for rad21
#########################################################################################

## Loop/TAD files
loops <- read.table(file.path(loop.dir,"ctrl_loops.bedpe"), sep= "\t", header= F)
tads <- read.table(file.path(loop.dir,"ctrl_tads.bedpe"), sep= "\t", header= F)

## Loop through tads and loops to get nearest distance and signal
df.lps <- data.frame()
conds <- c("h3k9me3_ctrl.bed", "h3k9me3_6hrsAux.bed", "h3k4me3_ctrl.bed", "h3k4me3_6hrsAux.bed","h3k27ac_ctrl.bed", "h3k27ac_6hrsAux.bed", "h3k27me3_6hrsAux.bed", "h3k27me3_ctrl.bed")
for (i in 1:length(conds)){
  
  cond <- conds[i]
  print(cond)
  
  x <- getOverlap(dir.bed = bed.dir, bed = cond,interval=loops, cond = paste0(cond, "_loops"),anchors = F)
  df.lps <- rbind(df.lps, x)
}

## Bring in TADs
df.tads <- data.frame()
for (i in 1:length(conds)){
  
  cond <- conds[i]
  print(cond)
  
  y <- getOverlap(dir.bed = bed.dir, bed = cond,interval=tads, cond = paste0(cond, "_tads"),anchors = F)
  df.tads <- rbind(df.tads, y)
}

## Reframe data
n.tds <- df.tads %>% dplyr::group_by(coord) %>% dplyr::reframe(signal = sum(signalValue))
n.lps <- df.lps %>% dplyr::group_by(coord) %>% dplyr::reframe(signal = sum(signalValue))
n.df <- rbind(n.tds,n.lps)

## Prepare data for plotting
n.df <- tidyr::separate_wider_delim(n.df,coord,delim = "_",names = c("chrom","start","end", "mark", "status","top"), cols_remove=F)
n.df <- n.df %>% dplyr::mutate(size = abs(as.numeric(start)-as.numeric(end)))
n.df$bin <- dplyr::ntile(n.df$size, as.numeric(7))
n.df$cond <- paste0(n.df$mark," ", n.df$status, " ", n.df$top, " ",n.df$bin)

##---------------------------------------- Generating SI and main text plots -##

## Filter on criteria
filt.lst <- c("h3k9me3","h3k4me3","h3k27ac", "h3k27me3")
data <- n.df[n.df$mark == filt.lst[3],]
data <- data %>% dplyr::filter(bin %in% c(1,4,5)) ## for maintext plots

## Plot parameters
lim = 15 ## 
title = c("H3K27ac Coverage", "H3K9me3 Coverage")

## Get mean bin size for loops and generate labels
size.lab <- data %>% group_by(bin) %>% reframe(mean.size = mean(size)) 
lab <- paste0(format(as.integer(rep(round_any(size.lab$mean.size,1000),4)),big.mark=",",scientific=F, trim = T), " bp") ## lab for facet plot
lab <- paste0(sapply(strsplit(as.character(lab), ","), `[`, 1), " KBP")
       
## Main text plot 1
# k9.micro.v1 # k27.micro.v1
p1 <- ggplot(data, aes(x = log2(signal), y = as.factor(bin), fill = stat(x))) +
  geom_density_ridges_gradient(scale = 2, rel_min_height = 0.01) +
  scale_fill_viridis_c(name = "Signal", option = "C") +
  labs(title = title[1]) + scale_x_continuous(limits = c(0,lim), breaks = seq(0,lim,5))+
  scale_y_discrete( labels = lab) + ylab("") +
  theme(axis.text.x = element_text(angle = 0,  hjust=1), legend.position="", plot.title = element_text(size=16), text = element_text(size=16, family="Arial"))
p2 <- p1 + facet_grid(fct_rev(status) ~ top,axes = "all", axis.labels = "all_x")
p2 + theme_minimal() + theme(axis.text.x = element_text(angle = 0,  hjust=1), legend.position="", plot.title = element_text(size=16), text = element_text(size=16, family="Arial"))

## Main text plot 2
# k9.micro.v2 # k27.micro.v2
p1 <- ggplot(data, aes(x = log2(signal), y=as.factor(bin), fill = factor(stat(quantile)))) +
  stat_density_ridges(
    geom = "density_ridges_gradient", calc_ecdf = TRUE,
    quantiles = 4, quantile_lines = TRUE
  ) + scale_fill_viridis_d(name = "Quartiles") +labs(title = title) +
  labs(title = title[1]) + scale_x_continuous(limits = c(0,lim), breaks = seq(0,lim,5))+
  scale_y_discrete( labels = lab) + ylab("")
p2 <- p1 + facet_grid(fct_rev(status) ~ top,axes = "all", axis.labels = "all_x")
p2 + theme_classic() + theme(axis.text.x = element_text(angle = 0,  hjust=1), legend.position="", plot.title = element_text(size=18), text = element_text(size=18, family="Arial"))

#########################################################################################
# Overlap Analysis of ChIP marks for rad21 : Get chromosome totals
#########################################################################################

## Loop through beds 
df.beds <- data.frame()
conds <- c("h3k9me3_ctrl.bed", "h3k9me3_6hrsAux.bed","h3k27ac_ctrl.bed", "h3k27ac_6hrsAux.bed")
for (i in 1:length(conds)){
  
  cond <- conds[i]
  print(cond)
  
  x <- as.data.frame(read.table(file.path(bed.dir,cond), sep= "\t", header= T))
  colnames(x) <- c("chr", "start", "end", "na", "score", "na", "signalValue", "pval", "qval", "peak") ## bed file
  
  x$cond <- paste0(cond)
  x$cond <- gsub(".bed", "", x$cond)
  df.beds <- rbind(df.beds, x)
}

df.beds$signalValue = rep(1, nrow(df.beds)) ## binarize PTM regions ## if we dont want to use signal

## Data wrangle
chr.df <- df.beds
chr.df <- chr.df[,-c(4:6)]

## Data wrangle more
chr.df <- chr.df  %>% dplyr::group_by(chr, cond) %>% dplyr::reframe(signal = sum(signalValue))
chr.df <- chr.df[-c(93:94),]
chr.df$factor<- as.numeric(ifelse(chr.df$chr == "chrX" , 23, sapply(strsplit(as.character(chr.df$chr), "chr"), `[`, 2)))
chr.df <- tidyr::separate_wider_delim(chr.df,cond,delim = "_",names = c( "mark", "status"), cols_remove=F)

## Labs
labs <- seq(1,22,1)
labs <- append(labs,"X")

## Plot per chrom change for SI
p<-chr.df %>%
  mutate(chr= forcats::fct_reorder(as.character(chr),as.numeric(factor))) %>% 
  ggplot(aes(x=chr, y=log2(signal),group = cond, color=cond)) +
  geom_line(aes(linetype= cond, color=cond))+
  geom_point(aes(color=cond),size = 2)+
  scale_linetype_manual(values=c("twodash", "solid", "twodash", "solid"))+ 
  scale_x_discrete(labels = labs)

p2 <- p + facet_grid( ~mark)
p2+ theme_classic() + theme(axis.text.x = element_text(angle = 90,  hjust=1), legend.position="bottom", plot.title = element_text(size=18), text = element_text(size=18, family="Arial"))

##---------------------------------------- Get delta -##

chr.ct <- chr.df[chr.df$status == "ctrl",]
chr.ko <- chr.df[chr.df$status == "6hrsAux",]

delta <- data.frame(chr = chr.ct$chr, mark = chr.ct$mark, ctrl.signal = log2(chr.ct$signal), ko.signal = log2(chr.ko$signal)) %>% 
  dplyr::mutate(diff = ko.signal-ctrl.signal) %>% 
  dplyr::mutate(delta = ((ctrl.signal+diff)/ctrl.signal))

##---------------------------------------- Plot delta for main text figure-##

pal <- c("#FF8C00", "#A034F0", "#159090")
theme_set(theme_classic(base_size = 18))
ggplot(delta, aes(mark, delta)) + 
  ggdist::stat_halfeye(aes(color = mark,
                           fill = after_scale(colorspace::lighten(color, .1))),
                       adjust = .4, width = .8, .width = 0, justification = -.3, point_colour = NA) + 
  geom_boxplot(width = .35, outlier.shape = NA) +
  #gghalves::geom_half_point(side = "l", range_scale = 0, shape = 95, size = 15, alpha = .3, width = 0.8)+
  scale_color_manual(values = pal, guide = "none") +
  scale_fill_manual(values = pal, guide = "none") + 
  scale_y_continuous(limit = c(0.8,1.35),breaks = seq(0,1.35,0.1),labels = scales::percent)+
  stat_summary(
    geom = "text",
    fun = "median",
    aes(label = paste0(round(..y.., 2)*100, "%"),
        color = mark,
        color = after_scale(colorspace::darken(color, .2, space = "HLS"))),
    size = 6,
    vjust = -3.5
  )+ scale_x_discrete(labels = c("H3K27ac", "H3K9me3"))+
  labs(
    x = NULL,
    y = "Δ%",
    title = "Δ% Per Chromosome"
  )


#########################################################################################
# Proseq: merge count replicates
#########################################################################################

bed.dir <- file.path(root.dir, "proseq")

name <- c("6hrsAux","ctrl")
bed.r1 <- as.data.frame(read.table(file.path(bed.dir,paste0(name[1],"_rep1.tsv")), sep= "\t", header= T))
bed.r2 <- as.data.frame(read.table(file.path(bed.dir,paste0(name[1],"_rep1.tsv")), sep= "\t", header= T))

bed.r1 <- as.data.frame(read.table(file.path(bed.dir,paste0(name[1],"_merged.tsv")), sep= "\t", header= T))
bed.r2 <- as.data.frame(read.table(file.path(bed.dir,paste0(name[2],"_merged.tsv")), sep= "\t", header= T))

colnames(bed.r1) <-c("chr", "start", "end", "GENEID", "TXNAME", "TXSTRAND", "signal") ## Pro seq
colnames(bed.r2) <-c("chr", "start", "end", "GENEID", "TXNAME", "TXSTRAND", "signal") ## Pro se

new.bed <- bed.r1 %>% mutate(signal = (bed.r1$signal + bed.r2$signal)/2)

library("readr")
write_tsv(new.bed, file = paste0(bed.dir, name[1],"_merged.tsv"), col_names=T) ## merge signal from each rep

#########################################################################################
# Proseq analysis
#########################################################################################

## Loop through tads and loops to get nearest distance and signal
df.lps <- data.frame()
conds <- c("6hrsAux_merged.tsv", "ctrl_merged.tsv")
for (i in 1:length(conds)){
  
  cond <- conds[i]
  print(cond)
  
  x <- getOverlap(dir.bed = bed.dir, bed = cond,interval=loops, cond = paste0(cond, "_loops"),anchors = F)
  df.lps <- rbind(df.lps, x)
}

df.tads <- data.frame()
for (i in 1:length(conds)){
  
  cond <- conds[i]
  print(cond)
  
  y <- getOverlap(dir.bed = bed.dir, bed = cond,interval=tads, cond = paste0(cond, "_tads"),anchors = F)
  df.tads <- rbind(df.tads, y)
}

## Reframe data
n.tds <- df.tads %>% dplyr::group_by(coord) %>% dplyr::reframe(signal = sum(signal))
n.lps <- df.lps %>% dplyr::group_by(coord) %>% dplyr::reframe(signal = sum(signal))
n.df <- rbind(n.tds,n.lps)

## Prepare data for plotting
n.df <- tidyr:: separate_wider_delim(n.df,coord,delim = "_",names = c("chrom","start","end", "status","rep","top"), cols_remove=F)
n.df <- n.df %>% mutate(size = abs(as.numeric(start)-as.numeric(end)))
n.df$bin <- dplyr::ntile(n.df$size, as.numeric(7))
n.df$cond <- paste0( n.df$status, " ", n.df$top, " ",n.df$bin)

## Get mean bin size for loops
size.lab <- n.df %>% dplyr::group_by(bin) %>% dplyr::reframe(mean.size = mean(size)) 

lim = 10000 ## limit for axes
ggplot(n.df, aes(x = log10(signal), y = cond, fill = stat(x))) +
  geom_density_ridges_gradient(scale = 3, rel_min_height = 0.01) +
  scale_fill_viridis_c(name = "Signal", option = "C") +
  labs(title = 'h3k27me3') + scale_x_continuous(limits = c(0,lim), breaks = seq(0,lim,1))

ggplot(n.df, aes(x = signal, y=cond, fill = factor(stat(quantile)))) +
  stat_density_ridges(
    geom = "density_ridges_gradient", calc_ecdf = TRUE,
    quantiles = 4, quantile_lines = TRUE
  ) + scale_x_continuous(limits = c(0,lim), breaks = seq(0,lim,100))+
  scale_fill_viridis_d(name = "Quartiles") + theme_classic()

#########################################################################################
# CHIA PET analysis
#########################################################################################

## bring in CHIA PET loops
loop.dir <- file.path(root.dir,"chiapet")
bed.dir <- file.path(root.dir,"rad21_chip")

name <- c("6hrsAux","ctrl")
loops.ctrl <- read.table(file.path(loop.dir,paste0("Pol2_",name[2],".bedpe")), sep= "\t", header= F)
loops.ko <- read.table(file.path(loop.dir,paste0("Pol2_",name[1],".bedpe")), sep= "\t", header= F)

## conds
conds <- c("h3k9me3_ctrl.bed", "h3k9me3_6hrsAux.bed")
conds <- c( "h3k27ac_ctrl.bed", "h3k27ac_6hrsAux.bed")

## Loop through tads and loops to get nearest distance and signal
df.ct <- data.frame()
for (i in 1:length(conds)){
  
  cond <- conds[i]
  print(cond)
  
  x <- getOverlap(dir.bed = bed.dir, bed = cond,interval=loops.ctrl, cond = paste0(cond, "_ct"),anchors = F)
  df.ct <- rbind(df.ct, x)
}

## Bring in TADs
df.ko <- data.frame()
for (i in 1:length(conds)){
  
  cond <- conds[i]
  print(cond)
  
  y <- getOverlap(dir.bed = bed.dir, bed = cond,interval=loops.ko, cond = paste0(cond, "_ko"),anchors = F)
  df.ko <- rbind(df.ko, y)
}

## Reframe data
n.ct <- df.ct %>% dplyr::group_by(coord) %>% dplyr::reframe(signal = sum(signalValue))
n.ko <- df.ko %>% dplyr::group_by(coord) %>% dplyr::reframe(signal = sum(signalValue))
n.df <- rbind(n.ct,n.ko)

## Prepare data for plotting
n.df <- tidyr::separate_wider_delim(n.df,coord,delim = "_",names = c("chrom","start","end", "mark", "mark_status","chia_status"), cols_remove=F)
n.df <- n.df %>% dplyr::mutate(size = abs(as.numeric(start)-as.numeric(end)))
n.df$bin <- dplyr::ntile(n.df$size, as.numeric(7))
n.df$cond <- paste0(n.df$mark," ", n.df$mark_status, " ", n.df$chia_status, " ",n.df$bin)


##---------------------------------------- Generating SI and main text plots -##

# n.df.k9 <- n.df
# n.df.k27 <- n.df

saveRDS(n.df.k9, file.path(loop.dir,"n.df.k9.rds"))
saveRDS(n.df.k27, file.path(loop.dir,"n.df.k27.rds"))

# Load the city object as city
n.df.k9 <- readRDS(file.path(loop.dir,"n.df.k9.rds"))
n.df.k27 <- readRDS(file.path(loop.dir,"n.df.k27.rds"))

data = n.df.k9

data <- data %>% filter(bin %in% c(1,4,5)) ## for maintext plots k27
title = "H3K27ac Coverage"

data <- data %>% filter(bin %in% c(2,3,4)) ## for maintext plots k9
title = "H3K9me3 Coverage"

## Modify labels for legibility
data$chia_status <- ifelse(data$chia_status == "ct", "CHIA-PET: DMSO", "CHIA-PET: 6 Hrs Aux")
data$mark_status <- ifelse(data$mark_status == "ctrl", "ChIP: DMSO", "ChIP: 6 Hrs Aux")

## plotting parameters
lim = 17

## Get mean bin size for loops and generate labels
size.lab <- data %>% group_by(bin) %>% reframe(mean.size = mean(size)) 
lab <- paste0(format(as.integer(rep(round_any(size.lab$mean.size,1000),4)),big.mark=",",scientific=F, trim = T), " bp") ## lab for facet plot
lab <- paste0(sapply(strsplit(as.character(lab), ","), `[`, 1), " KBP")


## SI 1
# k9.chia.v1 # k27.chia.v1
p1 <- ggplot(data, aes(x = log2(signal), y = as.factor(bin), fill = stat(x))) +
  geom_density_ridges_gradient(scale = 2, rel_min_height = 0.01) +
  scale_fill_viridis_c(name = "Signal", option = "C") +
  labs(title = title) + scale_x_continuous(limits = c(0,lim), breaks = seq(0,lim,5))+
  scale_y_discrete( labels = lab) + ylab("") +
  theme(axis.text.x = element_text(angle = 0,  hjust=1), legend.position="", plot.title = element_text(size=16), text = element_text(size=16, family="Arial"))
p2 <- p1 + facet_grid(fct_rev(mark_status) ~ fct_rev(chia_status),axes = "all", axis.labels = "all_x")
p2 + theme_classic() + theme(axis.text.x = element_text(angle = 0,  hjust=1), legend.position="", plot.title = element_text(size=18), text = element_text(size=18, family="Arial"))

## SI 2
# k9.chia.v2 # k27.chia.v2
p1 <- ggplot(data, aes(x = log2(signal), y=as.factor(bin), fill = factor(stat(quantile)))) +
  stat_density_ridges(
    geom = "density_ridges_gradient", calc_ecdf = TRUE,
    quantiles = 4, quantile_lines = TRUE
  ) + scale_fill_viridis_d(name = "Quartiles") +labs(title = title) +
  labs(title = title) + scale_x_continuous(limits = c(0,lim), breaks = seq(0,lim,5))+
  scale_y_discrete( labels = lab) + ylab("")
p2 <- p1 + facet_grid(fct_rev(mark_status) ~ fct_rev(chia_status),axes = "all", axis.labels = "all_x")
#p2 + theme_minimal() + theme(axis.text.x = element_text(angle = 0,  hjust=1), legend.position="", plot.title = element_text(size=16), text = element_text(size=16, family="Arial"))
p2 + theme_classic() + theme(axis.text.x = element_text(angle = 0,  hjust=1), legend.position="", plot.title = element_text(size=18), text = element_text(size=18, family="Arial"))

##---------------------------------------- Generating main text combined k9/k27 plot -##

library(plyr)

## Wranlge data
data.1 = n.df.k9
data.2 = n.df.k27

data <- rbind(data.1, data.2)

data$bin <- ntile(data$size, as.numeric(7))
data$cond <- paste0(data$mark," ", data$mark_status, " ", data$chia_status, " ",data$bin)

data <- data %>% filter(bin %in% c(1,4,5)) ## for maintext plots k27

## Modify labels for legibility
data$chia_status <- ifelse(data$chia_status == "ct", "CHIA-PET: DMSO", "CHIA-PET: 6 Hrs Aux")
data$mark_status <- ifelse(data$mark_status == "ctrl", "ChIP: DMSO", "ChIP: 6 Hrs Aux")

## plotting parameters
lim = 17
title = "CHIA-PET Loop PTM Coverage"

## Get mean bin size for loops and generate labels
size.lab <- data %>% group_by(bin) %>% reframe(mean.size = mean(size)) 
lab <- paste0(format(as.integer(rep(round_any(size.lab$mean.size,1000),4)),big.mark=",",scientific=F, trim = T), " bp") ## lab for facet plot
lab <- paste0(sapply(strsplit(as.character(lab), ","), `[`, 1), " KBP")

data <- data[data$chia_status == "CHIA-PET: DMSO",] ## exclude non control loops for final figure

## Main figure
# k9.chia.v2 # k27.chia.v2
p1 <- ggplot(data, aes(x = log2(signal), y=as.factor(bin), fill = factor(stat(quantile)))) +
  stat_density_ridges(
    geom = "density_ridges_gradient", calc_ecdf = TRUE,
    quantiles = 4, quantile_lines = TRUE
  ) + scale_fill_viridis_d(name = "Quartiles", option = "B") +labs(title = title) +
  labs(title = title) + scale_x_continuous(limits = c(0,lim), breaks = seq(0,lim,5))+
  scale_y_discrete( labels = lab) + ylab("")
p2 <- p1 + facet_grid(fct_rev(mark_status) ~ fct_rev(mark),axes = "all", axis.labels = "all_x")
#p2 + theme_minimal() + theme(axis.text.x = element_text(angle = 0,  hjust=1), legend.position="", plot.title = element_text(size=16), text = element_text(size=16, family="Arial"))
p2 + theme_classic() + theme(axis.text.x = element_text(angle = 0,  hjust=1), legend.position="", plot.title = element_text(size=18), text = element_text(size=18, family="Arial"))

#########################################################################################
# association analysis
#########################################################################################

dir <- file.path(root.dir, "2_coll_STORM")

## Wrangle association data for plotting ## we're using the 60nm results
assoc <- data.frame(t(read.table(file.path(dir, "rad21_edu_60nm.csv"), sep= ",", header= F)))

## drop columns and rows
assoc <- assoc[,-c(2:3)]
assoc <- assoc[-c(1),]
rownames(assoc) <- c(1:nrow(assoc))

## convert columns to numeric
i <- c(2:4)  
assoc[, i] <- apply(assoc[, i], 2, function(x) as.numeric(as.character(x)))
sapply(assoc, class) ## check class

##---------------------------------------- Generating SI plot for all groups  -##

assoc <- assoc %>% dplyr::mutate(Unassociated = (1- rowSums(assoc[ , c(2:4)], na.rm=TRUE))) %>% dplyr::mutate(X1 = rownames(assoc))
colnames(assoc) <- c("X", "Small Cluster", "Large Cluster", "Both", "Unassociated")
assoc <- assoc[,c("X",  "Unassociated","Large Cluster","Small Cluster",  "Both")]

assoc <- na.omit(assoc) %>% tidyr::gather(key=loc, value = value, 2:5) %>% dplyr::mutate(Cell = as.numeric(X)) %>%  dplyr::rename(factor=X)

## Set colors for individual plots
v_colors =  viridis(6, option = "D")
v_colors

assoc<- ddply(assoc , "Cell",transform, label_ypos=cumsum(value))
#assoc <- assoc[!assoc$loc == "Unassociated",] ## for unassociated removed
assoc$label_ypos <- ifelse(assoc$value < 0.06, NA, assoc$label_ypos) ## for full data

## SI plot
assoc %>% 
  mutate(Cell = forcats::fct_reorder(as.character(Cell), as.numeric(factor))) %>% 
  ggplot(aes(x = Cell, y = value, fill = loc)) + 
  geom_bar(stat="identity", position="stack", color = NA) + scale_fill_manual(values = rev(v_colors)) + theme_classic()+
  theme( legend.position="bottom",plot.title = element_text(size=16), text = element_text(size=16, family="Arial"), axis.text.x = element_text(angle = 0)) +
  ggtitle("SMLM: RAD21 & EDU") + xlab("Cell") + ylab("Rad21 Domain Association") +scale_y_continuous(limit = c(0,1),breaks = seq(0,1,0.1),labels = scales::percent)+
  geom_text(aes(y=label_ypos, label= paste0(round(value*100, 1),"%")), vjust=1.6, color="white", size=3.5)


##---------------------------------------- Generating SI plot - no Unassociated -##

## Wrangle association data for plotting
assoc2 <- data.frame(t(read.table(file.path(dir, "rad21_edu_60nm.csv"), sep= ",", header= F)))

## drop columns and rows
assoc2 <- assoc2[,-c(2:3)]
assoc2 <- assoc2[-c(1),]
rownames(assoc2) <- c(1:nrow(assoc2))

## convert columns to numeric
i <- c(2:4)  
assoc2[, i] <- apply(assoc2[, i], 2, function(x) as.numeric(as.character(x)))
sapply(assoc2, class) ## check class

assoc2 <- assoc2 %>% dplyr::mutate(cum = (rowSums(assoc2[ , c(2:4)], na.rm=TRUE))) %>% dplyr::mutate(X1 = rownames(assoc2))
assoc2[,c(2:4)] <- assoc2[,c(2:4)]/assoc2[,c(5)]
colnames(assoc2) <- c("X", "Small Cluster", "Large Cluster", "Both")
assoc2 <- assoc2[,c("X","Large Cluster","Small Cluster",  "Both")]

assoc2 <- na.omit(assoc2) %>% tidyr::gather(key=loc, value = value, 2:4) %>% dplyr::mutate(Cell = as.numeric(X)) %>%  dplyr::rename(factor=X)

## Set colors for individual plots
v_colors =  viridis(5, option = "D")# E= Civis
v_colors

plot<- ddply(assoc2 , "Cell",transform, label_ypos=cumsum(value))
plot$label_ypos <- ifelse(plot$value < 0.06, NA, plot$label_ypos) ## for full data
plot$factor <- as.numeric(plot$factor) + as.numeric(rep(c(0.1,0.2,0.3),12))

## SI.scaled.bar
plot%>% 
  dplyr::mutate(loc = forcats::fct_reorder(as.character(loc), as.numeric(factor))) %>% 
  dplyr::mutate(Cell = forcats::fct_reorder(as.character(Cell), as.numeric(factor))) %>% 
  ggplot(aes(x = Cell, y = value, fill = forcats::fct_rev(loc))) + 
  geom_bar(stat="identity", position="stack", color = NA) + scale_fill_manual(values = v_colors) + theme_classic()+
  ggtitle("SMLM: RAD21 & EDU") + xlab("Cell") + ylab("Scaled Rad21 Domain Association") +scale_y_continuous(limit = c(0,1),breaks = seq(0,1,0.1),labels = scales::percent)+
  geom_text(aes(y=label_ypos, label= paste0(round(value*100, 1),"%")), vjust=1.6, color="white", size=3.5) + 
  theme( legend.position="bottom",plot.title = element_text(size=16), text = element_text(size=16, family="Arial"), axis.text.x = element_text(angle = 0)) 

##---------------------------------------- Generating main text plots 1 -##

summary <- assoc %>% ## include Unassociated, unscaled values
  dplyr::group_by(loc) %>% dplyr::reframe(mean_N=mean(value),  
                                               sd_N=sd(value),  
                                               N_N=n(),  
                                               se=sd_N/sqrt(N_N),  
                                               upper_limit=mean_N+se,  
                                               lower_limit=mean_N-se)

summary <- summary[-c(4),] ## remove unassociated ## for assoc only

## main.bar.all
ggplot(summary, aes(x=loc, y=mean_N, fill = loc)) +  scale_fill_manual(values = v_colors)+scale_color_manual(values = v_colors)+
  geom_bar(stat="identity", color = "black", size = 2) +  scale_y_continuous(limit = c(0,.08),breaks = seq(0,0.08,0.01),labels = scales::percent) +
  geom_errorbar(aes(ymin=lower_limit, ymax=upper_limit, color = loc), width = 0.5, size = 1.5) + theme_classic() + ylab("Domain Association") + xlab("")+
  theme( legend.position="bottom",plot.title = element_text(size=18), text = element_text(size=18, family="Arial"), axis.text.x = element_text(angle = 0))

##---------------------------------------- Generating main text plots 2 -##
## Single bar plot summary

## Wrangle association data for plotting
assoc <- data.frame(t(read.table(file.path(dir, "rad21_edu_60nm.csv"), sep= ",", header= F)))

## drop columns and rows
assoc <- assoc[,-c(2:3)]
assoc <- assoc[-c(1),]
rownames(assoc) <- c(1:nrow(assoc))

## convert columns to numeric
i <- c(2:4)  
assoc[, i] <- apply(assoc[, i], 2, function(x) as.numeric(as.character(x)))
sapply(assoc, class) ## check class

assoc <- assoc %>% dplyr::mutate(cum = (rowSums(assoc[ , c(2:4)], na.rm=TRUE))) %>% dplyr::mutate(X1 = rownames(assoc))
colnames(assoc) <- c("X", "Small Cluster", "Large Cluster", "Both", "All")
assoc <- na.omit(assoc) %>% tidyr::gather(key=loc, value = value, 2:5) %>% dplyr::mutate(Cell = as.numeric(X)) %>%  dplyr::rename(factor=X)

summary <- assoc %>% 
  dplyr::group_by(loc) %>%## for no Unassociated, scaled values
  dplyr::reframe(mean_N=mean(value),
                            sd_N=sd(value),  
                            N_N=n(),  
                            se=sd_N/sqrt(N_N),  
                            upper_limit=mean_N+se,  
                            lower_limit=mean_N-se)
summary <- summary[c(1),]

## Main text plot
ggplot(summary, aes(x=loc, y=mean_N)) + 
  geom_bar(stat="identity", color = "black", fill = "white",size = 2, width = 0.5) +  scale_y_continuous(limit = c(0,1),breaks = seq(0,1,0.1),labels = scales::percent) +
  geom_errorbar(aes(ymin=lower_limit, ymax=upper_limit), width = 0.2, size = 1.5) + theme_classic() + ylab("Domain Association") + xlab("")+
  theme( legend.position="bottom",plot.title = element_text(size=16), text = element_text(size=16, family="Arial"), axis.text.x = element_text(angle = 0))+
  geom_text(aes(y=mean_N+0.12, label= paste0(round(mean_N*100, 2),"%")), vjust=1.6, color="black", size=5)

