gc()

# import packages
library(ggplot2)
library(viridis)
library(dplyr)
library(tidyverse)
library(forcats)
library(GENOVA)

## Broadly, colours of heatmaps can be changed by setting the GENOVA.colour.palette option.
## https://github.com/robinweide/GENOVA/issues/298
options("GENOVA.colour.palette" = "whitered")

## -- introduce data -- ##

## -- micro c data -- ##
dir <- "/projects/b1042/BackmanLab/HiC2/opt/juicer/work/022324_microC/cools"
list.files(paste0(dir), all.files=F, include.dirs = FALSE)

wt.5kb <- load_contacts(signal_path = "/projects/b1042/BackmanLab/HiC2/opt/juicer/work/022324_microC/juicer_analysis/CTRL/mega/aligned/inter_30_5k.cool",
                        sample_name = "WT",
                        resolution = 5e3,
                        balancing = F, # this is the default
                        colour = "black")


rad21.5kb <- load_contacts(signal_path = "/projects/b1042/BackmanLab/HiC2/opt/juicer/work/022324_microC/juicer_analysis/6hrs_Aux/mega/aligned/inter_30_5k.cool",
                           sample_name = "Rad21",
                           resolution = 5e3,
                           balancing = F, # this is the default
                           colour = "black")

wt.100kb <- load_contacts(signal_path = "/projects/b1042/BackmanLab/HiC2/opt/juicer/work/022324_microC/juicer_analysis/CTRL/mega/aligned/inter_30_100k.cool",
                          sample_name = "WT",
                          resolution = 100e3,
                          balancing = F, # this is the default
                          colour = "black")


rad21.100kb <- load_contacts(signal_path = "/projects/b1042/BackmanLab/HiC2/opt/juicer/work/022324_microC/juicer_analysis/6hrs_Aux/mega/aligned/inter_30_100k.cool",
                             sample_name = "Rad21",
                             resolution = 100e3,
                             balancing = F, # this is the default
                             colour = "black")
gc()
## -- Triangle plots -- ##

## These loci were all taken from Rao et al 2017 figure 1a

## Loci from chr8 133800000 134600000 (hg19) converted to chr8 132787754	133587757 (hg38)
pyramid(exp = wt.5kb, #rad21.chr8.133-134.pyr
        chrom = "8",
        start = 132787754,
        end=133587757)

pyramid(exp = rad21.5kb,
        chrom = "chr8",
        start = 132787754,
        end=133587757)

## Loci from chr4 40800000 42100000 (hg19) converted to chr4 40797983	42097983 (hg38)
pyramid(exp = wt.5kb, 
        chrom = "4",
        start = 40797983,
        end=42097983)

pyramid(exp = rad21.5kb, #wt.chr4.40-42.pyr
        chrom = "chr4",
        start = 40797983,
        end=42097983)

## Loci from chr1 91900000 95800000 (hg19) converted to chr1	91434443	95334444 (hg38)
pyramid(exp = wt.5kb, #rad21.chr1.91-95.pyr
        chrom = "1",
        start = 91434443,
        end=95334444)

pyramid(exp = rad21.5kb,
        chrom = "chr1",
        start = 91434443,
        end=95334444)

## KO pyramid using first loci

library(patchwork)

wt.pyr <- pyramid(exp = wt.5kb,
                  chrom = "8",
                  start = 132787754,
                  end=133587757) + ggplot2::ggtitle("Wildtype")

ko.pyr <- pyramid(exp = rad21.5kb,
                  chrom = "8",
                  start = 132787754,
                  end=133587757) + ggplot2::ggtitle("Rad21 KO")

ko.pyr / wt.pyr + plot_layout(guides = "collect")


## Loci from chr14:68200001-69500000 (hg19) converted to chr14	67733283	69033283 (hg38)
library(patchwork)

wt.pyr <- pyramid(exp = synced[[1]],
                  chrom = "14",
                  start = 67733283,
                  #colour = c(0,75),
                  end=69033283) + ggplot2::ggtitle("Wildtype")

ko.pyr <- pyramid(exp = synced[[2]],
                  chrom = "14",
                  start = 67733283,
                  #colour = c(0,75),
                  end=69033283) + ggplot2::ggtitle("Rad21 KO")

#pyr.chr14.6.7mb.chr6.9
ko.pyr / wt.pyr + plot_layout(guides = "collect")

## Difference pyramid
pyramid_difference(
  exp1 = synced[[2]],
  exp2 = synced[[1]],
  chrom = "chr14",
  start = 67733283,
  end=69033283)

pyramid_difference(
  exp1 = synced[[2]],
  exp2 = synced[[1]],
  chrom = "14",
  start = 67733283,
  end=69033283)


## -- domainogram -- ##

## Check ideal window size for optimal tad calling

ID <- insulation_domainogram(
  wt,
  chrom = '7',
  start = 25e6,
  end = 29e6,
  window_range = c(1, 101),
  step = 2
)
visualise(ID)

## -- TAD analysis -- ##

## sync indices between conditions
synced <- sync_indices(list(wt.5kb, rad21.5kb))
synced <- sync_indices(list(wt.100kb, rad21.100kb))
gc()

rm(wt.5kb, rad21.5kb)
rm(wt.25kb, rad21.25kb)

lapply(synced, attributes) ## Get all attributes from list
synced[[1]][1] ## Access first element in list of wt

rm(aux,wt) ## declutter workspace as needed

## Checking different insulation windows using synced list
insulation.w25 <- insulation_score(synced,window = 25) ## Window size 25 seems to be the most relevant when tested, for 5kb resolution
insulation.w40 <- insulation_score(synced,window = 40)

df <- data.frame(insulation.w40$insula_score)
#write.csv(df, "/home/lmc0633/rad21paper/cool.plots/insulation.csv", row.names=FALSE)
write.csv(df, "/home/lmc0633/PM paper/insulation.csv", row.names=FALSE)

## Inspect insulation windows
visualise(insulation.w25,
          chr = '4', start = 91434443, end=95334444,
          contrast = 2)

## Check INF and NA values
check <- data.frame(insulation.w25[1])
sum(is.infinite(check$insula_score.WT))
sum(is.na(check$insula_score.WT))
nrow(check)

## Call TADs on different window sizes #himatrix.91.4-95.3.diff.5kb.50w
TADcalls.25 <- call_TAD_insulation(insulation.w25)## use window size 25
TADcalls.40 <- call_TAD_insulation(insulation.w40) 


wt <- data.frame(TADcalls.40$WT, cond = rep("WT",nrow(TADcalls.40$WT)))
r21 <- data.frame(TADcalls$Rad21, cond = rep("Rad21",nrow(TADcalls$Rad21)))
df <- rbind(wt, actd)

#write.csv(df, "/home/lmc0633/rad21paper/cool.plots/GenovaTADs.csv", row.names=FALSE)
write.csv(df, "/home/lmc0633/PM paper/GenovaTADs.csv", row.names=FALSE)

## Read in TAD data
wt.domains= read.delim('/projects/b1042/BackmanLab/HiC2/opt/juicer/work/022324_microC/contact_domains/CTRL/mega/inter_30_contact_domains/domains.bedpe', h = F, skip = 1)
rad21.domains = read.delim('/projects/b1042/BackmanLab/HiC2/opt/juicer/work/022324_microC/contact_domains/6hrs_Aux/mega/inter_30_contact_domains/domains.bedpe', h = F, skip = 1)

## Prepare tad ddata
wt.domains <- wt.domains[-c(1),]
wt.domains <- wt.domains[,c(1:3)]

rad21.domains <- rad21.domains[-c(1),]
rad21.domains <- rad21.domains[,c(1:3)]

## Plot a test matrix
hic_matrixplot(exp1 = synced[[1]],
               exp2 = synced[[2]],
               chrom = '4',
               start = 91434443,
               end=95334444,
               #tads = list(TADcalls$Rad21, TADcalls$WT), # see ATA
               #tads = list(rad21.domains, wt.domains),
               #tads.type = list('lower','upper'),
               #tads.colour = c('red','#91cf60'), # green TAD-borders
               cut.off = 10, # upper limit of contacts
               skipAnn = T) # skip the outside annotation

## Plot a test matrix
hic_matrixplot(exp1 = wt.5kb,
               exp2 = rad21.5kb,
               chrom = "chr14",
               start = 67733283,
               end=69033283,
               #tads = list(TADcalls$Rad21, TADcalls$WT), 
               tads = list(wt.domains, rad21.domains),# see ATA
               tads.type = list('lower','upper'), # only plot in lower triangle
               tads.colour = c('red','#91cf60'), # green TAD-borders
               cut.off = 20, # upper limit of contacts
               skipAnn = T) # skip the outside annotation


## aggregate TAD analysis
ATA.ah <- ATA(list(WT= wt.5kb, Rad21 = rad21.5kb), bed = wt.domains) ## Arrowhead TADs
ATA.gn <- ATA(synced, bed = TADcalls$WT) ## Genova TADs
ATA

visualise(ATA.gn,
          colour_lim = c(0,20),
          colour_lim_contrast = c(-5,5),
          metric = "diff",
          focus = 2
)

## Custom colors ATA plot

## Colors for difference plot
colors <- colorRampPalette(c("#443A83FF", "white","black"))(12)

## Colors for ATA part
v_colors <- c("white","#F4EDCA", "#EAD357FF", "#D0BE67FF", "#B6A971FF", "#9E9677FF", "#878479FF", "#727374FF", "#5E626EFF", "#48526BFF", "#2A406CFF", "#00306FFF", "#00204DFF")

## https://github.com/robinweide/GENOVA/issues/298
visualise(ATA.ah, colour_lim = c(0,30),metric = "diff", focus = 2) +
  ggplot2::scale_fill_gradientn(
    colours = c(v_colors),limits= c(0,25), na.value="black"
  ) +
  ggplot2::scale_colour_gradientn(limits = c(-5,5),
                                  aesthetics = "altfill",
                                  colours = c(colors),na.value="#365C8DFF",
                                  guide = ggplot2::guide_colourbar(available_aes = "altfill")
  )


## -- TAD analysis - 2D Saddleplot -- ##
## Taken from https://github.com/robinweide/GENOVA/issues/272

library(data.table)
library(ggplot2)

ATA$signal_raw$Rad21
## function written for 2 samples
saddle.ata <- function(ata) {
  
  # Melt array
  df <- data.table(
    x = as.vector(slice.index(ata$signal, 1)),
    y = as.vector(slice.index(ata$signal, 2)),
    sample = as.vector(slice.index(ata$signal, 3)),
    value = as.vector(ata$signal)
  )
  # Calculate distance (unit is arbitrary due to ATA)
  df[, dist := x - y]
  # Take averages per off-diagonal band
  df <- df[, list(value = mean(value)), by = c("dist", "sample")]
  # Split by sample
  df <- split(df, df$sample)
  
  # Adjust sample 2 a bit so that a difference will show
  # Don't do this with real data!
  #df[[2]][, value := sqrt(value)]
  
  # Recombine the two samples, if you have more, this needs to be repeated for all samples
  df <- df[[1]][df[[2]], on = c("dist")]
  # Convert distance to TAD units
  df[, dist := scales::rescale(dist, to = c(-2, 2))]
  
  return(df)
}

df <- saddle.ata(ATA.ah)

## Prepare data for plotting
sample1 <- data.frame(dist = df$dist, cond = rep("wt", nrow(df)), value = df$value)
sample2 <- data.frame(dist = df$dist, cond = rep("rad21", nrow(df)), value = df$i.value)

plot.df <- rbind(sample1, sample2)

# Plot 2D saddle
p <- plot.df %>%
  ggplot(aes(x=dist , y=value)) + 
  geom_line(size = 0.5, aes(color = cond)) +theme_classic()+
  #scale_y_continuous(limits = c(-2.5,  2))+
  #scale_x_continuous(limits = c(5,8.5), breaks = c(4, 4.5,5, 5.5, 6,6.5,7,7.5, 8, 8.5))+
  scale_color_manual(values=c("#C70039", "black")) +
  labs(title="Aggregate TAD Analysis", subtitle="Wt vs Rad21", 
       caption="")+theme(text = element_text(size = 12),axis.text.x = element_text(angle = 0, vjust = 1, hjust=1, size = 12),
                         axis.text.y = element_text(angle = 0,size = 12))+ 
  xlab("Distance From Diagonal (AU)") + 
  ylab("Signal")
p

## -- Relative Contact Probability -- ##

rcp <- RCP(explist = list(wt.5kb, rad21.5kb),
           chromsToUse = NULL, maxDistance = 20e6)

gc()
## Plot RCP
visualise(rcp) +scale_color_manual(values=c("#440154FF" ,"#21908CFF"))
#visualise(rcp) +scale_color_manual(values=c("#f8a07e" ,"#a059a0"))

## PLot log2 differential of RCP between conditions
visualise(rcp, contrast = 2, metric = 'lfc') +scale_color_manual(values=c("#f8a07e"))

## -- trans interactions -- ##

## Matrix plot of all trans interactions
cm <- chromosome_matrix(list(WT = synced.100kb[[1]], KO = synced.100kb[[2]]),include_chr = "all",expected = "trans", sort_chr = TRUE)

labels <- as.character(c(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,"X"))

visualise(cm) + scale_y_discrete(labels = labels)+ scale_x_discrete(labels = labels)


## -- Compartments -- ##

## Group compartments into a single list
synced.100kb <- sync_indices(list(wt.100kb, rad21.100kb))
gc()

## Get eigens
cs = compartment_score(synced.100kb)

## comp strength - down: 3, 5, 10, 13, 19, 20  up: 7,8, 11, 14, 17,18 (MICRO-C rad21)
visualise(cs, chr = "chr22") +scale_color_manual(values=c("black" ,"#C70039"))

## comp strength - down: 2, 3, 4,5, 7, 9, 11, 13, 16, 17, 18, 19, 
visualise(cs, chr = "X") +scale_color_manual(values=c("black" ,"#C70039"))

## Saddle plot calculation
attr(cs, 'signed') <- T # Fake signed CS
saddle.plot = saddle(synced.100kb,
                     CS_discovery = cs,
                     bins = 50)

p <- visualise(saddle.plot)

## Colors for difference plot
colors <- colorRampPalette(c("#443A83FF", "white","black"))(12)

## Colors for saddle part
v_colors <- c("white","#F4EDCA", "#EAD357FF", "#D0BE67FF", "#B6A971FF", "#9E9677FF", "#878479FF", "#727374FF", "#5E626EFF", "#48526BFF", "#2A406CFF", "#00306FFF", "#00204DFF")


## Remap for plotting with different color scale
p$layers[[1]]$mapping <- ggplot2::aes(fill = log2(obsexp))
lims <- range(log2(p$data$obsexp))
values <- scales::rescale(c(lims[1], 0, lims[2]))

## Plot saddle
p + ggplot2::scale_fill_gradientn(
  colours = c(v_colors),
  values  = values,
  limits = lims,
  name = "Log2(FC)"
) +ggplot2::scale_colour_gradientn(limits = c(-1,1),
                                   aesthetics = "altfill",
                                   colours = c(colors),na.value="#365C8DFF",
                                   guide = ggplot2::guide_colourbar(available_aes = "altfill"))


## Rain plot of Act D compartment strength
library(ggridges)
library(ggrain)

CSS <- quantify(saddle.plot)
#write.csv(CSS, "/home/lmc0633/rad21paper/cool.plots/compartmentscores.csv", row.names=FALSE)
write.csv(CSS, "/home/lmc0633/PM paper/compartmentscores.csv", row.names=FALSE)

## Prepare data
compared <- tidyr::spread(unique(CSS[,-c(3,4)]), key = 'exp', value = 'strength')
plot.df <- compared %>% 
  gather(cond, value, 2:3) 

## Cut outliers
plot.df <- plot.df[plot.df$value < 10,]

plot.df %>%
  ggplot( aes(1, y=value, fill=cond, color = cond ))+
  geom_rain(alpha = .5, rain.side = 'l',trim = F,
            boxplot.args = list(color = "black", outlier.shape = NA),
            boxplot.args.pos = list(
              position = ggpp::position_dodgenudge(x = .095, width = 0.105), width = 0.075))  +
  theme_classic() + theme(axis.text.x=element_blank(), axis.ticks.x=element_blank())+
  scale_fill_manual(values = c("#21908CFF", "#440154FF" )) +
  scale_color_manual(values = c("#21908CFF", "#440154FF"))+   theme( plot.title = element_text(size=11)) +
  ggtitle("") + ylab("Compartment Strength") + xlab("")


## Plot differntial dot plot
with(compared, plot(WT, Rad21, xlim = c(0,5), ylim = c(0,5), pch = 20))
abline(a = 0, b = 1, lty = 5)

## chr21:32.4-39Mb; [D], chr1:167-177Mb lets check these, from rao 2017
#lifted chr21	31027681	37627698	chr21:32400001-39000000	1
#lifted chr1	167030763	177030864	chr1:167000001-177000000	1

## Plot a test matrix ## Need higher resolution for these
hic_matrixplot(exp1 = synced.100kb[[1]],
               exp2 = synced.100kb[[2]],
               chrom = '21',
               start = 15000000,
               end=45000000,
               cut.off = 500, # upper limit of contacts
               skipAnn = F) # skip the outside annotation

gc()

## comp strength - down: 3, 5, 10, 13, 19, 20  up: 7,8, 11, 14, 17,18 (MICRO-C rad21)
## Pearson correlation matrix
## chr17.comp.matrix.plot.r21.rw.UP chr11.comp.matrix.plot.r21.corr.UP
compartment_matrixplot(
  exp1 = synced.100kb[[1]],
  CS_discovery = cs,
  chrom = "chr11", arm = "q",
  metric = "correlation"
)


## Test plots chr13.comp.matrix.plot.wt.rw.DN chr13.comp.matrix.plot.r21.corr.DN
compartment_matrixplot(
  exp1 = synced.100kb[[2]],
  CS_discovery = cs,
  chrom = "chr13", 
  arm = "q",
  metric = "correlation"
)

# chr14.comp.matrix.plot.wt
compartment_matrixplot(
  #exp1 = synced[[1]],
  exp1 =  synced[[2]],
  CS_discovery = cs,
  chrom = "chr14", arm = "q",
  #metric = "obsexp"
  colour_lim = c(0, 80)
  
)

# chr1.comp.matrix.plot.wt
compartment_matrixplot(
  exp1 = wt.100kb,
  #exp1 = rad21.100kb,
  CS_discovery = cs,
  chrom = "chr1", arm = "p",
  start = 167030763,
  end=177030864,
  #metric = "obsexp",
  #colour_bar = T
)


## Read in TAD data
wt.loops= read.delim('/projects/b1042/BackmanLab/HiC2/opt/juicer/work/022324_microC/contact_domains/CTRL/mega/loop_domains/loops.bedpe', h = F, skip = 1)
rad21.loops = read.delim('/projects/b1042/BackmanLab/HiC2/opt/juicer/work/022324_microC/contact_domains/6hrs_Aux/mega/loop_domains/loops.bedpe', h = F, skip = 1)

wt.loops <- wt.loops[-c(1),]
rad21.loops  <- rad21.loops[-c(1),]

wt.loops <- wt.loops[,c(1,22,23)]
rad21.loops  <- rad21.loops[,c(1,22,23)]

colnames(wt.loops) <- c("chrom", "start", "end")
colnames(rad21.loops) <- c("chrom", "start", "end")

hic_matrixplot(exp1 = synced[[1]],
               exp2 = synced[[2]],
               chrom = 'chr7',
               start = 25e6,
               end=30e6,
               loops = wt.loops, # see APA
               loops.colour = 'blue', # blue loops
               loops.type = 'upper', # only plot in upper triangle
               loops.radius = 20e3, # expand for visibility
               #tads = WT_TADs, # see ATA
               #tads.type = 'lower', # only plot in lower triangle
               #tads.colour = 'limegreen', # green TAD-borders
               cut.off = 25)

hic_matrixplot(exp1 = synced[[1]],
               exp2 = synced[[2]],
               chrom = '7',
               start = 25e6,
               end=30e6,
               loops = list(rad21.loops, wt.loops),  # see APA
               loops.colour = c('red','#91cf60'), # blue loops
               loops.type = c('lower','upper'), # only plot in upper triangle
               loops.radius = c(20e3, 20e3), # expand for visibility
               cut.off = 25, # upper limit of contacts
               skipAnn = T) 


hic_matrixplot(exp1 = synced[[1]],
               chrom = 'chr7',
               start = 25e6,
               end=30e6,
               tads = wt.loops, # see ATA
               tads.type = 'upper', # only plot in lower triangle
               tads.colour = '#91cf60', # green TAD-borders
               cut.off = 10, # upper limit of contacts
               skipAnn = T) # skip the outside annotation

hic_matrixplot(exp1 = synced[[1]],
               exp2 = synced[[2]],
               chrom = '7',
               start = 25e6,
               end=30e6,
               loops = actd.loops,  # see APA
               loops.colour = 'red', # blue loops
               loops.type = 'lower', # only plot in upper triangle
               #loops.radius = c(20e3, 20e3), # expand for visibility
               cut.off = 20, # upper limit of contacts
               skipAnn = T) 

APA<- APA(list("WT" = wt.5kb,'Rad21' = rad21.5kb),
          dist_thres = c(200e3, Inf),
          bedpe = wt.loops)

## Colors for difference plot
colors <- colorRampPalette(c("#443A83FF", "white","black"))(12)

## Colors for ATA part
v_colors <- c("white","#F4EDCA", "#EAD357FF", "#D0BE67FF", "#B6A971FF", "#9E9677FF", "#878479FF", "#727374FF", "#5E626EFF", "#48526BFF", "#2A406CFF", "#00306FFF", "#00204DFF")


## https://github.com/robinweide/GENOVA/issues/298
visualise(APA, colour_lim = c(0,30),metric = "diff", focus = 2) +
  ggplot2::scale_fill_gradientn(
    colours = c(v_colors),limits= c(0,20), na.value="black"
  ) +
  ggplot2::scale_colour_gradientn(limits = c(-5,5),
                                  aesthetics = "altfill",
                                  colours = c(colors),na.value="#365C8DFF",
                                  guide = ggplot2::guide_colourbar(available_aes = "altfill")
  )


## ------------------------------- ##
## Alt analysis of one Rep         ##
## ------------------------------- ##


rm(list = ls())

# import packages
library(ggplot2)
library(dplyr)
library(tidyverse)
library(hrbrthemes)
library(forcats)
library(GENOVA)
library(strawr)


## directories and variables 
dir <- "/projects/b1042/BackmanLab/HiC2/opt/juicer/work/112123_HiC/juicer_analysis/"
list.files(paste0(dir), pattern=".hic", all.files=F, include.dirs = FALSE)

exp.wt = "WT_HCT116_CTRL"
exp.r21 = "Rad21_6hrs_Aux"
hic ="inter_30.hic"

## -- introduce data -- ##

## contact probability resolution - 5kb
wt.5kb <- load_contacts(signal_path = paste0(dir,exp.wt,"/Rep2/aligned/",hic),
                        sample_name = "WT",
                        resolution = 5e3,
                        balancing = 'KR', # this is the default
                        colour = "black")

rad21.5kb <- load_contacts(signal_path = paste0(dir,exp.r21,"/Rep2/aligned/",hic),
                           sample_name = "Rad21",
                           resolution = 5e3,
                           balancing = 'KR', # this is the default
                           colour = "black")

## contact probability resolution - 10 kb
wt_1 <- load_contacts(signal_path = paste0("/projects/b1042/BackmanLab/juicer/work/112123_HiC/juicer_analysis/WT_HCT116_CTRL/Rep1/aligned/",hic),
                      sample_name = "WT_1",
                      resolution = 10e3,
                      balancing = 'KR', # this is the default
                      colour = "black")
gc()
wt_2 <- load_contacts(signal_path = paste0("/projects/b1042/BackmanLab/juicer/work/112123_HiC/juicer_analysis/WT_HCT116_CTRL/Rep2/aligned/",hic),
                      sample_name = "WT_2",
                      resolution = 10e3,
                      balancing = 'KR', # this is the default
                      colour = "black")

rad21_2 <- load_contacts(signal_path = paste0(dir,exp.r21,"/Rep2/aligned/",hic),
                         sample_name = "Rad21",
                         resolution = 10e3,
                         balancing = 'KR', # this is the default
                         colour = "black")

rad21.10kb.1 <- load_contacts(signal_path = paste0(dir,exp.r21,"/Rep1/aligned/",hic),
                              sample_name = "Rad21",
                              resolution = 10e3,
                              balancing = 'KR', # this is the default
                              colour = "black")

## TAD resolution - 25kb
wt.25kb <- load_contacts(signal_path = paste0(dir,exp.wt,"/Rep2/aligned/",hic),
                         sample_name = "WT",
                         resolution = 25e3,
                         balancing = 'KR', # this is the default
                         colour = "black")

rad21.25kb <- load_contacts(signal_path = paste0(dir,exp.r21,"/Rep2/aligned/",hic),
                            sample_name = "Rad21",
                            resolution = 25e3,
                            balancing = 'KR', # this is the default
                            colour = "black")

## Compartment resolution - 100kb
wt.100kb <- load_contacts(signal_path = paste0(dir,exp.wt,"/Rep2/aligned/",hic),
                          sample_name = "WT",
                          resolution = 1e5,
                          balancing = 'KR', # this is the default
                          colour = "black")

rad21.100kb <- load_contacts(signal_path = paste0(dir,exp.r21,"/Rep2/aligned/",hic),
                             sample_name = "Rad21",
                             resolution = 1e5,
                             balancing = 'KR', # this is the default
                             colour = "black")

actd.1 <- load_contacts(signal_path = paste0(dir,exps[1],"/Rep1/aligned/",hic),
                        sample_name = "ActD_1",
                        resolution = 25e3,
                        balancing = 'KR', # this is the default
                        colour = "black")

actd.2 <- load_contacts(signal_path = paste0(dir,exps[1],"/Rep2/aligned/",hic),
                        sample_name = "ActD_2",
                        resolution = 25e3,
                        balancing = 'KR', # this is the default
                        colour = "black")
gc()
## -- Triangle plots -- ##

## These loci were all taken from Rao et al 2017 figure 1a

wt <- wt.5kb
aux <- rad21.5kb

wt <- wt.10kb
aux <- rad21.10kb

## Loci from chr8 133800000 134600000 (hg19) converted to chr8 132787754	133587757 (hg38)
pyramid(exp = wt_1,
        chrom = "8",
        start = 132787754,
        end=133587757)

pyramid(exp = aux,
        chrom = "8",
        start = 132787754,
        end=133587757)

## Loci from chr4 40800000 42100000 (hg19) converted to chr4 40797983	42097983 (hg38)
pyramid(exp = wt_1,
        chrom = "4",
        start = 40797983,
        end=42097983)

pyramid(exp = aux,
        chrom = "4",
        start = 40797983,
        end=42097983)

## Loci from chr1 91900000 95800000 (hg19) converted to chr1	91434443	95334444 (hg38)
pyramid(exp = wt_1,
        chrom = "1",
        start = 91434443,
        end=95334444)

pyramid(exp = aux,
        chrom = "1",
        start = 91434443,
        end=95334444)

## KO pyramid using first loci

library(patchwork)

wt.pyr <- pyramid(exp = wt,
                  chrom = "8",
                  start = 132787754,
                  end=133587757) + ggplot2::ggtitle("Wildtype")

ko.pyr <- pyramid(exp = aux,
                  chrom = "8",
                  start = 132787754,
                  end=133587757) + ggplot2::ggtitle("Rad21 KO")

ko.pyr / wt.pyr + plot_layout(guides = "collect")

## Difference pyramid
pyramid_difference(
  exp1 = synced[[1]],
  exp2 = synced[[2]],
  chrom = "8",
  start = 132787754,
  end=133587757)

## -- domainogram -- ##

## Check ideal window size for optimal tad calling

ID <- insulation_domainogram(
  wt,
  chrom = '7',
  start = 25e6,
  end = 29e6,
  window_range = c(1, 101),
  step = 2
)
visualise(ID)

## -- TAD analysis -- ##

## sync indices between conditions
synced <- sync_indices(list(wt,aux))

lapply(synced, attributes) ## Get all attributes from list
synced[[1]][1] ## Access first element in list of wt

rm(aux,wt) ## declutter workspace as needed

## Checking different insulation windows using synced list
insulation.w25 <- insulation_score(synced,window = 25)
insulation.w40 <- insulation_score(synced,window = 40)
insulation.w50 <- insulation_score(synced,window = 50) ## Window size 50 seems to be the most relevant when tested, for 5kb resolution
insulation.w60 <- insulation_score(synced,window = 60)

## Inspect insulation windows
visualise(insulation.w50,
          chr = '4', start = 91434443, end=95334444,
          contrast = 2)

## Check INF and NA values
check <- data.frame(insulation.w25[1])
sum(is.infinite(check$insula_score.WT))
sum(is.na(check$insula_score.WT))
nrow(check)

## Call TADs on different window sizes #himatrix.91.4-95.3.diff.5kb.50w
TADcalls <- call_TAD_insulation(insulation.w25)
TADcalls <- call_TAD_insulation(insulation.w40)
TADcalls <- call_TAD_insulation(insulation.w50) ## use window size 50
TADcalls <- call_TAD_insulation(insulation.w60) 

## Plot a test matrix
hic_matrixplot(exp1 = synced[[1]],
               exp2 = synced[[2]],
               chrom = '4',
               start = 91434443,
               end=95334444,
               tads = list(TADcalls$Rad21, TADcalls$WT), # see ATA
               tads.type = list('lower','upper'), # only plot in lower triangle
               tads.colour = c('red','#91cf60'), # green TAD-borders
               cut.off = 7.5, # upper limit of contacts
               skipAnn = T) # skip the outside annotation

## aggregate TAD analysis
ATA <- ATA(synced, bed = TADcalls$WT)
ATA

list.files(dir)
visualise(ATA,
          colour_lim = c(0,20),
          colour_lim_contrast = c(-5,5),
          metric = "diff",
          focus = 2
)


## -- TAD analysis - 2D Saddleplot -- ##
## Taken from https://github.com/robinweide/GENOVA/issues/272

library(data.table)
library(ggplot2)

## function written for 2 samples
saddle.ata <- function(ata) {
  
  # Melt array
  df <- data.table(
    x = as.vector(slice.index(ata$signal, 1)),
    y = as.vector(slice.index(ata$signal, 2)),
    sample = as.vector(slice.index(ata$signal, 3)),
    value = as.vector(ata$signal)
  )
  # Calculate distance (unit is arbitrary due to ATA)
  df[, dist := x - y]
  # Take averages per off-diagonal band
  df <- df[, list(value = mean(value)), by = c("dist", "sample")]
  # Split by sample
  df <- split(df, df$sample)
  
  # Adjust sample 2 a bit so that a difference will show
  # Don't do this with real data!
  #df[[2]][, value := sqrt(value)]
  
  # Recombine the two samples, if you have more, this needs to be repeated for all samples
  df <- df[[1]][df[[2]], on = c("dist")]
  # Convert distance to TAD units
  df[, dist := scales::rescale(dist, to = c(-2, 2))]
  
  return(df)
}

df.saddle <- saddle.ata(ATA)

## Prepare data for plotting
sample1 <- data.frame(dist = df$dist, cond = rep("wt", nrow(df)), value = df$value)
sample2 <- data.frame(dist = df$dist, cond = rep("rad21", nrow(df)), value = df$i.value)

plot.df <- rbind(sample1, sample2)

# Plot 2D saddle
p <- plot.df %>%
  ggplot(aes(x=dist , y=value)) + 
  geom_line(size = 0.5, aes(color = cond)) +theme_classic()+
  #scale_y_continuous(limits = c(-2.5,  2))+
  #scale_x_continuous(limits = c(5,8.5), breaks = c(4, 4.5,5, 5.5, 6,6.5,7,7.5, 8, 8.5))+
  scale_color_manual(values=c("#C70039", "black")) +
  labs(title="Aggregate TAD Analysis", subtitle="Wt vs Rad21", 
       caption="")+theme(text = element_text(size = 12),axis.text.x = element_text(angle = 0, vjust = 1, hjust=1, size = 12),
                         axis.text.y = element_text(angle = 0,size = 12))+ 
  xlab("Distance From Diagonal (AU)") + 
  ylab("Signal")
p


## -- ---------------------- -- ##

## -- Relative Contact Probability -- ##

rcp <- RCP(explist = list(actd.1, actd.2),
           chromsToUse = NULL)
gc()
rm(pol2.5kb)
saveRDS(rcp, "/home/lmc0633/rad21paper/test.plots/rcpwt.rds")

rcp <- RCP(explist = synced,
           chromsToUse = NULL)
## Plot RCP
visualise(rcp) +scale_color_manual(values=c("black" ,"#C70039"))
dev.off()
## PLot log2 differential of RCP between conditions
visualise(rcp, contrast = 1, metric = 'lfc')

## -- trans interactions -- ##

## Matrix plot of all trans interactions
cm <- chromosome_matrix(list(WT = wt.100kb, KO = rad21.100kb),include_chr = "all",sort_chr = TRUE)

labels <- as.character(c(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,"X"))
visualise(cm)
#visualise(cm) + scale_y_discrete(labels = labels)+ scale_x_discrete(labels = labels)


## -- Compartments -- ##

## Group compartments into a single list
synced.100kb <- sync_indices(list(wt.100kb, rad21.100kb))

## Get eigens
cs = compartment_score(synced.100kb)

## Lets check all chroms for which ones show significant effect
## 4, 5, 6, 11, 13(really good), 16(really good),17, 18, 20
visualise(cs, chr = "X") +scale_color_manual(values=c("black" ,"#C70039"))

## Saddle plot calculation
attr(cs, 'signed') <- F # Fake signed CS
saddle.plot = saddle(synced.100kb,
                     CS_discovery = cs,
                     bins = 60)

p <- visualise(saddle.plot)

## Remap for plotting with different color scale
p$layers[[1]]$mapping <- ggplot2::aes(fill = log2(obsexp))
lims <- range(log2(p$data$obsexp))
values <- scales::rescale(c(lims[1], 0, lims[2]))

## Plot saddle
p + ggplot2::scale_fill_gradientn(
  colours = c("blue", "white", "red"),
  values  = values,
  limits = lims,
  name = "Log2(FC)"
)

library(PupillometryR)

CSS <- quantify(saddle.plot)
compared <- tidyr::spread(unique(CSS[,-c(3,4)]), key = 'exp', value = 'strength')

plot.df <- compared %>% 
  gather(cond, value, 2:3) 
max(plot.df$value)
min(plot.df$value)
mean(plot.df$value)
plot.df%>%
  ggplot( aes(x=cond, y=value, fill=cond)) +
  geom_flat_violin(color = NA, position = position_nudge(x = .15), trim=F)+ coord_flip() +
  scale_fill_manual(values=c("#C70039", "black"))+
  stat_boxplot(geom = "errorbar", width = 0.15, aes(color = "grey")) +
  theme_classic() + geom_boxplot(width = 0.2,notch=TRUE, notchwidth = 0.3, outlier.shape = NA, color = "white")+
  scale_y_continuous(limits = c(0,5))+
  theme(legend.position="none", plot.title = element_text(size=11), axis.text.x = element_text(angle = 0)) +
  ggtitle("Compartment Insulation") + xlab("Condition") + ylab("Strength")

## Plot differntial dot plot
with(compared, plot(WT, Rad21, xlim = c(0,5), ylim = c(0,5), pch = 20))
abline(a = 0, b = 1, lty = 1)

## chr21:32.4-39Mb; [D], chr1:167-177Mb lets check these, from rao 2017
#lifted chr21	31027681	37627698	chr21:32400001-39000000	1
#lifted chr1	167030763	177030864	chr1:167000001-177000000	1

## Plot a test matrix ## Need higher resolution for these
hic_matrixplot(exp1 = synced.100kb[[1]],
               exp2 = synced.100kb[[2]],
               chrom = '21',
               start = 31027681,
               end=37627698,
               cut.off = 1500, # upper limit of contacts
               skipAnn = F) # skip the outside annotation

## Pearson correlation matrix
compartment_matrixplot(
  #exp1 = wt.100kb,
  exp1 = rad21.100kb,
  CS_discovery = cs,
  chrom = "13", arm = "all",
  metric = "correlation",
  colour_bar = T
)

## Test plots
compartment_matrixplot(
  #exp1 = wt.100kb,
  exp1 = rad21.100kb,
  CS_discovery = cs,
  chrom = "21", arm = "all",
  start = 31027681,
  end=37627698,
  metric = "correlation",
  colour_bar = T
)



compartment_matrixplot(
  exp1 = wt.100kb,
  #exp1 = rad21.100kb,
  CS_discovery = cs,
  chrom = "1", arm = "all",
  start = 167030763,
  end=177030864,
  metric = "correlation",
  colour_bar = T
)


rm(compartment_matrixplot)

## -- check reproducibility -- ##

## directories and variables 
dir.root <- "/projects/b1042/BackmanLab/juicer/work/112123_HiC/juicer_analysis/"
dir.hic <- "/aligned/"
reps <- c("Rep1", "Rep2")
exps <- c("1hr_ActD",  "Pol2_6hrs_Aux")

list.files(paste0(dir.root, exps[1],"/",reps[1], dir.hic), all.files=F, include.dirs = FALSE)

## -- check reproducibility -- ##
library(hicrep)
paste0(dir.root, exps[1],"/",reps[1], dir.hic,"inter_30.hic")
## Compute a single chromosome
hic.r1 <- paste0(dir.root, exps[1],"/",reps[1], dir.hic,"inter_30.hic")
hic.r2 <- paste0(dir.root, exps[1],"/",reps[2], dir.hic,"inter_30.hic")

hic.r1 <- paste0(dir.root, exps[2],"/",reps[1], dir.hic,"inter_30.hic")
hic.r2 <- paste0(dir.root, exps[2],"/",reps[2], dir.hic,"inter_30.hic")

mat1 <- hic2mat(hic.r1, chromosome1 = '18', chromosome2 = '18', resol =   25000, method = "NONE")
mat2 <- hic2mat(hic.r2, chromosome1 = '18', chromosome2 = '18', resol =   25000, method = "NONE")
gc()

scc.out = get.scc(mat1, mat2, resol = 25000, h = 5, lbr = 0, ubr = 5000000)
scc.out

## Compute all chromosomes
res = 10000

all.scc <- list()
for (i in paste0(c(as.character(1:22), "X"))){
  
  cat(paste0("processing chrom: ", i))
  
  mat1.chr <- hic2mat(hic.r1, chromosome1 = i, chromosome2 = i, resol =   res, method = "NONE")
  mat2.chr <- hic2mat(hic.r2, chromosome1 = i, chromosome2 = i, resol =   res, method = "NONE")
  all.scc[[i]] = get.scc(mat1.chr, mat2.chr, res, 5, lbr = 0, ubr = 5000000)
  
}

gc()

## extract values
scc.df <- data.frame()
for(i in c(1:length(all.scc))){
  
  scc.corr <- all.scc[[i]][3]
  name <- names(all.scc[i])
  std <- all.scc[[i]][4]
  
  df <- data.frame(name,scc.corr, std)
  scc.df <- rbind(scc.df, df)
}

## Add factors
scc.df$factor<- ifelse(scc.df$name == "X" , 23, scc.df$name)
scc.df$factor <- as.numeric(scc.df$factor)

library(viridis)

## dotplot with eror bars
plot <- scc.df %>% 
  mutate(name = fct_reorder(name, as.numeric(factor))) %>% 
  ggplot(aes(x=name, y=scc, color = scc)) + scale_color_viridis( option = "D")+
  geom_point(size=2) + theme_classic() +
  geom_errorbar(aes(ymin=scc-std, ymax=scc+std))

plot
