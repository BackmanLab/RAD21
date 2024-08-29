## Rad21 Micro-C Analysis

### Notes

**Built in normalizations**
micro C data has SCALE and NONE normalizations

**Available resolutions**
[2500000, 1000000, 500000, 250000, 100000, 50000, 25000, 10000, 5000, 2000, 1000, 500, 200, 100, 50, 20, 10]

### 022624: initial submission

 Notes on converting micro-c .hic data from Encode to .cool format, to use [cooltools](https://cooltools.readthedocs.io/en/latest/) and generate smaller file size. [GENOVA](https://rdrr.io/github/robinweide/GENOVA/) stalls when loading HiC data over 10 gb in size. Cooltools was used by this [group](https://epigeneticsandchromatin.biomedcentral.com/articles/10.1186/s13072-022-00473-4#Sec12) that did micro C for visialization. Cooler conversion also allows me to access specific resolutions in a smaller file. Data is below

 * [Micro-c Control data](https://www.encodeproject.org/experiments/ENCSR958BEA/)

* [Micro-c 6 hours aux](https://www.encodeproject.org/experiments/ENCSR087JOM/)

* [All experiments](https://www.encodeproject.org/search/?searchTerm=5-Phenyl-1H-indole-3-acetic+acid+rad21&type=Dataset&limit=50)

* [RNAseq](https://www.encodeproject.org/experiments/ENCSR094RQC/)

**Create an environment on Quest Analytics for running Jupiter Notebooks**

I created a conda environment using the guidance [here](https://kb.northwestern.edu/page.php?id=94116) so I could install [HiCStraw](https://pypi.org/project/hic-straw/) and [Cool Tools](https://cooltools.readthedocs.io/en/latest/) on NU's HPC for the conversion


```shell
module purge all
module load python-miniconda3/4.12.0
conda create --name jupyter-kernel-cooltools -c anaconda python=3.9 numpy pandas matplotlib ipykernel --yes
conda activate jupyter-kernel-cooltools
pip install hic2cool --user
pip install cooltools  --user
pip install cooler --user

```
Use conda in the future when installing for environment. Pip may bypass dependency check. Use pip command below if needed.

``` conda init ``` will rewrite your bashrc for you so you can access those applications

```shell pip install --no-deps <LIB>```

 to reactivate the environment and use these packages use the following commands:

 ```shell
module purge all
module load python-miniconda3/4.12.0
conda activate jupyter-kernel-cooltools
```

When finished, deactivate using the following command ```shell conda deactivate```. .hic files can be parsed in this environment using the following python code in Jupiter Notebooks using one of the Quest analytic nodes ordered through [Quest On Demand](https://services.northwestern.edu/TDClient/30/Portal/KB/ArticleDet?ID=2234), taken from two github issues sections, [one](https://github.com/deeptools/HiCExplorer/issues/821) and [two](https://github.com/deeptools/HiCExplorer/issues/798).

```python
## Import packages
import hicstraw
import numpy as np
import os
import pandas as pd
import cooler

## Get working directory
dir="/projects/b1042/BackmanLab/HiC2/opt/juicer/work/022324_microC/juicer_analysis/CTRL/mega/aligned/"
#dir="/projects/b1042/BackmanLab/HiC2/opt/juicer/work/022324_microC/juicer_analysis/6hrs_Aux/mega/aligned/"
os.chdir(dir)

## Set workign directory
cwd = os.getcwd()
print(cwd)

## Set variables
resolution = 100000
data_type = 'observed' # (previous default / "main" data) or 'oe' (observed/expected)
normalization = "NONE"  # , VC, VC_SQRT, KR, SCALE, etc.
cond = "CTRL"

hic_file = '/projects/b1042/BackmanLab/HiC2/opt/juicer/work/022324_microC/juicer_analysis/'+cond+'/mega/aligned/inter_30.hic'
cool_file = '/projects/b1042/BackmanLab/HiC2/opt/juicer/work/022324_microC/juicer_analysis/'+cond+'/mega/aligned/inter_30_'+str(resolution/1000)+'k.cool'

##  bring in HiC file using Straw
hic = hicstraw.HiCFile(hic_file)

## Set resolutions
assert resolution in hic.getResolutions(), \
    f"{resolution} is not part of the possible resolutions {','.join(hic.getResolutions())}"

## Declare chromsizes data frame and write the chromosome sizes
chrom_sizes = pd.Series({chrom.name: chrom.length for chrom in hic.getChromosomes() if chrom.name != "All"})
# First write the chromosome sizes:
with open(hic.getGenomeID() + '.size', 'w') as fsize:
    for chrom in hic.getChromosomes():
        if chrom.name != "All":
            fsize.write(f"{chrom.name}\t{chrom.length}\n")

## Then write the counts in text file:
with open(cool_file.replace('.cool', ".txt"), 'w') as fo:
    for i in range(len(chrom_sizes)):
        for j in range(i, len(chrom_sizes)):
            chrom1 = chrom_sizes.index[i]
            chrom2 = chrom_sizes.index[j]
            result = hicstraw.straw(data_type, normalization, hic_file, chrom1, chrom2, 'BP', resolution)
            for k in range(len(result)):
                start1 = result[k].binX
                start2 = result[k].binY
                value = result[k].counts
                fo.write(f"{chrom1}\t{start1}\t{start1}\t{chrom2}\t{start2}\t{start2}\t{value}\n")

## Call command for cooler to convert .txt to .hic
## This call does not work in the environment for some reason
os.system(f"cooler load -f bg2 {hic.getGenomeID()}.size:{resolution} {cool_file.replace('.cool', '.txt')} {cool_file}")
```

My /.bashrc file in cd~ is configured for use with python-miniconda3/4.12.0. Use of cooler only requires loading ```shell module load python-miniconda3/4.12.0``` This is the usage for running [cooler load](https://cooler.readthedocs.io/en/latest/cli.html#cooler-load) on the bedgraph txt file generated by the python code above

```shell
cooler load <options> <chromsizesfile>:<resolution> <bedgraphfile> <coolfilename>
```
Here's an example of the code I used

```shell
cooler load -f bg2 hg38.size:5000 inter_30_5000.txt inter_30_5k.cool
cooler load -f bg2 hg38.size:5000 inter_30_5.0k.txt inter_30_5k.cool
```

Individual replicates also need to be merged. To merge .cool reps that were generated, use the following command in the CLI

```shell 
cooler merge -c <chunksize> <mergedoutfile> <infilerep> 
```
an example: 

Use '::' to specify a group path", to the inside of the mcool file. can't directly merge the entire mcool file, but can merge replicates at a single resolution. 

```shell 
cooler merge -c 10000 MB_5kb_merged.cool path/inter_30_5k_r1.cool path/inter_30_5k_r2.cool
```

The above files can be run on Quest without memory allocation as long as the resolution is low (50k-100k). For higher resolution data, use a slurm script like below:

```shell
#!/bin/bash
#SBATCH -A p32171
#SBATCH -p short
#SBATCH -N 1
#SBATCH --ntasks-per-node=1
#SBATCH --mem=180G
#SBATCH --mail-user=lucascarter2025@u.northwestern.edu
#SBATCH --mail-type=BEGIN
#SBATCH --mail-type=END
#SBATCH -t 4:00:00
#SBATCH --job-name=run_cooler

# Load modules
module purge all
module load python-miniconda3/4.12.0

# Set your working directory
cd /projects/b1042/BackmanLab/HiC2/opt/juicer/work/022324_microC/juicer_analysis/CTRL/mega/aligned/

# Run python script to convert contacts to contact matrix
cooler load -f bg2 hg38.size:5000 inter_30_5.0k.txt inter_30_5k.cool
```

These coolers then need to be balanced before analysis. To balance .cool files that were generated, use the following command in the CLI

```shell 
cooler balance -p <nodes> -c <chunksize> <coolfile>
```

an example: 

```shell 
cooler balance -p 1 -c 10000 inter_30_100k.cool 
```

Once coolers are generated and balanced for each resolution of interest from .hic files, these can be analyzed in R using [GENOVA](https://rdrr.io/github/robinweide/GENOVA/).

color scales for GENOVA from [here](https://github.com/robinweide/GENOVA/issues/190)

```shell

inferno     -> bezier_corrected_hot
blue-red    -> bezier_corrected_divergent 
white-red   -> bezier_corrected_whiteRed 
ARA_diff    -> bezier_corrected_greenPink 

```
### GENOVA analysis in R is below

```R
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


## Rain plot of compartment strength
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


```

### 080924 RAD21 paper resubmission 

 We're using publicly available RAD21 AID2 Mint-ChIP data to compare TADs and loops for the distribution of these marks. All Mint-ChIP files can be found on encode [here](https://www.encodeproject.org/search/?searchTerm=5-Phenyl-1H-indole-3-acetic+acid+rad21&type=Experiment&target.label=H3K27ac&target.label=H3K9me3&target.label=H3K4me3&target.label=H3K27me3) 

* [h3k27ac control](https://www.encodeproject.org/experiments/ENCSR323XBZ/)

* [h3k27ac treated](https://www.encodeproject.org/experiments/ENCSR871TIY/)

* [h3k27me3 control](https://www.encodeproject.org/experiments/ENCSR899DBQ/)

* [h3k27me3 treated](https://www.encodeproject.org/experiments/ENCSR314MXA/)

* [h3k9me3 control](https://www.encodeproject.org/experiments/ENCSR412FVA/)

* [h3k9me3 treated](https://www.encodeproject.org/experiments/ENCSR808WHQ/)

* [h3k4me3 control](https://www.encodeproject.org/experiments/ENCSR175TRD/)

* [h3k4me3 treated](https://www.encodeproject.org/experiments/ENCSR356EWB/)

TF chip for CTCF:

* [ctcf control](https://www.encodeproject.org/experiments/ENCSR216SIK/)

* [ctcf treated](https://www.encodeproject.org/experiments/ENCSR107GWZ/)

CHIA PET

* [POLR2A treated](https://www.encodeproject.org/experiments/ENCSR771PKY/)

* [POLR2A control](https://www.encodeproject.org/experiments/ENCSR353XNY/)

micro c data:

 * [Micro-c Control data](https://www.encodeproject.org/experiments/ENCSR958BEA/)

* [Micro-c 6 hours aux](https://www.encodeproject.org/experiments/ENCSR087JOM/)

other data:

* [All experiments](https://www.encodeproject.org/search/?searchTerm=5-Phenyl-1H-indole-3-acetic+acid+rad21&type=Dataset&limit=50)

* [RNAseq](https://www.encodeproject.org/experiments/ENCSR094RQC/)

* [PROseq](https://www.encodeproject.org/gene-silencing-series/ENCSR968VYH/)


## Resubmission analysis in R

```R

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


```