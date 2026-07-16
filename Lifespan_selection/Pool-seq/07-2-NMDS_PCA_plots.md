# Make NMDS and PCA pots from population allele frequencies

# 2021-04-01

I will use a previous R script developed for this purpose, and made a few modifications so it receives parameters from the command line.

To reformat plots, I will download all the generated files:
```bash
cat color_by_spawning-season_15pools_NWHer.txt

SampleName	CategoryLabel	CategoryColor	CategoryShape
EI-F	Early-I-F	#fff7bc	0
EI-M	Early-I-M	#fee391	0
EII-F	Early-II-F	#fec44f	1
EII-M	Early-II-M	#fe9929	1
EIII-F	Early-III-F	#ec7014	2
EIII-M	Early-III-M	#cc4c02	2
EIV-F	Early-IV-F	#993404	5
EIV-M	Early-IV-M	#662506	5
LV-F	Late-LV-F	#f7fcb9	15
LV-M	Late-LV-M	#d9f0a3	15
LVI-F	Late-LVI-F	#addd8e	16
LVI-M	Late-LVI-M	#78c679	16
LVII-F	Late-LVII-F	#41ab5d	17
LVII-M	Late-LVII-M	#238443	17
LVIII-F	Late-LVIII-F	#006837	18
LVIII-M	Late-LVIII-M	#004529	18
```

I modified the R script to read in parameters from the command line, the script is called `07-02-NMDS-PCA-plots.sh`:

```bash
#!/bin/bash
#SBATCH -A def-ruzza
#SBATCH --cpus-per-task=4
#SBATCH --mem-per-cpu=8G
#SBATCH -t 20:00:00
#SBATCH -J NMDS-PCA-plots
#SBATCH -e NMDS-PCA-plots_%J_%A_%a.err
#SBATCH -o NMDS-PCA-plots_%J_%A_%a.out
#SBATCH --mail-type=all
#SBATCH --mail-user=angela.fuentespardo@gmail.com

# Load required software.
module load gcc/7.3.0 r/3.6.1

# Set environment variables.
# Files.
WORK_DIR='/home/afuentes/project/afuentes/analysis/07-popgen-stats/NMDS-PCA/'
INPUT_FILE='/Dropbox/PostDoc_UU/Projects/Beetles/analysis/06-apply-Neff-calc-AF/OTT_16_pools.UG.SNPs.hf.GQ20.mono.miss20.mac3.DPq5-99.AD.Neff.AF.txt'
COLORS_FILE='/home/afuentes/project/afuentes/analysis/07-popgen-stats/NMDS-PCA/color_by_spawning-season_15pools_NWHer.txt'
# Programs.
R_CODE='/home/afuentes/project/afuentes/code/utility-code/07-2-make-NMDS-PCA-plots.R'

# Run the R script, passing on the working directory and input file path (IN THIS ORDER!!)
Rscript $R_CODE $WORK_DIR $INPUT_FILE $COLORS_FILE

```


prcomp(P, center=FALSE, scale=FALSE) or, equivalently, princomp(P, cor=FALSE)


OLD --------
R script `07-2-make-NMDS-PCA-plots.R`:
```R
#####################################################################
#                                                                   #
# Perform Non-metric Multidimensional Scaling based on              #
# population allele frequencies of pool-seq data                    #
#                                                                   #
# E-mail: apfuentesp@gmail.com, Uppsala University                  #
# Date: 2020-04-02                                                  #
#                                                                   #
#####################################################################

## Clean environment space.
rm(list=ls())

## Receive arguments from the bash command.
#args <- commandArgs(trailingOnly=TRUE)
#cat(args, sep='\n')

## Assign these arguments to variables.
#working_dir <- args[1]  # /home/afuentes/project/afuentes/analysis/07-popgen-stats/nMDS
#input_file <- args[2]   # /home/afuentes/project/afuentes/06-apply-Neff-calc-AF/NWAtlantic_herring_15_pools.SNPs.hf.DP20-300.GQ10.mono.miss20.maf0.05.AD.Neff.AF.tp
#colors_file <- args[3]   # /home/afuentes/project/afuentes/analysis/07-popgen-stats/nMDS/color_by_spawning-season_15pools_NWHer.txt

## For debugging.
working_dir <- '~/Dropbox/PhD_Dal/2019-11_Reanalysis/analysis/07-1-popgen-stats/nMDS-PCA-plots'
#input_file <- '~/Dropbox/PhD_Dal/2019-11_Reanalysis/analysis/06-apply-Neff-calc-AF/test_NWAtlantic_herring_15_pools.SNPs.hf.DP20-300.GQ10.mono.miss20.maf0.05.AD.Neff.AF.100K.tp'
#input_file <- '~/Dropbox/PhD_Dal/2019-11_Reanalysis/analysis/06-apply-Neff-calc-AF/NWAtlantic_herring_15_pools.SNPs.hf.DP20-300.GQ10.mono.miss20.maf0.05.AD.Neff.AF.tp'
#input_file <- '~/Dropbox/PhD_Dal/2019-11_Reanalysis/analysis/06-apply-Neff-calc-AF/NWAtlantic_herring_13_pools.SNPs.hf.DP20-300.GQ10.mono.miss20.maf0.05.AD.Neff.AF.tp'
input_file <- '~/Dropbox/PhD_Dal/2019-11_Reanalysis/analysis/06-apply-Neff-calc-AF/NWAtlantic_herring_13_pools.SNPs.hf.DP20-300.GQ10.mono.miss20.maf0.05.AD.Neff.AF.txt'
#input_file <- '~/Dropbox/PhD_Dal/2019-11_Reanalysis/analysis/06-apply-Neff-calc-AF/NWAtlantic_herring_13_pools.SNPs.hf.DP20-300.GQ10.mono.miss20.mac4.AD.Neff.AF.txt'
#colors_file <- '~/Dropbox/PhD_Dal/2019-11_Reanalysis/analysis/07-1-popgen-stats/nMDS/color_by_spawning-season_15pools_NWHer.txt'
colors_file <- '~/Dropbox/PhD_Dal/2019-11_Reanalysis/analysis/07-1-popgen-stats/nMDS-PCA-plots/color_by_spawning-season_13pools_NWHer.txt'

print('# ------------- Set R objects from BASH command line: ------------- #')
working_dir
input_file
colors_file

## Set working directory.
setwd(working_dir)

## ------------------------------------------------------------------
## Data loading and preprocessing
## ------------------------------------------------------------------

## Load the pool allele frequency data.
path_to_file <- input_file

## Save filename.
library(tools)
filename <- file_path_sans_ext(basename(path_to_file))
filename

library(data.table)
poolData <- fread(path_to_file, data.table=FALSE, header=TRUE, sep='\t', stringsAsFactors=FALSE)

head(poolData)  # Verify modification.
dim(poolData)
str(poolData)
class(poolData)


## ------------------------------------------------------------------
## Generate nMDS plot
## ------------------------------------------------------------------

library("vegan")

## Remove rows (loci) with missing data,as NMDS cannot handle those.
sum(is.na(poolData))  # Count number of NAs
#[1] 74  # mac 4
#[1] 86  # maf 0.05
poolData_subset <- na.omit(poolData)
head(poolData_subset)
dim(poolData_subset)
class(poolData_subset)
#rm(poolData)  # Clean up space.

## Standardize allele frequency (AF) to the most common allele of one pop for visualization. I chose 1st pop for convenience.
## If AF of the first pop < 0.5, 1-AF for all columns in the row, if not, leave the AF as they are, check next row...
AF_df <- poolData_subset
head(AF_df)
rownames(AF_df) <- AF_df[ , 1]
AF_df <- AF_df[ , -1]
head(AF_df)
rm(poolData_subset)

AF_df[1:ncol(AF_df)] <- lapply(AF_df[1:ncol(AF_df)], function(x) ifelse(AF_df[,4] < 0.5, 1-x, x)) # AF_df[,4] AF will be standardized respect to Col 4.
AF_df_t <- t(AF_df)
dim(AF_df_t)
AF_df_t[,1:10]
rm(AF_df)  # Clean up space.

## Generate the NMDS object (excluding missing data).
set.seed(2)
NMDS <- metaMDS(AF_df_t, k=2) # k: The number of reduced dimensions. Rows: samples, Columns: loci.

class(NMDS)
str(NMDS)

## Obtain a stressplot to assess how good...# Stress < ?? good!
stressplot(NMDS)

## Load file listing colors por samples.
col_pop <- read.table(colors_file, header=TRUE, stringsAsFactors=FALSE)

print(col_pop)
class(col_pop)
target_order <- col_pop$Pop
target_order

NMDS_mt <- NMDS$points[match(target_order, row.names(NMDS$points)), ]
class(NMDS_mt)

## Make plot.
set.seed(77)

#pdf("file_name.pdf",width=10, height=10)
#plot(NMDS_mt, type="n", xlim = c(-0.03, 0.03), ylim = c(-0.0135,0.0125), # mac 4
plot(NMDS_mt, type="n", xlim = c(-0.02, 0.02), ylim = c(-0.02,0.02), # maf 0.05
     ylab="NMDS2", xlab="NMDS1",
     main=paste0("nMDS plot based on population allele frequencies \n File:",filename))
points(NMDS_mt, pch=19, col=col_pop$Color)
text(NMDS_mt[,1], NMDS_mt[,2], labels=row.names(NMDS_mt), pos=4, cex=0.7)
legend('topright', legend=c("Spring", "Fall"), fill=c("blue", "red"),
       box.col="black", cex=0.8, title="Spawning season")
abline(h=0, v=0, col="black")

#dev.off()

# http://huboqiang.cn/2016/03/03/RscatterPlotPCA


## ------------------------------------------------------------------
## Generate a PCA plot
## ------------------------------------------------------------------

## If there is missing data, it should be replaced by the mean allele frequency using scaleGen:
#x.lat <- tab(x, freq=TRUE, NA.method="mean")  # Following http://adegenet.r-forge.r-project.org/files/PRstats/practical-MVAintro.1.0.pdf

## Perform PCA.
library("adegenet")

pca1 <- dudi.pca(AF_df_t, center=TRUE, scale=FALSE, scannf = FALSE, nf = 12)
pca1

# http://huboqiang.cn/2016/03/03/RscatterPlotPCA

## Generate the plot.
library(ggplot2)
library(ggrepel)

percentage <- round(pca1$eig / sum(pca1$eig) * 100, 2)
percentage <- paste("PC",seq(1:length(pca1$eig))," (",paste(as.character(percentage),"% of genetic variation",") ",sep=""),sep="")


## Make all plots.
for (i in seq(1:(length(pca1$eig)-1))) {
  #i=1
  p <- ggplot(pca1$li, aes_string(x=paste0("Axis",i), y=paste0("Axis",i+1))) +
    geom_point(size=3, alpha = 1/1.4, aes(color=as.factor(col_pop$Color))) + #scale_size_manual(values=c(3,3)) + #alpha = 1/2
    #scale_shape_manual(values=c(16,17)) +
    #scale_fill_discrete(name = "Collection season", labels = c("Spring", "Fall")) +
    geom_text_repel(label=rownames(pca1$li)) +
    theme(panel.background = element_blank(),
          panel.border=element_rect(fill=NA),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          strip.background=element_blank(),
          axis.text.x=element_text(colour="black"),
          axis.text.y=element_text(colour="black"),
          axis.ticks=element_line(colour="black"),
          plot.margin=unit(c(1,1,1,1),"line")) +
    xlab(percentage[i]) + ylab(percentage[i+1]) +
    scale_color_identity(name = "Collection season",
                         breaks = c("blue","red"),
                         labels = c("Spring", "Fall"),
                         guide = "legend") +
    theme(legend.position = "bottom")
    #theme(legend.position = c(0.9, 0.08), legend.background = element_rect(#fill="lightblue",
    #  size=0.3, linetype="solid",
    #  colour ="black"))

  print(p)
}


## Make only plot PC1 and PC3.
ggplot(pca1$li, aes(x=Axis1, y=Axis3)) +
  geom_point(size=3, alpha = 1/1.4, aes(color=as.factor(col_pop$Color))) + #scale_size_manual(values=c(3,3)) + #alpha = 1/2
  #scale_shape_manual(values=c(16,17)) +
  #scale_fill_discrete(name = "Collection season", labels = c("Spring", "Fall")) +
  geom_text_repel(label=rownames(pca1$li)) +
  theme(panel.background = element_blank(),
        panel.border=element_rect(fill=NA),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.background=element_blank(),
        axis.text.x=element_text(colour="black"),
        axis.text.y=element_text(colour="black"),
        axis.ticks=element_line(colour="black"),
        plot.margin=unit(c(1,1,1,1),"line")) +
  xlab(percentage[1]) + ylab(percentage[3]) +
  scale_color_identity(name = "Collection season",
                       breaks = c("blue","red"),
                       labels = c("Spring", "Fall"),
                       guide = "legend") +
  theme(legend.position = "bottom")
  #theme(legend.position = c(0.9, 0.08), legend.background = element_rect(#fill="lightblue",
  #                                                                       size=0.3, linetype="solid",
  #                                                                       colour ="black"))

#ggsave("PCA_LAT_baseline_pure_Northern-Southern.pdf", scale=1)


## Make plot only for fall spawners.
to_remove <- c("BDO-S","NTS-S","SIL-S","SPH-S")
pca1_li_onlyFall <- pca1$li[!rownames(pca1$li) %in% to_remove, ]
pca1_li_onlyFall

color_fall <- rep("red", nrow(pca1_li_onlyFall))

ggplot(pca1_li_onlyFall, aes(x=Axis1, y=Axis2)) +
  geom_point(size=3, alpha = 1/1.4, aes(color=as.factor(color_fall))) + #scale_size_manual(values=c(3,3)) + #alpha = 1/2
  #scale_shape_manual(values=c(16,17)) +
  #scale_fill_discrete(name = "Collection season", labels = c("Spring", "Fall")) +
  geom_text_repel(label=rownames(pca1_li_onlyFall)) +
  theme(panel.background = element_blank(),
        panel.border=element_rect(fill=NA),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.background=element_blank(),
        axis.text.x=element_text(colour="black"),
        axis.text.y=element_text(colour="black"),
        axis.ticks=element_line(colour="black"),
        plot.margin=unit(c(1,1,1,1),"line")) +
  xlab(percentage[1]) + ylab(percentage[2]) +
  scale_color_identity(name = "Collection season",
                       breaks = c("blue","red"),
                       labels = c("Spring", "Fall"),
                       guide = "legend") +
  theme(legend.position = "bottom")
  #theme(legend.position = c(0.9, 0.08), legend.background = element_rect(#fill="lightblue",
  #  size=0.3, linetype="solid",
  #  colour ="black"))

#ggsave("PCA_LAT_baseline_pure_Northern-Southern.pdf", scale=1)

```
