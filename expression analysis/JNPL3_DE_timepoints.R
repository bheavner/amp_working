# Goal: DE analysis of JNPL3 Tau mice (transgenic (+) vs wild type (-) ) at each timepoint.

# I think this means building a dispersion model with all samples, then applying that model to DE testing for pairwise comparisons at each timepoint (similar to limma "makecontrasts"). I'll use timepoints as a factor variable (t = 2, 6, or 12 mo). As a future refinement, I'd like to use gender as a block effect for the rTG samples.

# output should be a table of top most differentially expressed genes, with the following columns: gene ID; logFc' logCPM; pValue' FDR' external_gene_name; mean gene level (tg); mean gene level (wt);  sd (tg); sd (wt); kurtosis (tg); kurtosis (wt); (I decided not to include 95% CI of effect size or SD of effect size - see https://support.bioconductor.org/p/61640/); histogram of top 10

## load libraries
library(dplyr) # for subsetting data
library(R.utils) # for unzipping data
library(synapseClient) # to download data
library(edgeR) # for DE analysis
library(biomaRt) # for gene name lookups
library(ggplot2) # for better boxplots
library(reshape2) # for melt for ggplot
library(PerformanceAnalytics) # for kertosis calculations

## Get data and define groups

#Login to Synapse using credentials saved in .synapseConfig file
synapseLogin()

# get the transposed readcount file and covariates file from synapse
countFile <- synGet('syn3192634') 
covariatesFile <- synGet('syn2875343') # the working dir copy

# unzip count file and load for processing
localCountFilePath <- getFileLocation(countFile)

if(!file.exists(sub('.gz', '', localCountFilePath))) {
  gunzip(localCountFilePath)
}

localCountFilePath <- sub('.gz', '', localCountFilePath) #trim the .gz suffix

counts <- read.table(localCountFilePath, header = TRUE, stringsAsFactors = FALSE)

#PROBLEM: - an rTGMinus sample
#in count data, sample “LP62_4”
#in covariates, sample “LP_62_4”
# for now, change column name in counts - TODO: fix covariates file
colnames(counts)[1] <- "LP_62_4"

# load covariates file to have handy
covariates <- read.table(getFileLocation(covariatesFile), header = TRUE, stringsAsFactors = FALSE)

# get JNPL3 subset of counts from counts file
JNPL3Counts <- dplyr::select(counts, 
                             one_of(covariates[(covariates$Experiment == "MAPT_P301L"), ]$Mouse_ID))

# define groups
targets <- data.frame(Sample = covariates[(covariates$Experiment == "MAPT_P301L"), ]$Mouse_ID, 
                      Genotype = covariates[(covariates$Experiment == "MAPT_P301L"), ]$Genotype, 
                      Age = as.factor(covariates[(covariates$Experiment == "MAPT_P301L"), ]$Age_months))

group <- factor(paste(targets$Genotype, targets$Age, sep = "."))

## Make DGEList object to start working on DE analysis
JNPL3 <- DGEList(counts = JNPL3Counts, group = group)

## filter data 
# require minimum of 100 counts per million for at least 2 samples
dim(JNPL3)
#39179    24

keep <- rowSums(cpm(JNPL3)>100) >= 2
d <- JNPL3[keep,]
dim(d)
#3133   24

# reset library sizes after filtering
d$samples$lib.size <- colSums(d$counts)

# normalize the data using TMM
d <- calcNormFactors(d, method = "TMM")

plotMDS(d, method="bcv", col=as.numeric(d$samples$group))
legend("bottomleft", as.character(unique(d$samples$group)), col=1:3, pch=20)
# note: now the clustering is more interesting: -2 cluster well and separate from all others, 2 of the +2 do (and one doesn't), one of the -6 is far away, the -12 cluster well, etc...

# GLM estimate of dispersion

# consider all levels of time for WT and TG separately: NPL3+ (transgenic), JNPL3- (WT), timepoints 2, 6, and 12 months
# see edgeR user guide p. 25 & 26 - I'm doing nested interaction approach

design.mat <- model.matrix(~0 + group) 
## WHY DOESN'T THIS MAKE targets$Genotype-:targets$Ageq2 and targets$Genotype+:targets$Age2 ?
# not needed - per seth "the contrast for the main effect of genotype (combined data from all time points) would be specified as (2+ + 6+ + 12+)/3 - (2- + 6- + 12-)/3"


#design.mat <- model.matrix(~0 + d$samples$group)
colnames(design.mat) <- levels(d$samples$group)
d2 <- estimateGLMCommonDisp(d,design.mat)
d2 <- estimateGLMTrendedDisp(d2,design.mat) #method = "power" also a possibility
# You can change method to "auto", "bin.spline", "power", "spline", "bin.loess".
# The default is "auto" which chooses "bin.spline" when > 200 tags and "power" otherwise.
d2 <- estimateGLMTagwiseDisp(d2, design.mat)
plotBCV(d2) # seems nice - blue line looks like it follows points...

# glm EDGER ANALYSIS WITH glmLRT TEST
fit <- glmFit(d2, design.mat) 
# Andrew used fit <- glmFit(d2, design.mat, dispersion=dge$tagwise.dispersion) but since "If NULL will be extracted from y, with order of precedence: tagwise dispersion, trended dispersions, common dispersion.", I'll leave it out and let it get it from d2.
de_glm  <- glmLRT(fit, contrast=c(-1, 0, 0, 1 ,0, 0)) # -1*JNPL3Minus.12mo 1*JNPL3Plus.12mo 
#so a positive logFC means transgenic (JNPL3Plus) is higher expression than WT (JNPL3Minus).

de_glm  <- glmLRT(fit, contrast=c(-1, -1, -1, 1, 1, 1)) # -1*JNPL3Minus 1*JNPL3Plus 

summary(decideTestsDGE(de_glm, p.value=0.01))
#[,1]
#-1    1
#0  3132
#1     0

topTags(de_glm)
# Look at DE with glmFit.
#Coefficient:  -1*-.12 -1*-.2 -1*-.6 1*+.12 1*+.2 1*+.6 
#logFC    logCPM        LR        PValue           FDR
#ENSMUSG00000060938  13.386082  6.756046 659.54313 1.877984e-145 5.883722e-142
#ENSMUSG00000062353 -31.111692  5.914889 167.17428  3.064507e-38  4.800550e-35
#ENSMUSG00000023087   2.200460  8.354657  77.98370  1.038945e-18  1.085005e-15
#ENSMUSG00000021091   4.490502  8.823266  64.27616  1.081463e-15  8.470558e-13
#ENSMUSG00000029816  11.811166  6.650852  58.55945  1.972398e-14  1.235905e-11
#ENSMUSG00000073411   3.492203  7.965863  52.16643  5.099034e-13  2.662546e-10
#ENSMUSG00000027014  -4.150732  6.378709  49.69279  1.798059e-12  8.047598e-10
#ENSMUSG00000036887   4.411974  7.506238  46.26312  1.033917e-11  4.049076e-09
#ENSMUSG00000068099   1.406636  7.707540  45.16059  1.815215e-11  6.318963e-09
#ENSMUSG00000079037   2.930920 10.942642  41.45994  1.203115e-10  3.769360e-08

# box plot of top differentially expressed genes across samples
topTen <- rownames(topTags(de_glm, n = 10))
toPlot <- d[topTen] #10 rows, 24 columns
top <- topTags(de_glm, n=50)

# get gene names for x axis labels
ensembl=useMart("ensembl", dataset="mmusculus_gene_ensembl")
geneNames <- getBM(c("ensembl_gene_id", "external_gene_name"), 
                   filters = "ensembl_gene_id", 
                   values = rownames(top)[1:10], 
                   ensembl)

# get label order right for plot
labels <- geneNames$external_gene_name[(order(match(geneNames$ensembl_gene_id, rownames(top)[1:10])))]

# munge to use with ggplot and plot
test  <- data.frame(toPlot$samples$group, t(toPlot$counts))
toPlot2 <- melt(test)
ggplot(toPlot2, aes(factor(variable), value)) + 
  geom_boxplot(aes(fill = toPlot.samples.group)) +
  scale_x_discrete(labels = labels)

## FINALLY, Build summary table with desired outputs (using the glmLRT results here)
top <- data.frame(topTags(de_glm, n=99)) # the 100th isn't in ensembl, so fix later..

ensembl=useMart("ensembl", dataset="mmusculus_gene_ensembl")
geneNames <- getBM(c("ensembl_gene_id", "external_gene_name"), 
                   filters = "ensembl_gene_id", 
                   values = rownames(top), 
                   ensembl)

# get label order right
names <- geneNames$external_gene_name[(order(match(geneNames$ensembl_gene_id, rownames(top))))]

top <- cbind(top, names)
colnames(top)[6] <- "Gene_Name"

## add mean gene level, SD, kurtosis, 95% CI of effect size or SD of effect size to summary table

means <- rowMeans(d2$counts[rownames(top), ])
top <- cbind(top, means)
colnames(top)[7] <- "Mean_count"

stdevs <- apply(d2$counts[rownames(top), ], 1, sd)
top <- cbind(top, stdevs)
colnames(top)[8] <- "Count_Std_Dev"

kurtosises  <- apply(d2$counts[rownames(top), ], 1, kurtosis)
top <- cbind(top, kurtosises)
colnames(top)[9] <- "Count_Kurtosis"

write.table(top, file="JNPL3_tg_vs_wt_timepoints.tsv", quote=F)