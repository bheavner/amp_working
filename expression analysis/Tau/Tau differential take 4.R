# perhaps normalization isn't needed for DE analysis - that's for quantification? I remain unclear on the point.

# meanwhile, follow this: https://web.stanford.edu/class/bios221/labs/rnaseq/lab_4_rnaseq.html

# Want to do DE analysis within tissue of Tau mice (JNPL3 and rTG4510), focusing on TG vs not TG.
# output should be a table of top 100 most differentially expressed genes, with the following columns: gene ID; logFc' logCPM; pValue' FDR' external_gene_name; mean gene level; sd; kurtosis; 95% CI of eeffect size or SD of effect size; histogram of valuse

# Next refinement -want covariates: age, sex, RIN, flow cell

library(dplyr) # for subsetting data
library(R.utils) # for unzipping data
library(synapseClient) # to download data
library(edgeR) # for DE analysis
library(biomaRt) # for gene name lookups
library(ggplot2) # for better boxplots
library(reshape2) # for melt for ggplot
library(PerformanceAnalytics) # for kertosis calculations

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

## define groups: JNPL3+ (transgenic), JNPL3- (WT), rTG+ rTG-

# JNPL3+:
JNPL3Plus <- dplyr::select(filter(covariates, Experiment == "MAPT_P301L" & Genotype == "+"), Mouse_ID)
# 15 samples

# JNPL3-:
JNPL3Minus <- dplyr::select(filter(covariates, Experiment == "MAPT_P301L" & Genotype == "-"), Mouse_ID)
# 9 samples

# rTG+:
rTGPlus <- dplyr::select(filter(covariates, Experiment == "rTG4510" & Genotype == "+"), Mouse_ID)
# 18 samples

# rTG-:
rTGMinus <- dplyr::select(filter(covariates, Experiment == "rTG4510" & Genotype == "-"), Mouse_ID)
# 18 samples

JNPL3Plus <- as.character(JNPL3Plus$Mouse_ID)
JNPL3Minus <- as.character(JNPL3Minus$Mouse_ID)

rTGPlus <- as.character(rTGPlus$Mouse_ID)
rTGMinus <- as.character(rTGMinus$Mouse_ID)

JNPL3PlusCols <- is.element(as.character(colnames(counts)), JNPL3Plus)
JNPL3MinusCols <- is.element(as.character(colnames(counts)), JNPL3Minus)
rTGPlusCols <- is.element(as.character(colnames(counts)), rTGPlus)
rTGMinusCols <- is.element(as.character(colnames(counts)), rTGMinus)

groups <- c(rep(0, length(counts[1,])))
groups[JNPL3PlusCols] <- "JNPL3Plus"
groups[JNPL3MinusCols] <- "JNPL3Minus"
groups[rTGPlusCols] <- "rTGPlus"
groups[rTGMinusCols] <- "rTGMinus"

#PROBLEM: - an rTGMinus sample
#in count data, sample “LP62_4”
#in covariates, sample “LP_62_4”
# -- now fixed on line 38
#rTGMinusCols[1] <- TRUE
#groups[1] <- "rTGMinus"

## now make DGEList object for edgeR
d <- DGEList(counts, group=factor(groups))

## filter data - require minimum of 100 counts per million for at least 2 samples
d.full <- d # keep the old one in case we mess up
dim(d.full)
#39179    60

keep <- rowSums(cpm(d)>100) >= 2
d <- d[keep,]
dim(d)
#3462   60

# reset library sizes after filtering
d$samples$lib.size <- colSums(d$counts)

# normalize the data using TMM
d <- calcNormFactors(d)


## Data Exploration
plotMDS(d, method="bcv", col=as.numeric(d$samples$group))
legend("bottomleft", as.character(unique(d$samples$group)), col=1:3, pch=20)

# some observations - X374367 and X307983 don't cluster with anything (this is different than before). 
# JNPL3Minus seems to separate into two groups: X370764, X370767, X370773, X370763 vs others. 
# This plot is worth spending more time with.

## Estimate dispersion. Since this is trying to get an idea of overall variability across the genome for the dataset, I'll go back, subset counts by tissue (JNPL3 = spinal, rTG = forebrain), and renormalize, etc.

JNPL3Samples <- c(JNPL3Plus, JNPL3Minus)
JNPL3Cols <- is.element(as.character(colnames(counts)), JNPL3Samples)
JNPL3 <- DGEList(counts = dplyr::select(counts, one_of(JNPL3Samples)),
                 group = groups[JNPL3Cols])

rTGSamples <- c(rTGPlus, rTGMinus)
rTGCols <- is.element(as.character(colnames(counts)), rTGSamples)

#PROBLEM: - an rTGMinus sample
#in count data, sample “LP62_4”
#in covariates, sample “LP_62_4”
# -- now fixed on line 38
#rTGCols[1] <- "TRUE"

rTG <- DGEList(counts = dplyr::select(counts, one_of(rTGSamples)),
        group = groups[rTGCols])

## PROCEED WITH JNPL3 ANALYSIS - COME BACK TO rTG

## filter data - require minimum of 100 counts per million for at least 2 samples
d.full <- JNPL3 # keep the old one in case we mess up
dim(d.full)
#39179    24

keep <- rowSums(cpm(JNPL3)>100) >= 2
d <- JNPL3[keep,]
dim(d)
#3133   24

# reset library sizes after filtering
d$samples$lib.size <- colSums(d$counts)

# normalize the data using TMM
d <- calcNormFactors(d)


## Data Exploration
plotMDS(d, method="bcv", col=as.numeric(d$samples$group))
legend("bottomleft", as.character(unique(d$samples$group)), col=1:3, pch=20)

# This plot is worth spending more time with.
# Some observations: JNPL3Plus are generally closer together than the Minus (except X370773 and x374367, which are further away)

## Estimate dispersion.
d1 <- estimateCommonDisp(d, verbose=T) #assume all same for this pass, GLM later
# Disp = 0.09527 , BCV = 0.3087 

d1 <- estimateTagwiseDisp(d1)

plotBCV(d1) #plots the tagwise biological coefficient of variation (square root of dispersions) against log2-CPM.

# observation - looks like a discontinuity in the dispersions around logCPM ~ 7

# use GLM estimate of dispersion
design.mat <- model.matrix(~ 0 + d$samples$group)
colnames(design.mat) <- levels(d$samples$group)
d2 <- estimateGLMCommonDisp(d,design.mat)
d2 <- estimateGLMTrendedDisp(d2,design.mat) #method = "power" also a possibility
# You can change method to "auto", "bin.spline", "power", "spline", "bin.loess".
# The default is "auto" which chooses "bin.spline" when > 200 tags and "power" otherwise.
d2 <- estimateGLMTagwiseDisp(d2, design.mat)
plotBCV(d2) #nice!

# Look at DE with exact test.
de_common <- exactTest(d1) # naive method where we only fit a common dispersion.
de_glm  <- exactTest(d2) # DE with GLM dispersion estimate

# Look at DE with glmFit.

topTags(de_common, n = 10)
topTags(de_glm, n = 10)

de2 <- decideTestsDGE(de_glm, adjust.method="BH", p.value=0.05)
summary(de2)
# The total number of differentially expressed genes at FDR < 0.05 is:
#[,1]
#-1    2 # down-regulated
#0  3131 # not differentially expressed
#1     0 # up-regulated

# box plot of top differentially expressed genes across samples
# I think I want to plot the raw counts (d) for a given gene
# e.g. plot genes corresponding to rownames(topTags(de_glm, n = 10))

topTen <- rownames(topTags(de_glm, n = 10))

toPlot <- d[topTen] #10 rows, 24 columns

boxplot(toPlot$counts[,JNPL3Plus], use.cols = FALSE)
boxplot(toPlot$counts[,JNPL3Minus], use.cols = FALSE)

# munge to use with ggplot
test  <- data.frame(toPlot$samples$group, t(toPlot$counts))

toPlot2 <- melt(test)

ggplot(toPlot2, aes(factor(variable), value)) + geom_boxplot(aes(fill = toPlot.samples.group))

###
#plotSmear generates a plot of the tagwise log-fold-changes against log-cpm (analogous to an MA-plot for microarray #data). DE tags are highlighted on the plot -- but maybe only appropriate for de1 set?:
  
de2tags <- rownames(d2)[as.logical(de2)] 
plotSmear(de_glm, de.tags=de2tags)
abline(h = c(-2, 2), col = "blue")

#likelihood  ratio test should be for de2 set- something is amiss
design.mat
fit <- glmFit(d2, design.mat)
lrt <- glmLRT(fit, contrast = c(1,0))
topTags(lrt, n=10) # all have FDR = 0... seems odd...

de2 <- decideTestsDGE(lrt, adjust.method="BH", p.value = 0.05)
de2tags <- rownames(d2)[as.logical(de2)]
plotSmear(lrt, de.tags=de2tags) # this is not what I expected...
abline(h = c(-2, 2), col = "blue")

## Build summary table
top <- topTags(de_glm, n=50)

ensembl=useMart("ensembl", dataset="mmusculus_gene_ensembl")
geneNames <- getBM("external_gene_name", filters = "ensembl_gene_id", values = rownames(top), ensembl)

top <- cbind(top, geneNames)
colnames(top)[5] <- "Gene_Name"

## add mean gene level, SD, kurtosis, 95% CI of effect size or SD of effect size to summary table

means <- rowMeans(d2$counts[rownames(top), ])
top <- cbind(top, means)
colnames(top)[6] <- "Mean_count"

stdevs <- apply(d2$counts[rownames(top), ], 1, sd)
top <- cbind(top, stdevs)
colnames(top)[7] <- "Count_Std_Dev"

kurtosises  <- apply(d2$counts[rownames(top), ], 1, kurtosis)
top <- cbind(top, kurtosises)
colnames(top)[8] <- "Count_Kurtosis"

# todo: add 95% CI of effect size or SD of effect size to summary table - NO: https://support.bioconductor.org/p/61640/


# key question: how is this done: "Known covariate factors, including sex, race, age, RIN, PMI, batch and site, were corrected using a linear model to remove the confounding effects." ? Is the linear model when dispersion is calculated? 

# google "correcting rnaseq for known covariates" 

# first, ID covariates of interest - 
# https://www.biostars.org/p/97624/
# http://stats.stackexchange.com/questions/59879/logistic-model-what-is-more-important-anova-chi-sq-test-or-significance-of-coe

# then build linear model with that covariate as factor

