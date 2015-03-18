# Goal: DE analysis of JNPL3 Tau mice (transgenic vs wild type).

# subgoals: 
# impact of more or less stringent criteria (particularly removing low readcount and dispersion models)?
# how to best analyze and account for covariates?

# future goal: 
# more exploratory data analysis - visualization of readcounts before and after various processing steps.

# output should be a table of top most differentially expressed genes, with the following columns: gene ID; logFc' logCPM; pValue' FDR' external_gene_name; mean gene level; sd; kurtosis; (I decided not to include 95% CI of effect size or SD of effect size - see https://support.bioconductor.org/p/61640/); histogram of values

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

## define groups: JNPL3+ (transgenic), JNPL3- (WT), (ignore rTG+ rTG- for now)
# (there MUST be a cleaner way to do this)

# JNPL3+:
JNPL3Plus <- dplyr::select(filter(covariates, Experiment == "MAPT_P301L" & Genotype == "+"), Mouse_ID)
# 15 samples

# JNPL3-:
JNPL3Minus <- dplyr::select(filter(covariates, Experiment == "MAPT_P301L" & Genotype == "-"), Mouse_ID)
# 9 samples

JNPL3Plus <- as.character(JNPL3Plus$Mouse_ID)
JNPL3Minus <- as.character(JNPL3Minus$Mouse_ID)

JNPL3PlusCols <- is.element(as.character(colnames(counts)), JNPL3Plus)
JNPL3MinusCols <- is.element(as.character(colnames(counts)), JNPL3Minus)

groups <- c(rep(0, length(counts[1,])))
groups[JNPL3PlusCols] <- "JNPL3Plus"
groups[JNPL3MinusCols] <- "JNPL3Minus"

JNPL3Samples <- c(JNPL3Plus, JNPL3Minus)

JNPL3Cols <- is.element(as.character(colnames(counts)), JNPL3Samples)


## Make DGEList object to start working on DE analysis
JNPL3 <- DGEList(counts = dplyr::select(counts, one_of(JNPL3Samples)),
                 group = groups[JNPL3Cols])

## filter data 
# require minimum of 100 counts per million for at least 2 samples
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
d <- calcNormFactors(d, method = "TMM")

plotMDS(d, method="bcv", col=as.numeric(d$samples$group))
legend("bottomleft", as.character(unique(d$samples$group)), col=1:3, pch=20)

# This plot is worth spending more time with.
# Some observations: JNPL3Plus are generally closer together than the Minus (except X370773 and x374367, which are further away)

## FIRST, CLASSIC EDGER ANALYSIS
# Estimate dispersion - first simple way.
d1 <- estimateCommonDisp(d, verbose=T) #assume all same for this pass, GLM later
# Disp = 0.09527 , BCV = 0.3087 

d1 <- estimateTagwiseDisp(d1)

plotBCV(d1) #plots the tagwise biological coefficient of variation (square root of dispersions) against log2-CPM.
# observation - looks like a discontinuity in the dispersions around logCPM ~ 7

# Look at DE with exact test
de.tgw <- exactTest(d1)
summary(decideTestsDGE(de.tgw, p.value=0.01))
# Here the entries for -1, 0 and 1 are for down-regulated, non-differentially expressed and up-regulated tags respectively.
#   [,1]
#-1    1
#0  3132
#1     0

topTags(de.tgw)
#Comparison of groups:  JNPL3Plus-JNPL3Minus 
#logFC    logCPM       PValue         FDR
#ENSMUSG00000061808 -5.7919888  7.345933 8.194772e-07 0.002567422
#ENSMUSG00000079037  0.7948772 10.942645 1.897399e-05 0.026610982
#ENSMUSG00000042109 -0.9310217  7.857671 2.548131e-05 0.026610982
#ENSMUSG00000013275  0.4216274  7.707057 8.475993e-05 0.066388217
#ENSMUSG00000004187 -0.2979301  7.847114 3.764522e-04 0.235884923
#ENSMUSG00000025780  0.6013725  7.271255 4.554743e-04 0.237833476
#ENSMUSG00000020932  0.9282346 12.041427 1.338830e-03 0.599222195

# box plot of top differentially expressed genes across samples
topTen <- rownames(topTags(de.tgw, n = 10))
toPlot <- d[topTen] #10 rows, 24 columns
top <- topTags(de.tgw, n=50)

# get gene names for x axis labels
ensembl=useMart("ensembl", dataset="mmusculus_gene_ensembl")
geneNames <- getBM("external_gene_name", 
                   filters = "ensembl_gene_id", 
                   values = rownames(top)[1:10], 
                   ensembl)

# munge to use with ggplot and plot
test  <- data.frame(toPlot$samples$group, t(toPlot$counts))
toPlot2 <- melt(test)
ggplot(toPlot2, aes(factor(variable), value)) + 
  geom_boxplot(aes(fill = toPlot.samples.group)) +
  scale_x_discrete(labels = geneNames[1:10, ])



## NEXT, glm EDGER ANALYSIS WITH EXACT TEST
# To compare, use GLM estimate of dispersion
design.mat <- model.matrix(~0 + d$samples$group)
colnames(design.mat) <- levels(d$samples$group)
d2 <- estimateGLMCommonDisp(d,design.mat)
d2 <- estimateGLMTrendedDisp(d2,design.mat) #method = "power" also a possibility
# You can change method to "auto", "bin.spline", "power", "spline", "bin.loess".
# The default is "auto" which chooses "bin.spline" when > 200 tags and "power" otherwise.
d2 <- estimateGLMTagwiseDisp(d2, design.mat)
plotBCV(d2) # seems nice - blue line looks like it follows points...

# Look at DE with exact test (same test as with simple model)
de_common <- exactTest(d2) # DE with GLM dispersion estimate with exact Test.
summary(decideTestsDGE(de_common, p.value=0.01))
#[,1]
#-1    1
#0  3132
#1     0

topTags(de_common)
#Comparison of groups:  JNPL3Plus-JNPL3Minus 
#logFC    logCPM       PValue          FDR
#ENSMUSG00000061808 -5.7919834  7.346404 1.967464e-08 6.164065e-05
#ENSMUSG00000042109 -0.9310217  7.857690 2.548227e-05 3.991797e-02
#ENSMUSG00000013275  0.4216278  7.707029 8.200863e-05 8.564435e-02
#ENSMUSG00000079037  0.7948759 10.942642 2.462941e-04 1.929099e-01
#ENSMUSG00000004187 -0.2979300  7.847129 3.920177e-04 2.344835e-01
#ENSMUSG00000025780  0.6013726  7.271208 4.490588e-04 2.344835e-01
#ENSMUSG00000061740  0.5970873  6.940286 1.551455e-03 6.943871e-01

# box plot of top differentially expressed genes across samples
topTen <- rownames(topTags(de_common, n = 10))
toPlot <- d[topTen] #10 rows, 24 columns
top <- topTags(de_common, n=50)

# get gene names for x axis labels
ensembl=useMart("ensembl", dataset="mmusculus_gene_ensembl")
geneNames <- getBM("external_gene_name", 
                   filters = "ensembl_gene_id", 
                   values = rownames(top)[1:10], 
                   ensembl)

# munge to use with ggplot and plot
test  <- data.frame(toPlot$samples$group, t(toPlot$counts))
toPlot2 <- melt(test)
ggplot(toPlot2, aes(factor(variable), value)) + 
  geom_boxplot(aes(fill = toPlot.samples.group)) +
  scale_x_discrete(labels = geneNames[1:10, ])


# THIRD, glm EDGER ANALYSIS WITH glmLRT TEST
fit <- glmFit(d2, design.mat) # Andrew used fit <- glmFit(d2, design.mat, dispersion=dge$tagwise.dispersion)
de_glm  <- glmLRT(fit, coef = 2) # DE with GLM dispersion estimate with GLMRT test

#de_glm  <- glmLRT(fit, contrast=c(1,-1,0)) # Andrew's version - doesn't work for me

summary(decideTestsDGE(de_glm, p.value=0.01))
#[,1]
#-1 3133
#0     0
#1     0

topTags(de_glm)
# Look at DE with glmFit.
#Coefficient:  JNPL3Minus 
#logFC   logCPM       LR PValue FDR
#ENSMUSG00000024998 -13.48575 6.758922 2214.383      0   0
#ENSMUSG00000059602 -13.48327 6.633828 1781.801      0   0
#ENSMUSG00000068735 -13.41580 6.638248 1770.110      0   0
#ENSMUSG00000024816 -13.40593 6.803815 2041.364      0   0
#ENSMUSG00000061740 -13.39296 6.940286 1634.186      0   0
#ENSMUSG00000029765 -13.38953 6.653501 1844.233      0   0
#ENSMUSG00000039474 -13.35474 6.795300 2883.819      0   0

# why all 0 PValue and 0 FDR? Seems like I may need a different test?

# box plot of top differentially expressed genes across samples
topTen <- rownames(topTags(de_glm, n = 10))
toPlot <- d[topTen] #10 rows, 24 columns
top <- topTags(de_glm, n=50)

# get gene names for x axis labels
ensembl=useMart("ensembl", dataset="mmusculus_gene_ensembl")
geneNames <- getBM("external_gene_name", 
                   filters = "ensembl_gene_id", 
                   values = rownames(top)[1:10], 
                   ensembl)

# munge to use with ggplot and plot
test  <- data.frame(toPlot$samples$group, t(toPlot$counts))
toPlot2 <- melt(test)
ggplot(toPlot2, aes(factor(variable), value)) + 
  geom_boxplot(aes(fill = toPlot.samples.group)) +
  scale_x_discrete(labels = geneNames[1:10, ])


## FINALLY,
## Build summary table with desired outputs
top <- topTags(de_glm, n=50)

ensembl=useMart("ensembl", dataset="mmusculus_gene_ensembl")
geneNames <- getBM("external_gene_name", filters = "ensembl_gene_id", values = rownames(top), ensembl)

top <- cbind(top, geneNames)
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