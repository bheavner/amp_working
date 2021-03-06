# Goal: DE analysis of JNPL3 Tau mice (transgenic (+) vs wild type (-) ).

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
#3133   24 #seth says this is too stringent - 70% of genes in genome are expressed in brain - ~20k features normal

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

# get gene names for x axis labels -- ORDER ISN'T PRESERVED FROM BIOMART QURY!
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


# THIRD, glm EDGER ANALYSIS WITH glmLRT TEST
fit <- glmFit(d2, design.mat) 
# Andrew used fit <- glmFit(d2, design.mat, dispersion=dge$tagwise.dispersion) but since "If NULL will be extracted from y, with order of precedence: tagwise dispersion, trended dispersions, common dispersion.", I'll leave it out and let it get it from d2.
de_glm  <- glmLRT(fit, contrast=c(-1,1)) # -1*JNPL3Minus 1*JNPL3Plus 
#so a positive logFC means transgenic (JNPL3Plus) is higher expression than WT (JNPL3Minus).

summary(decideTestsDGE(de_glm, p.value=0.01))
#[,1]
#-1    1
#0  3132
#1     0

topTags(de_glm)
# Look at DE with glmFit.
#Coefficient:  -1*JNPL3Minus 1*JNPL3Plus 
#logFC    logCPM        LR       PValue          FDR
#ENSMUSG00000061808 -5.7919834  7.346404 32.212142 1.382243e-08 4.330566e-05
#ENSMUSG00000042109 -0.9310217  7.857690 17.670991 2.626011e-05 4.113646e-02
#ENSMUSG00000013275  0.4216278  7.707029 15.571870 7.942754e-05 8.294882e-02
#ENSMUSG00000079037  0.7948759 10.942642 13.563288 2.306527e-04 1.806587e-01
#ENSMUSG00000004187 -0.2979300  7.847129 12.543855 3.975113e-04 2.236823e-01
#ENSMUSG00000025780  0.6013726  7.271208 12.404182 4.283733e-04 2.236823e-01
#ENSMUSG00000061740  0.5970873  6.940286 10.106146 1.477756e-03 6.614016e-01
#ENSMUSG00000024164  0.9961478  6.799560  9.496849 2.058251e-03 6.975576e-01
#ENSMUSG00000053702  0.2586839  7.857723  9.400092 2.169745e-03 6.975576e-01
#ENSMUSG00000027479  0.3252234  7.901499  9.087314 2.573886e-03 6.975576e-01

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


# FOURTH, glm EDGER ANALYSIS WITH glmQLFTest TEST - see examples in ?glmQLFit for futher refinement

QLFTfit <- glmQLFit(d2, design.mat) 
de_glmQLFT  <- glmQLFTest(QLFTfit, contrast=c(-1,1)) # -1*JNPL3Minus 1*JNPL3Plus 
#so a positive logFC means transgenic (JNPL3Plus) is higher expression than WT (JNPL3Minus).

summary(decideTestsDGE(de_glmQLFT, p.value=0.01))
#[,1]
#-1    0
#0  3133
#1     0

topTags(de_glmQLFT)
#Coefficient:  -1*JNPL3Minus 1*JNPL3Plus 
#logFC    logCPM         F       PValue       FDR
#ENSMUSG00000061808 -5.7908017  7.346404 17.074116 0.0003596785 0.6418135
#ENSMUSG00000013275  0.4216076  7.707029 17.831728 0.0004097117 0.6418135
#ENSMUSG00000042109 -0.9309660  7.857690 14.350161 0.0008651960 0.9035530
#ENSMUSG00000025780  0.6014224  7.271208 11.305397 0.0025171312 0.9998911
#ENSMUSG00000079037  0.7948747 10.942642 20.396186 0.0026141586 0.9998911
#ENSMUSG00000061740  0.5971518  6.940286  9.374968 0.0052459760 0.9998911
#ENSMUSG00000020932  0.9282341 12.041424  9.028230 0.0060183845 0.9998911
#ENSMUSG00000027479  0.3252075  7.901499 10.264608 0.0065121841 0.9998911
#ENSMUSG00000031700  0.3239644  7.258742  9.930533 0.0071769211 0.9998911
#ENSMUSG00000031596  0.4337531  7.173476  8.483408 0.0074956168 0.9998911

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

write.table(top, file="JNPL3_tg_vs_wt.tsv", quote=F)