# Goal: first-pass DE analysis of CRND8 APP mice (transgenic (+) vs wild type (-) ).

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

# get the readcount file and covariates file from synapse
countFile <- synGet('syn3483519') 
covariatesFile <- synGet('syn3483880') # the working dir copy

# unzip count file and load for processing
localCountFilePath <- getFileLocation(countFile)

if(!file.exists(sub('.gz', '', localCountFilePath))) {
  gunzip(localCountFilePath)
}

localCountFilePath <- sub('.gz', '', localCountFilePath) #trim the .gz suffix

counts <- read.table(localCountFilePath, header = TRUE, stringsAsFactors = FALSE)

# load covariates file to have handy
covariates <- read.table(getFileLocation(covariatesFile), header = TRUE, stringsAsFactors = FALSE)

# Make DGEList object with CRND8 counts and groups defined by covariates$Genotype

# easiest just to subset out the CRND8 samples at the beginning
CRND8Covariates  <- covariates[covariates$Experiment == "CRND8",]
CRND8counts  <- counts[, colnames(counts) %in% CRND8Covariates$Mouse_ID]
orderedGenotype  <- covariates$Genotype[order(match(CRND8Covariates$Mouse_ID, colnames(CRND8counts)))]

CRND8  <- DGEList(counts = CRND8counts, group = orderedGenotype)

## filter data 
# require minimum of 100 counts per million for at least 2 samples
d.full <- CRND8 # keep the old one in case we mess up
dim(d.full)
#39179    88

keep <- rowSums(cpm(CRND8)>100) >= 2
d <- CRND8[keep,]
dim(d)
#2796   88 #seth says this is too stringent - 70% of genes in genome are expressed in brain - ~20k features normal. In this case, 70% is 27425 features. So, try lower minimum

keep <- rowSums(cpm(CRND8)>10) >= 2
d <- CRND8[keep,]
dim(d)
# 11252 88

keep <- rowSums(cpm(CRND8)>1) >= 2
d <- CRND8[keep,]
dim(d)
# 16175 88

# I'll stick with the strict one
keep <- rowSums(cpm(CRND8)>100) >= 2
d <- CRND8[keep,]
dim(d)

# or use raw readcount of >10 for at least 2 samples:
keep <- rowSums(CRND8$counts>10) >= 2
d <- CRND8[keep,]
dim(d)
# 19394    88


# reset library sizes after filtering
d$samples$lib.size <- colSums(d$counts)

# normalize the data using TMM
d <- calcNormFactors(d, method = "TMM")

plotMDS(d, method="bcv", col=as.numeric(d$samples$group))
legend("bottomleft", as.character(unique(d$samples$group)), col=1:3, pch=20)

# This plot is worth spending more time with.
# Some observations: looks like there's a key covariate separating samples into 2 groups more than genotype (I suspect male/female)
# So, let's segregate data and just do females for DE analysis. Can do differently later.


CRND8FemaleCovariates  <- CRND8Covariates[CRND8Covariates$Sex == "F",]
CRND8FemaleCounts  <- CRND8counts[, colnames(CRND8counts) %in% CRND8FemaleCovariates$Mouse_ID]
orderedFemaleGenotype  <- covariates$Genotype[order(match(CRND8FemaleCovariates$Mouse_ID, colnames(CRND8FemaleCounts)))]

CRND8Female  <- DGEList(counts = CRND8FemaleCounts, group = orderedFemaleGenotype)

## filter data 
# require minimum of 100 counts per million for at least 2 samples
d.full <- CRND8Female # keep the old one in case we mess up
dim(d.full)
#39179    46

keep <- rowSums(cpm(CRND8Female)>100) >= 2
d <- CRND8Female[keep,]
dim(d)
#2729   46

# reset library sizes after filtering
d$samples$lib.size <- colSums(d$counts)

# normalize the data using TMM
d <- calcNormFactors(d, method = "TMM")

plotMDS(d, method="bcv", col=as.numeric(d$samples$group))
legend("bottomleft", as.character(unique(d$samples$group)), col=1:3, pch=20)

# This plot is worth spending more time with.
# Some observations: Removing gender did a lot, but it looks like there's another covariate separating samples into 2 groups. Proceed anyway.


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
#-1    0
#0  2729
#1     0

topTags(de.tgw)
# Comparison of groups:  +-- 
#   logFC   logCPM      PValue       FDR
# ENSMUSG00000034730 -0.2109890 9.213325 0.003371372 0.9474738
# ENSMUSG00000026185  0.2494932 8.110462 0.003544733 0.9474738
# ENSMUSG00000040479 -0.1793679 9.135856 0.005921796 0.9474738
# ENSMUSG00000044060  0.1529706 7.832681 0.007384821 0.9474738
# ENSMUSG00000003352 -0.1527321 8.094938 0.008676303 0.9474738
# ENSMUSG00000019254 -0.1526298 7.243193 0.010696263 0.9474738
# ENSMUSG00000015094 -0.1493345 7.671115 0.010860216 0.9474738
# ENSMUSG00000029361  0.1737343 7.039045 0.012108308 0.9474738

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
#Comparison of groups:  CRND8Plus-CRND8Minus 
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
de_glm  <- glmLRT(fit, contrast=c(-1,1)) # -1*CRND8Minus 1*CRND8Plus 
#so a positive logFC means transgenic (CRND8Plus) is higher expression than WT (CRND8Minus).

summary(decideTestsDGE(de_glm, p.value=0.01))
#[,1]
#-1    1
#0  3132
#1     0

topTags(de_glm)
# Look at DE with glmFit.
#Coefficient:  -1*CRND8Minus 1*CRND8Plus 
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
de_glmQLFT  <- glmQLFTest(QLFTfit, contrast=c(-1,1)) # -1*CRND8Minus 1*CRND8Plus 
#so a positive logFC means transgenic (CRND8Plus) is higher expression than WT (CRND8Minus).

summary(decideTestsDGE(de_glmQLFT, p.value=0.01))
#[,1]
#-1    0
#0  3133
#1     0

topTags(de_glmQLFT)
#Coefficient:  -1*CRND8Minus 1*CRND8Plus 
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

write.table(top, file="CRND8_tg_vs_wt.tsv", quote=F)