# Starting to explore normalized rnaseq data from tau mice. Goal: diferential expression analysis.

library(synapseClient)
library(R.utils)
library(dplyr)
library(edgeR)

#Login to Synapse using credentials saved in .synapseConfig file
synapseLogin()

# get the normalized readcount and covariates files from synapse
countFile <- synGet('syn3219030') # start with the normalized gene counts - can play with transcript counts later
covariatesFile <- synGet('syn3219041')

# unzip count file and load for processing
localCountFilePath <- getFileLocation(countFile)

if(!file.exists(sub('.gz', '', localCountFilePath))) {
    gunzip(localCountFilePath)
}

localCountFilePath <- sub('.gz', '', localCountFilePath) #trim the .gz suffix

normalizedCounts <- read.table(localCountFilePath, header = TRUE)

# load covariates file to have handy
covariates <- read.table(getFileLocation(covariatesFile), header = TRUE, sep = ",")

# As my first comparison, I'll do rTg4510 (F) vs JNPL3 using DESeq2. First, use dplyr to find sample numbers.
# plan - use covariates file to select counts from countfile.
dim(covariates) #(60 x 11)

row.names(covariates)
#"1"  "2"  "3"  "4"  "5"  "6"  "7"  "8"  "9"  "10" "11" "12" "13" "14" "15" "16" "17" "18" "19" "20" "21" "22" "23" "24" "25" "26" "27" "28" "29" "30" "31" "32" "33" "34" "35" "36" "37" "38" "39" "40" "41" "42" "43" "44" "45" "46" "47" "48" "49" "50" "51" "52" "53" "54" "55" "56" "57" "58" "59" "60"

colnames(covariates)
# "Mouse_ID" "Experiment" "RIN" "Genotype" "Sex" "Age_months" "RLIMS.ID" "Seq.Run.ID" "Lane.Number" "Clusters" "Raw_RNAseq_file_name"

# this works:
filter(covariates, Sex == "F")

# how many JNPL3 females in covariates file?
dim(filter(covariates, Sex == "F" & Experiment == "MAPT_P301L"))
# 24 x 11

# how many rTg4510 females in covariates file?
dim(filter(covariates, Sex == "F" & Experiment == "rTG4510"))
# 18 x 11

# what are the sample IDs for the rTg4510 females?
rTg4510FemaleIds <- filter(covariates, Sex == "F" & Experiment == "rTG4510")$Mouse_ID
JPNL3FemaleIds <- filter(covariates, Sex == "F" & Experiment == "MAPT_P301L")$Mouse_ID

# prepend "X" to JNPL3FemaleIds to match ID in readcount colnames


# want to select readcounts for those Mouse_ID values

# from http://www.bioconductor.org/packages/release/bioc/vignettes/edgeR/inst/doc/edgeRUsersGuide.pdf, edgeR analysis could be as easy as:
group <- c(1,1,2,2)
D <- DGEList(counts=y, group=group)
D <- estimateCommonDisp(D)
D <- estimateTagwiseDisp(D)

# p. 18 - how to define design matrix? -- see p25, p27, 29
D <- estimateGLMCommonDisp(D, design)
D <- estimateGLMTrendedDisp(D, design)
D <- estimateGLMTagwiseDisp(D, design)

#As a brief example, consider a situation in which are three treatment groups, each with
#two replicates, and the researcher wants to make pairwise comparisons between them. A
#linear model representing the study design can be tted to the data with commands such as:
group <- factor(c(1,1,2,2,3,3))
design <- model.matrix(~group)
fit <- glmFit(y,design)

#To compare 2 vs 1:
lrt.2vs1 <- glmLRT(fit,coef=2)
topTags(lrt.2vs1)

#To compare 3 vs 1:
lrt.3vs1 <- glmLRT(fit,coef=3)
#To compare 3 vs 2:
lrt.3vs2 <- glmLRT(fit,contrast=c(0,-1,1))
#The contrast argument in this case requests a statistical test of the null hypothesis that
#coeffcient3 - coeffcient2 is equal to zero.
#To find genes dierent between any of the three groups:
lrt <- glmLRT(fit,coef=2:3)
topTags(lrt)


# With that reading, I changed my mind on design. I'll start with looking at te rTg4510 - vs + genotype first, following example 4.3 (p. 50) of the vignette at http://www.bioconductor.org/packages/release/bioc/vignettes/edgeR/inst/doc/edgeRUsersGuide.pdf

# list of rTg4510 samples:
select(normalizedCounts, c(LP_8_3, LP_9_6, R.f3_8, R695_8, R695_9, R698_1, R708_1, R716_4, R719_6, R720_5, R749_4, R753_2, R760_4, R763_5, R773_1, R791.3, R803_3, R804_3, R804_4, R816_4, R817_6, R819_3, R821_1))

rtgCounts <- DGEList(counts=select(normalizedCounts, c(LP_8_3, LP_9_6, R.f3_8, R695_8, R695_9, R698_1, R708_1, R716_4, R719_6, R720_5, R749_4, R753_2, R760_4, R763_5, R773_1, R791.3, R803_3, R804_3, R804_4, R816_4, R817_6, R819_3, R821_1)))

# add design matrix
genotype <- c("-","-","+","+","+","+","+","-","-","-","-","-","-","-","+","-","+","-","+","+","+","+","-")

plotMDS(rtgCounts)


# what if do MDS plot of both rTg4510 and JNPL3 mice?
allCounts <- DGEList(counts=normalizedCounts) #17.3 Mb - big?
plotMDS(allCounts)

# make MDS plots prettier - dots instead of labels, colored by rTg or JNPL3
# see ?plotMDS.DGEList on methods

plotMDS(rtgCounts, method="bcv", col=as.numeric(rtgCounts$samples$name))
legend("bottomleft", as.character(unique(dge$samples$name)), col=1:4, pch=20)

plotMDS(allCounts, method="bcv", col=as.numeric(rtgCounts$samples$name))
legend("bottomleft", as.character(unique(dge$samples$name)), col=1:4, pch=20)
# Should be able to run Cory's EdgeR script here... Also generates prettier MDS plot
