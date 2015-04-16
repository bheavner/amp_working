# Goal: differential expression analysis within tissue of Tau mice, focusing on oldest mice - I want to make a list of the most differentially expressed genes between the oldest JNPL3 mice with + genotypes and those with - genotypes. Then I'll apply the same analysis to the oldest rTG4510 samples (some of which are missing from the normalized RNASeq data at the moment).

### move this section to a function - I generally copy and paste it?
library(dplyr)
library(R.utils)
library(synapseClient)
library(edgeR)
library(biomaRt)

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

normalizedCounts <- read.table(localCountFilePath, header = TRUE, stringsAsFactors = FALSE)

# load covariates file to have handy
covariates <- read.table(getFileLocation(covariatesFile), header = TRUE, sep = ",", stringsAsFactors = FALSE)

## finished function section?

# what are the JNPL3 samples I'm interested in?

colnames(covariates)
# 12 month F -:
minusSamples <- dplyr::select(filter(covariates, Sex == "F" & Experiment == "MAPT_P301L" & Age_months == 12 & Genotype == "-"), Mouse_ID)
# 112-10, 112-11, 186848

# 12 month F +:
plusSamples <- dplyr::select(filter(covariates, Sex == "F" & Experiment == "MAPT_P301L" & Age_months == 12 & Genotype == "+"), Mouse_ID)
# 307983, 307987, 309581, 334446, 334455, 334456

minusSamples <- as.character(minusSamples$Mouse_ID)
plusSamples <- as.character(plusSamples$Mouse_ID)

samplesOfInterest <- c(minusSamples, plusSamples)

# colnames(normalizedCounts) are different than covariate Mouse_IDs, so need to transform 112-10 to 112.10, 112-11 to 112.11, and prepend all with X
samplesOfInterest <- sub('-', '.', samplesOfInterest)
samplesOfInterest <- paste("X", samplesOfInterest, sep="")

# starting with example 4.3 (p. 50) of the vignette at http://www.bioconductor.org/packages/release/bioc/vignettes/edgeR/inst/doc/edgeRUsersGuide.pdf
# I need to assign groups in the DGEList object
group <- c(rep("-", 3), rep("+", 6))

# build a DGEList object with these samples
JNPL3Samples <- DGEList(counts=dplyr::select(normalizedCounts, one_of(samplesOfInterest)), group = group)

# Now a detour: "classic edgeR for differential expression" via http://davetang.org/muse/2012/01/19/the-dgelist-object-in-r/
d <- estimateCommonDisp(JNPL3Samples,verbose=T)
# Disp = 0.03638, BCV = 0.1907

d <- estimateTagwiseDisp(d)
de.tgw <- exactTest(d)
summary(decideTestsDGE(de.tgw, p.value=0.01))
# Here the entries for -1, 0 and 1 are for down-regulated, non-differentially expressed and up-regulated tags respectively.
#   [,1]
#-1    35
#0  39117
#1     27

topTags(de.tgw)
#Comparison of groups:  +-- (BH: WHY ARE THERE 3?? B/C the middle "-" is a separator)
#    logFC   logCPM       PValue          FDR
#ENSMUSG00000047676 -7.071383 4.410832 1.283144e-48 5.027230e-44
#ENSMUSG00000027014 -4.346475 7.093599 6.870511e-45 1.345899e-40
#ENSMUSG00000081824 -5.361503 3.252118 2.072516e-25 2.706637e-21
#ENSMUSG00000027559 -4.412237 3.759666 5.494069e-14 5.381303e-10
#ENSMUSG00000046623 -5.339283 2.485383 7.461687e-13 5.846829e-09
#ENSMUSG00000062353 -7.784012 3.783870 1.917639e-12 1.252186e-08
#ENSMUSG00000028452 -1.340927 5.793819 4.375459e-12 2.448944e-08
#ENSMUSG00000022193  1.731153 5.462969 7.028475e-11 3.442108e-07
#ENSMUSG00000091400 -5.962422 2.228828 2.127539e-10 9.261648e-07
#ENSMUSG00000062515 -3.577321 2.805075 2.889561e-10 1.132101e-06

#now say if you wanted to know the tagwise dispersions for each gene
names(d)

normalizedCounts$twd <- d$tagwise.dispersion
head(normalizedCounts)

#plot the tagwise dispersions
#histogram at the end of the post
hist(normalizedCounts$twd, breaks=20, xlim=c(0,3))

#now if you wanted to save the fold changes and p values per gene
names(de.tgw)
#[1] "table"      "comparison" "genes"
head(de.tgw$table)
normalizedCounts <- cbind(normalizedCounts, de.tgw$table)
head(normalizedCounts)

#I want the fdr
normalizedCounts$PValue_fdr <- p.adjust(method="fdr",p=normalizedCounts$PValue)
head(normalizedCounts)

#sanity check with the decideTestsDGE() call
table(normalizedCounts$PValue_fdr<0.01)

#FALSE  TRUE
#39117    62

#to write out this useful data frame
write.table(normalizedCounts, file="JNPL3_12_M.tsv", quote=F)
sessionInfo()

# so

# from https://web.stanford.edu/class/bios221/labs/rnaseq/lab_4_rnaseq.html
# plot the MDS
plotMDS(d, method="bcv", col=as.numeric(d$samples$group))
legend("bottomleft", as.character(unique(d$samples$group)), col=1:3, pch=20)

# I note that the + and - don't separate very nicely - x307983 and x309581 are far away from other + samples.

# plot the tagwise biological coefficient of variation (square root of dispersions) against log2-CPM.
plotBCV(d)
# a single estimate for the coefficient of variation might not be a bad model since tagwise dispersion does not increase much as the counts per million (CPM) increases (at least, once above a certain average log CPM). But we can fit other models to dispersion

# I want to get list of top differentially expressed genes with gene names (not just Ensembl IDs)
top <- topTags(de.tgw, n=50)

ensembl=useMart("ensembl", dataset="mmusculus_gene_ensembl")
geneNames <- getBM("external_gene_name", filters = "ensembl_gene_id", values = rownames(top), ensembl)

top <- cbind(top, geneNames)


# back to example 4.3 (p. 50) of the vignette at http://www.bioconductor.org/packages/release/bioc/vignettes/edgeR/inst/doc/edgeRUsersGuide.pdf


## NEXT, DO WITH DESeq (and DESeq2?) (see https://web.stanford.edu/class/bios221/labs/rnaseq/lab_4_rnaseq.html)
