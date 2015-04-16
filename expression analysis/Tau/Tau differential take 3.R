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
#countFile <- synGet('syn3219030') # start with the normalized gene counts - can play with transcript counts later
countFile <- synGet('syn3207189') # the working dir copy
#covariatesFile <- synGet('syn3219041')
covariatesFile <- synGet('syn2875343') # the working dir copy

# unzip count file and load for processing
localCountFilePath <- getFileLocation(countFile)

if(!file.exists(sub('.gz', '', localCountFilePath))) {
    gunzip(localCountFilePath)
}

localCountFilePath <- sub('.gz', '', localCountFilePath) #trim the .gz suffix

normalizedCounts <- read.table(localCountFilePath, header = TRUE, stringsAsFactors = FALSE)

# load covariates file to have handy
covariates <- read.table(getFileLocation(covariatesFile), header = TRUE, stringsAsFactors = FALSE)

## finished function section?

# what are the JNPL3 samples I'm interested in?

colnames(covariates)
# 12 month F -:
minusSamples <- dplyr::select(filter(covariates, Sex == "F" & Experiment == "MAPT_P301L" & Age_months == 12 & Genotype == "-"), Mouse_ID)
# X112.10, X112.11, X186848

# 12 month F +:
plusSamples <- dplyr::select(filter(covariates, Sex == "F" & Experiment == "MAPT_P301L" & Age_months == 12 & Genotype == "+"), Mouse_ID)
# X307983, X307987, X309581, X334446, X334455, X334456

minusSamples <- as.character(minusSamples$Mouse_ID)
plusSamples <- as.character(plusSamples$Mouse_ID)

samplesOfInterest <- c(minusSamples, plusSamples)

# colnames(normalizedCounts) are different than covariate Mouse_IDs, so need to transform 112-10 to 112.10, 112-11 to 112.11, and prepend all with X
#samplesOfInterest <- sub('-', '.', samplesOfInterest)
#samplesOfInterest <- paste("X", samplesOfInterest, sep="")

# starting with example 4.3 (p. 50) of the vignette at http://www.bioconductor.org/packages/release/bioc/vignettes/edgeR/inst/doc/edgeRUsersGuide.pdf
# I need to assign groups in the DGEList object
group <- c(rep("-", 3), rep("+", 6))

# build a DGEList object with these samples
JNPL3Samples <- DGEList(counts=dplyr::select(normalizedCounts, one_of(samplesOfInterest)), group = group)

# Now a detour: "classic edgeR for differential expression" via http://davetang.org/muse/2012/01/19/the-dgelist-object-in-r/
d <- estimateCommonDisp(JNPL3Samples,verbose=T)
# Disp = 0.03631, BCV = 0.1906

d <- estimateTagwiseDisp(d)
de.tgw <- exactTest(d)
summary(decideTestsDGE(de.tgw, p.value=0.01))
# Here the entries for -1, 0 and 1 are for down-regulated, non-differentially expressed and up-regulated tags respectively.
#   [,1]
#-1    35
#0  39118
#1     26

topTags(de.tgw)
#Comparison of groups:  +-- (BH: WHY ARE THERE 3?? B/C the middle "-" is a separator)
#    logFC   logCPM       PValue          FDR
#ENSMUSG00000047676 -7.063692 4.411022 4.741478e-48 1.857664e-43
#ENSMUSG00000027014 -4.345718 7.093495 1.482383e-44 2.903915e-40
#ENSMUSG00000081824 -5.355373 3.255806 4.074659e-25 5.321369e-21
#ENSMUSG00000027559 -4.410820 3.759977 8.791768e-14 8.611317e-10
#ENSMUSG00000046623 -5.330143 2.493234 7.201330e-13 5.642818e-09
#ENSMUSG00000062353 -7.763553 3.787814 2.107697e-12 1.376291e-08
#ENSMUSG00000028452 -1.341611 5.794217 5.237312e-12 2.931323e-08
#ENSMUSG00000022193  1.730715 5.464443 7.877320e-11 3.857819e-07
#ENSMUSG00000091400 -5.936785 2.237956 2.008242e-10 8.742322e-07
#ENSMUSG00000048087 -5.858780 2.193862 4.172765e-10 1.634847e-06

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
#39118    61

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
# a single estimate for the coefficient of variation might not be a bad model since tagwise dispersion does not increase much as the counts per million (CPM) increases (at least, once above a certain average log CPM). But we can fit other models to dispersion later.

# I want to get list of top differentially expressed genes with gene names (not just Ensembl IDs)
top <- topTags(de.tgw, n=50)

ensembl=useMart("ensembl", dataset="mmusculus_gene_ensembl")
geneNames <- getBM("external_gene_name", filters = "ensembl_gene_id", values = rownames(top), ensembl)

top <- cbind(top, geneNames)


# back to example 4.3 (p. 50) of the vignette at http://www.bioconductor.org/packages/release/bioc/vignettes/edgeR/inst/doc/edgeRUsersGuide.pdf


## NEXT, DO WITH DESeq (and DESeq2?) (see https://web.stanford.edu/class/bios221/labs/rnaseq/lab_4_rnaseq.html)
