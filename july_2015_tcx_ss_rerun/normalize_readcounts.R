# This script normalizes the merged count_files produced by running merge_file_counts.R
# on SNAPR output. It does so by dowloading the merged read count files, processing them,
# and uploading new synapse objects with appropriate provenance.

library(synapseClient)
library(R.utils)
library(edgeR)

#Login to Synapse using credentials saved in .synapseConfig file
synapseLogin()

codeFile <- ("https://github.com/PriceLab/AMP-ad/tree/0.2/MBB/TCX/normalize_readcounts.R")

# The files to normalize are:

# AMP-AD_MayoBB_UFL-Mayo-ISB_IlluminaHiSeq2000_TCX_GeneCounts_JulyRerun.txt.gz ('syn4650257')
# AMP-AD_MayoBB_UFL-Mayo-ISB_IlluminaHiSeq2000_TCX_TranscriptCounts_JulyRerun.txt.gz ('syn4650258')

countFileSynapseIDs <- c('syn4650257', 'syn4650258')

for (mergedCountFile in countFileSynapseIDs) {
    message("Normalizing ", mergedCountFile)

    # Download file from Synapse
    originalCountFile <- synGet(mergedCountFile)

    # unzip file and load for processing
    localFilePath <- getFileLocation(originalCountFile)

    if(!file.exists(substr(localFilePath, 1, nchar(localFilePath) - 3))) {
        gunzip(localFilePath)
    }

    localFilePath <- sub('.gz', '', localFilePath) #trim the .gz suffix

    transposedCounts <- read.table(localFilePath, header = TRUE, check.names = FALSE)

    # make DGEList object
    expr <- DGEList(transposedCounts, group = rep(1, ncol(transposedCounts)))

    # calculate normalization factors
    normFactors <- calcNormFactors(expr, method = ("TMM"))
    
    # If you get "error: Error in quantile.default(x, p = p) : 
    # missing values and NaN's not allowed if 'na.rm' is FALSE"

    # expr$samples$lib.size shows a library size of 0 for problematic samples

    # use normaliztion factors to calculate cpm -
    # per https://www.biostars.org/p/84087/, that's calculated as
    # count / (library size * normalization factor))

    normalizedCpm <- cpm(normFactors)

    # write the data to local dir
    newFileName <- sub('_JulyRerun.txt.gz', '', originalCountFile$properties$name)
    newFileName <- paste0(newFileName, "_Normalized_JulyRerun.txt", sep="")

    write.table(format(normalizedCpm, scientific = FALSE, digits = 5), 
                newFileName, quote = FALSE, sep = "\t", row.names = TRUE)
    
    # package it up, then create a Synapse object for the output file and upload with provenance

    gzip(newFileName)

    newFileName <- paste0(newFileName, ".gz", sep="")

    parentId <- originalCountFile$properties$parentId

    normalizedCountFile <- File(newFileName, parentId = parentId)

    normalizedCountFile <- synStore(normalizedCountFile,
                                    activityName="CPM (using TMM) from edgeR normalization",
                                    used=list(list(name = "normalize_readcounts.R",
                                                   url = codeFile, wasExecuted = T),
                                              list(entity=originalCountFile,
                                                   wasExecuted=F)))
}

