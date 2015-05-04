# This script normalizes the merged count_files produced by running merge_file_counts.R
# on SNAPR output. It does so by dowloading the merged read count files, processing them,
# and uploading new synapse objects with appropriate provenance.

library(synapseClient)
library(R.utils)
library(edgeR)

#Login to Synapse using credentials saved in .synapseConfig file
synapseLogin()

codeFile <- ("https://github.com/TODO")

# The files for this batch are:

# AMP-AD_TAUAPPms_UFL-Mayo-ISB_IlluminaHiSeq2000_Tau_GeneCounts_Updated.txt.gz - syn3910122
# AMP-AD_TAUAPPms_UFL-Mayo-ISB_IlluminaHiSeq2000_Tau_TranscriptCounts_Updated.txt.gz - syn3910523


countFileSynapseIDs <- c('syn3910122', 'syn3910523')

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

    Counts <- read.table(localFilePath, header = TRUE, check.names = FALSE)

    # make DGEList object
    expr <- DGEList(Counts, group = rep(1, ncol(Counts)))

    # calculate normalization factors
    normFactors <- calcNormFactors(expr, method = ("TMM"))

    # use normaliztion factors to calculate cpm -
    # per https://www.biostars.org/p/84087/, that's calculated as
    # count / (library size * normalization factor))

    normalizedCpm <- cpm(normFactors)

    # write the data to local dir

    newFileName <- sub('_transposed.txt.gz', '', originalCountFile$properties$name) #legacy?
    newFileName <- paste0(newFileName, "_normalized.txt", sep="")

    write.table(format(normalizedCpm, scientific = FALSE, digits = 5), newFileName, quote = FALSE, sep = "\t", row.names = TRUE, col.names = NA)    
    
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
