# This script normalizes the merged count_files produced by running merge_file_counts.R
# on SNAPR output. It does so by dowloading the merged read count files, processing them,
# and uploading new synapse objects with appropriate provenance.

library(synapseClient)
library(R.utils)
library(edgeR)

#Login to Synapse using credentials saved in .synapseConfig file
synapseLogin()

codeFile <- ("https://github.com/bheavner/ampSynapseProjects/blob/5c2fea4f6c4930e623df97eef96e6b0224aaeec1/rnaseqAnalysis/normalize_readcounts.R")

# The files for the first batch are:

# ad_pilot_rnaseq_gene_id_counts_transposed.txt.gz ('syn3191070')
# ad_pilot_rnaseq_transcript_id_counts_transposed.txt.gz ('syn3191083')

# psp_pilot_rnaseq_gene_id_counts_transposed.txt.gz ('syn3191085')
# psp_pilot_rnaseq_transcript_id_counts_transposed.txt.gz ('syn3191122')

# mouse_tau_rnaseq_gene_id_counts_transposed.txt.gz ('syn3192634')
# mouse_tau_rnaseq_transcript_id_counts_transposed.txt.gz ('syn3192651')

countFileSynapseIDs <- c('syn3192634', 'syn3192651')#'syn3191070', 'syn3191083', 'syn3191085', 'syn3191122', 'syn3192634', 'syn3192651')

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

    Counts <- read.table(localFilePath, header = TRUE)

    # make DGEList object
    expr <- DGEList(Counts, group = rep(1, ncol(transposedCounts)))

    # calculate normalization factors
    normFactors <- calcNormFactors(expr, method = ("TMM"))

    # use normaliztion factors to calculate cpm -
    # per https://www.biostars.org/p/84087/, that's calculated as
    # count / (library size * normalization factor))

    normalizedCpm <- cpm(normFactors)

    # write the data to local dir

    newFileName <- sub('_transposed.txt.gz', '', originalCountFile$properties$name)
    newFileName <- paste0(newFileName, "_normalized.txt", sep="")

    write.table(normalizedCpm, newFileName, quote = FALSE, sep = "\t", row.names = TRUE)

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
