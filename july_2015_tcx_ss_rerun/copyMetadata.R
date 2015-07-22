# This script copies metadata and provenance from older synapse objects

# syn4650257 - counts by gene ID was syn3818244
# syn4650265 - normalized counts by gene ID was syn3818429
# syn4650258 - counts by transcript ID was syn3818611
# syn4650430 - normalized counts by transcript ID was syn3818905

library(synapseClient)

#Login to Synapse using credentials saved in .synapseConfig file
synapseLogin()

# metatata for gene_ID
oldGeneID <- "syn3818244"
newGeneID <- "syn4650257"

oldGene <- synGet(oldGeneID, downloadFile=F)
newGene <- synGet(newGeneID, downloadFile=F)

oldAnnotations <- as.list(synGetAnnotations(oldGene))
synSetAnnotations(newGene) <- oldAnnotations

# inspect old activity
synGetActivity(oldGene)

# define new activity
codeFile <- list("https://github.com/PriceLab/AMP-ad/blob/v0.4/MBB/TCX/merge_MayoTCX_rnaseq_counts.R",
                 "https://github.com/PriceLab/AMP-ad/blob/v0.4/MBB/TCX/merge_count_files.R")   

act <- Activity(name='Merge output count files from SNAPR',
                used = as.list('syn3632530'),
                executed = as.list(codeFile))

generatedBy(newGene) = act

synStore(newGene, forceVersion = T)

## clear environment here for safety
rm(list = ls())

# metatata for normalized gene_ID
oldNormGeneID <- "syn3818429"
newNormGeneID <- "syn4650265"

oldNormGene <- synGet(oldNormGeneID, downloadFile=F)
newNormGene <- synGet(newNormGeneID, downloadFile=F)

oldNormAnnotations <- as.list(synGetAnnotations(oldNormGene))
synSetAnnotations(newNormGene) <- oldNormAnnotations

# inspect old activity
synGetActivity(oldNormGene)

# define new activity
codeFile <- list("https://github.com/PriceLab/AMP-ad/blob/v0.4/MBB/TCX/normalize_readcounts.R")   

act <- Activity(name='CPM (using TMM) from edgeR normalization',
                used = as.list('syn4650257'),
                executed = as.list(codeFile))

generatedBy(newNormGene) = act

synStore(newNormGene, forceVersion = T)


## clear environment here for safety
rm(list = ls())

# metatata for transcript_ID
oldTranscriptID <- "syn3818611"
newTranscriptID <- "syn4650258"

oldTranscript <- synGet(oldTranscriptID, downloadFile=F)
newTranscript <- synGet(newTranscriptID, downloadFile=F)

oldAnnotations <- as.list(synGetAnnotations(oldTranscript))
synSetAnnotations(newTranscript) <- oldAnnotations

# inspect old activity
synGetActivity(oldTranscript)

# define new activity
codeFile <- list("https://github.com/PriceLab/AMP-ad/blob/v0.4/MBB/TCX/merge_MayoTCX_rnaseq_counts.R",
                 "https://github.com/PriceLab/AMP-ad/blob/v0.4/MBB/TCX/merge_count_files.R")   

act <- Activity(name='Merge output count files from SNAPR',
                used = as.list('syn3632530'),
                executed = as.list(codeFile))

generatedBy(newTranscript) = act

synStore(newTranscript, forceVersion = T)

## clear environment here for safety
rm(list = ls())

# metatata for normalized gene_ID
# syn4650430 - normalized counts by transcript ID was syn3818905
oldNormTranscriptID <- "syn3818905"
newNormTranscriptID <- "syn4650430"

oldNormTranscript <- synGet(oldNormTranscriptID, downloadFile=F)
newNormTranscript <- synGet(newNormTranscriptID, downloadFile=F)

oldNormAnnotations <- as.list(synGetAnnotations(oldNormTranscript))
synSetAnnotations(newNormTranscript) <- oldNormAnnotations

# inspect old activity
synGetActivity(oldNormTranscript)

# define new activity
codeFile <- list("https://github.com/PriceLab/AMP-ad/blob/v0.4/MBB/TCX/normalize_readcounts.R")   

act <- Activity(name='CPM (using TMM) from edgeR normalization',
                used = as.list('syn4650258'),
                executed = as.list(codeFile))

generatedBy(newNormTranscript) = act

synStore(newNormTranscript, forceVersion = T)
