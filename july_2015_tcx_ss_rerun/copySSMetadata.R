# This script copies metadata and provenance from older synapse objects

# syn4719606 - counts by gene ID was syn3801377
# syn4719692 - normalized counts by gene ID was syn3801467
# syn4719609 - counts by transcript ID was syn3801548
# syn4719707 - normalized counts by transcript ID was syn3801636

library(synapseClient)

#Login to Synapse using credentials saved in .synapseConfig file
synapseLogin()

# metatata for gene_ID
oldGeneID <- "syn3801377"
newGeneID <- "syn4719606"

oldGene <- synGet(oldGeneID, downloadFile=F)
newGene <- synGet(newGeneID, downloadFile=F)

oldAnnotations <- as.list(synGetAnnotations(oldGene))
synSetAnnotations(newGene) <- oldAnnotations

# inspect old activity
synGetActivity(oldGene)

# define new activity
codeFile <- list("https://github.com/PriceLab/AMP-ad/blob/v0.3/SS/RushBroad/merge_count_files.R",
                 "https://github.com/PriceLab/AMP-ad/blob/v0.3/SS/RushBroad/merge_RushBroadSS_rnaseq_counts.R")   

act <- Activity(name='Merge output count files from SNAPR',
                used = as.list('syn3619668'),
                executed = as.list(codeFile))

generatedBy(newGene) = act

synStore(newGene, forceVersion = T)

## clear environment here for safety
rm(list = ls())

# metatata for normalized gene_ID
oldNormGeneID <- "syn3801467"
newNormGeneID <- "syn4719692"

oldNormGene <- synGet(oldNormGeneID, downloadFile=F)
newNormGene <- synGet(newNormGeneID, downloadFile=F)

oldNormAnnotations <- as.list(synGetAnnotations(oldNormGene))
synSetAnnotations(newNormGene) <- oldNormAnnotations

# inspect old activity
synGetActivity(oldNormGene)

# define new activity
codeFile <- list("https://github.com/PriceLab/AMP-ad/blob/v0.3/SS/RushBroad/normalize_readcounts.R")   

act <- Activity(name='CPM (using TMM) from edgeR normalization',
                used = as.list('syn4719606'),
                executed = as.list(codeFile))

generatedBy(newNormGene) = act

synStore(newNormGene, forceVersion = T)


## clear environment here for safety
rm(list = ls())

# metatata for transcript_ID
oldTranscriptID <- "syn3801548"
newTranscriptID <- "syn4719609"

oldTranscript <- synGet(oldTranscriptID, downloadFile=F)
newTranscript <- synGet(newTranscriptID, downloadFile=F)

oldAnnotations <- as.list(synGetAnnotations(oldTranscript))
synSetAnnotations(newTranscript) <- oldAnnotations

# inspect old activity
synGetActivity(oldTranscript)

# define new activity
codeFile <- list("https://github.com/PriceLab/AMP-ad/blob/v0.3/SS/RushBroad/merge_count_files.R",
                 "https://github.com/PriceLab/AMP-ad/blob/v0.3/SS/RushBroad/merge_RushBroadSS_rnaseq_counts.R")   

act <- Activity(name='Merge output count files from SNAPR',
                used = as.list('syn3619668'),
                executed = as.list(codeFile))

generatedBy(newTranscript) = act

synStore(newTranscript, forceVersion = T)

## clear environment here for safety
rm(list = ls())

# metatata for normalized gene_ID
# syn4719707 - normalized counts by transcript ID was syn3801636
oldNormTranscriptID <- "syn3801636"
newNormTranscriptID <- "syn4719707"

oldNormTranscript <- synGet(oldNormTranscriptID, downloadFile=F)
newNormTranscript <- synGet(newNormTranscriptID, downloadFile=F)

oldNormAnnotations <- as.list(synGetAnnotations(oldNormTranscript))
synSetAnnotations(newNormTranscript) <- oldNormAnnotations

# inspect old activity
synGetActivity(oldNormTranscript)

# define new activity
codeFile <- list("https://github.com/PriceLab/AMP-ad/blob/v0.3/SS/RushBroad/normalize_readcounts.R")   

act <- Activity(name='CPM (using TMM) from edgeR normalization',
                used = as.list('syn4719609'),
                executed = as.list(codeFile))

generatedBy(newNormTranscript) = act

synStore(newNormTranscript, forceVersion = T)
