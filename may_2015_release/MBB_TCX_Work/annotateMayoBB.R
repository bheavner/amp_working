# goal: add required metadata, filenames, and provenance to: 

# syn3163039    			folder
# syn3817650  				AMP-AD_MayoBB_UFL-Mayo-ISB_IlluminaHiSeq2000_TCX_Covariates.csv
# syn3818244					AMP-AD_MayoBB_UFL-Mayo-ISB_IlluminaHiSeq2000_TCX_GeneCounts.txt.gz
# syn3818429					AMP-AD_MayoBB_UFL-Mayo-ISB_IlluminaHiSeq2000_TCX_GeneCounts_Normalized.txt.gz
# syn3818611					AMP-AD_MayoBB_UFL-Mayo-ISB_IlluminaHiSeq2000_TCX_TranscriptCounts.txt.gz
# syn3818905					AMP-AD_MayoBB_UFL-Mayo-ISB_IlluminaHiSeq2000_TCX_TranscriptCounts_Normalized.txt.gz

library(synapseClient) # for synapse data exchange
library(tools) # for file path changes

#Login to Synapse using credentials saved in .synapseConfig file
synapseLogin()

### folder ###
synID <- "syn3163039"

# get info about current synapse object to modify
folderObject <- synGet(synID, downloadFile=F)

# define annotations as a list
folderAnnotations <- list(consortium = "AMP-AD",
                          study = "MayoRNAseq",
                          center = "UFL-Mayo-ISB",
                          disease = "Alzheimers Disease",
                          platform = "IlluminaHiSeq2000",
                          tissueTypeAbrv = "TCX",
                          tissueType = "Temporal Cortex",
                          organism = "Homo sapiens",
                          dataContact = "Ben Heavner ben.heavner@systemsbiology.org")
synSetAnnotations(folderObject) <- folderAnnotations

folderObject <- synStore(folderObject, forceVersion = T)

### covariates ###
file <- "AMP-AD_MayoBB_UFL-Mayo-ISB_IlluminaHiSeq2000_TCX_Covariates.csv"
synID <- "syn3817650"

# get info about current synapse object to modify
covariatesObject <- synGet(synID, downloadFile=F)

# define annotations as a list
covariatesObjectAnnotations <- list(consortium = "AMP-AD",
                                    study = "MayoRNAseq",
                                    center = "UFL-Mayo-ISB",
                                    disease = "Alzheimers Disease",
                                    dataType = "Covariates",
                                    platform = "IlluminaHiSeq2000",
                                    tissueTypeAbrv = "TCX",
                                    tissueType = "Temporal Cortex",
                                    organism = "Homo sapiens",
                                    fileType = "csv",
                                    dataContact = "Ben Heavner ben.heavner@systemsbiology.org")

synSetAnnotations(covariatesObject) <- covariatesObjectAnnotations

# no activity

covariatesObject <- synStore(covariatesObject, forceVersion = T)

##############
file <- "AMP-AD_MayoBB_UFL-Mayo-ISB_IlluminaHiSeq2000_TCX_GeneCounts.txt.gz"
synID <- "syn3818244"

# get info about current synapse object to modify
geneCountObject <- synGet(synID, downloadFile=F)

# define annotations as a list
geneCountObjectAnnotations <- list(consortium = "AMP-AD",
                                   study = "MayoRNAseq",
                                   center = "UFL-Mayo-ISB",
                                   disease = "Alzheimers Disease",
                                   dataType = "mRNA",
                                   platform = "IlluminaHiSeq2000",
                                   tissueTypeAbrv = "TCX",
                                   tissueType = "Temporal Cortex",
                                   organism = "Homo sapiens",
                                   fileType = "count",
                                   dataContact = "Ben Heavner ben.heavner@systemsbiology.org")

synSetAnnotations(geneCountObject) <- geneCountObjectAnnotations

# define activity
codeFile <- list("https://github.com/PriceLab/AMP-ad/blob/v0.4/MBB/TCX/merge_count_files.R",
                 "https://github.com/PriceLab/AMP-ad/blob/v0.4/MBB/TCX/merge_MayoTCX_rnaseq_counts.R")   

act <- Activity(name='Merge output count files from SNAPR',
                used = as.list('syn3632530'),
                executed = as.list(codeFile))

generatedBy(geneCountObject) = act

geneCountObject <- synStore(geneCountObject, forceVersion = T)

###############
file <- "AMP-AD_MayoBB_UFL-Mayo-ISB_IlluminaHiSeq2000_TCX_GeneCounts_Normalized.txt.gz"
synID <- "syn3818429"

# get info about current synapse object to modify
normalizedGeneCountObject <- synGet(synID, downloadFile=F)

# define annotations as a list
normalizedGeneCountObjectAnnotations <- list(consortium = "AMP-AD",
                                             study = "MayoRNAseq",
                                             center = "UFL-Mayo-ISB",
                                             disease = "Alzheimers Disease",
                                             dataType = "mRNA",
                                             platform = "IlluminaHiSeq2000",
                                             tissueTypeAbrv = "TCX",
                                             tissueType = "Temporal Cortex",
                                             organism = "Homo sapiens",
                                             fileType = "count",
                                             dataContact = "Ben Heavner ben.heavner@systemsbiology.org")

synSetAnnotations(normalizedGeneCountObject) <- normalizedGeneCountObjectAnnotations

# define activity
codeFile <- ("https://github.com/PriceLab/AMP-ad/blob/v0.4/MBB/TCX/normalize_readcounts.R")   

act <- Activity(name='CPM (using TMM) from edgeR normalization',
                used = as.list('syn3818244'),
                executed = as.list(codeFile))

generatedBy(normalizedGeneCountObject) = act

normalizedGeneCountObject <- synStore(normalizedGeneCountObject, forceVersion = T)




##############
file <- "AMP-AD_MayoBB_UFL-Mayo-ISB_IlluminaHiSeq2000_TCX_TranscriptCounts.txt.gz"
synID <- "syn3818611"

# get info about current synapse object to modify
transcriptCountObject <- synGet(synID, downloadFile=F)

# define annotations as a list
transcriptCountObjectAnnotations <- list(consortium = "AMP-AD",
                                         study = "MayoRNAseq",
                                         center = "UFL-Mayo-ISB",
                                         disease = "Alzheimers Disease",
                                         dataType = "mRNA",
                                         platform = "IlluminaHiSeq2000",
                                         tissueTypeAbrv = "TCX",
                                         tissueType = "Temporal Cortex",
                                         organism = "Homo sapiens",
                                         fileType = "count",
                                         dataContact = "Ben Heavner ben.heavner@systemsbiology.org")

synSetAnnotations(transcriptCountObject) <- transcriptCountObjectAnnotations

# define activity
codeFile <- list("https://github.com/PriceLab/AMP-ad/blob/v0.4/MBB/TCX/merge_count_files.R",
                 "https://github.com/PriceLab/AMP-ad/blob/v0.4/MBB/TCX/merge_MayoTCX_rnaseq_counts.R") 

act <- Activity(name='Merge output count files from SNAPR',
                used = as.list('syn3632530'),
                executed = as.list(codeFile))

generatedBy(transcriptCountObject) = act

transcriptCountObject <- synStore(transcriptCountObject, forceVersion = T)

###############
file <- "AMP-AD_MayoBB_UFL-Mayo-ISB_IlluminaHiSeq2000_TCX_TranscriptCounts_Normalized.txt"
synID <- "syn3818905"

# get info about current synapse object to modify
normalizedTranscriptCountObject <- synGet(synID, downloadFile=F)

# define annotations as a list
normalizedTranscriptCountObjectAnnotations <- list(consortium = "AMP-AD",
                                                   study = "MayoRNAseq",
                                                   center = "UFL-Mayo-ISB",
                                                   disease = "Alzheimers Disease",
                                                   dataType = "mRNA",
                                                   platform = "IlluminaHiSeq2000",
                                                   tissueTypeAbrv = "TCX",
                                                   tissueType = "Temporal Cortex",
                                                   organism = "Homo sapiens",
                                                   fileType = "count",
                                                   dataContact = "Ben Heavner ben.heavner@systemsbiology.org")

synSetAnnotations(normalizedTranscriptCountObject) <- normalizedTranscriptCountObjectAnnotations

# define activity
codeFile <- ("https://github.com/PriceLab/AMP-ad/blob/v0.4/MBB/TCX/normalize_readcounts.R")   

act <- Activity(name='CPM (using TMM) from edgeR normalization',
                used = as.list('syn3818611'),
                executed = as.list(codeFile))

generatedBy(normalizedTranscriptCountObject) = act

normalizedTranscriptCountObject <- synStore(normalizedTranscriptCountObject, forceVersion = T)

