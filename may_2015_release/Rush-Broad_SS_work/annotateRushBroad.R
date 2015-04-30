# goal: add required metadata, filenames, and provenance to: 

# syn3537579    			folder
# syn3630140  				AMP-AD_SampleSwap_UFL-Mayo-ISB_IlluminaHiSeq2000_RushBroadSS_Covariates.csv
# syn3801377					AMP-AD_SampleSwap_UFL-Mayo-ISB_IlluminaHiSeq2000_dIPFC_Rush-Broad-SS_GeneCounts.txt.gz
# syn3801467					AMP-AD_SampleSwap_UFL-Mayo-ISB_IlluminaHiSeq2000_dIPFC_Rush-Broad-SS_geneCounts_normalized.txt.gz
# syn3801548					AMP-AD_SampleSwap_UFL-Mayo-ISB_IlluminaHiSeq2000_dIPFC_Rush-Broad-SS_TranscriptCounts.txt.gz
# syn3801636					AMP-AD_SampleSwap_UFL-Mayo-ISB_IlluminaHiSeq2000_dIPFC_Rush-Broad-SS_TranscriptCounts_Normalized.txt.gz

library(synapseClient) # for synapse data exchange
library(tools) # for file path changes

#Login to Synapse using credentials saved in .synapseConfig file
synapseLogin()

### folder ###
synID <- "syn3537579"

# get info about current synapse object to modify
folderObject <- synGet(synID, downloadFile=F)

# define annotations as a list
folderAnnotations <- list(consortium = "AMP-AD",
                          study = list("ROSMAP", "SampleSwap"),
                          center = list("Broad-Rush", "UFL_Mayo-ISB"),
                          disease = "Alzheimers Disease",
                          platform = "IlluminaHiSeq2000",
                          tissueTypeAbrv = "PFC",
                          tissueType = "Dorsolateral Prefrontal Cortex",
                          organism = "  Homo sapiens",
                          dataContact = "Ben Heavner ben.heavner@systemsbiology.org")
synSetAnnotations(folderObject) <- folderAnnotations

folderObject <- synStore(folderObject, forceVersion = T)

### covariates ###
synID <- "syn3630140"

# get info about current synapse object to modify
covariatesObject <- synGet(synID, downloadFile=F)

# define annotations as a list
covariatesObjectAnnotations <- list(consortium = "AMP-AD",
                                    study = list("ROSMAP", "SampleSwap"),
                                    center = list("Broad-Rush", "UFL_Mayo-ISB"),
                                    disease = "Alzheimers Disease",
                                    platform = "IlluminaHiSeq2000",
                                    tissueTypeAbrv = "PFC",
                                    tissueType = "Dorsolateral Prefrontal Cortex",
                                    organism = "  Homo sapiens",
                                    dataContact = "Ben Heavner ben.heavner@systemsbiology.org",
                                    dataType = "Covariates",
                                    fileType = "csv")
synSetAnnotations(covariatesObject) <- covariatesObjectAnnotations

# no activity

covariatesObject <- synStore(covariatesObject, forceVersion = T)

## genecounts ##
synID <- "syn3801377"

# get info about current synapse object to modify
geneCountObject <- synGet(synID, downloadFile=F)

# define annotations as a list
geneCountObjectAnnotations <- list(consortium = "AMP-AD",
                                   study = list("ROSMAP", "SampleSwap"),
                                   center = list("Broad-Rush", "UFL_Mayo-ISB"),
                                   disease = "Alzheimers Disease",
                                   platform = "IlluminaHiSeq2000",
                                   tissueTypeAbrv = "PFC",
                                   tissueType = "Dorsolateral Prefrontal Cortex",
                                   organism = "  Homo sapiens",
                                   dataContact = "Ben Heavner ben.heavner@systemsbiology.org",
                                   dataType = "mRNA",
                                   fileType = "count")

synSetAnnotations(geneCountObject) <- geneCountObjectAnnotations

# define activity
codeFile <- list("https://github.com/PriceLab/AMP-ad/blob/v0.3/SS/RushBroad/merge_count_files.R",
                 "https://github.com/PriceLab/AMP-ad/blob/v0.3/SS/RushBroad/merge_RushBroadSS_rnaseq_counts.R")   

act <- Activity(name='Merge output count files from SNAPR',
                used = as.list('syn3619668'),
                executed = as.list(codeFile))

generatedBy(geneCountObject) = act

geneCountObject <- synStore(geneCountObject, forceVersion = T)

###############
## normalized genecounts ##
synID <- "syn3801467"

# get info about current synapse object to modify
normalizedGeneCountObject <- synGet(synID, downloadFile=F)

# define annotations as a list
normalizedGeneCountObjectAnnotations <- list(consortium = "AMP-AD",
                                             study = list("ROSMAP", "SampleSwap"),
                                             center = list("Broad-Rush", "UFL_Mayo-ISB"),
                                             disease = "Alzheimers Disease",
                                             platform = "IlluminaHiSeq2000",
                                             tissueTypeAbrv = "PFC",
                                             tissueType = "Dorsolateral Prefrontal Cortex",
                                             organism = "  Homo sapiens",
                                             dataContact = "Ben Heavner ben.heavner@systemsbiology.org",
                                             dataType = "mRNA",
                                             fileType = "count")
synSetAnnotations(normalizedGeneCountObject) <- normalizedGeneCountObjectAnnotations

# define activity
codeFile <- ("https://github.com/PriceLab/AMP-ad/blob/v0.3/SS/RushBroad/normalize_readcounts.R")   

act <- Activity(name='Normalize readcounts',
                used = as.list('syn3801377'),
                executed = as.list(codeFile))

generatedBy(normalizedGeneCountObject) = act

normalizedGeneCountObject <- synStore(normalizedGeneCountObject, forceVersion = T)




##############
## TranscriptCounts ##
synID <- "syn3801548"

# get info about current synapse object to modify
transcriptCountObject <- synGet(synID, downloadFile=F)

# define annotations as a list
transcriptCountObjectAnnotations <- list(consortium = "AMP-AD",
                                   study = list("ROSMAP", "SampleSwap"),
                                   center = list("Broad-Rush", "UFL_Mayo-ISB"),
                                   disease = "Alzheimers Disease",
                                   platform = "IlluminaHiSeq2000",
                                   tissueTypeAbrv = "PFC",
                                   tissueType = "Dorsolateral Prefrontal Cortex",
                                   organism = "  Homo sapiens",
                                   dataContact = "Ben Heavner ben.heavner@systemsbiology.org",
                                   dataType = "mRNA",
                                   fileType = "count")
synSetAnnotations(transcriptCountObject) <- transcriptCountObjectAnnotations

# define activity
codeFile <- list("https://github.com/PriceLab/AMP-ad/blob/v0.3/SS/RushBroad/merge_count_files.R",
                 "https://github.com/PriceLab/AMP-ad/blob/v0.3/SS/RushBroad/merge_RushBroadSS_rnaseq_counts.R")   

act <- Activity(name='Merge output count files from SNAPR',
                used = as.list('syn3619668'),
                executed = as.list(codeFile))

generatedBy(transcriptCountObject) = act

transcriptCountObject <- synStore(transcriptCountObject, forceVersion = T)

###############
## normalized transcript counts ##
synID <- "syn3801636"

# get info about current synapse object to modify
normalizedTranscriptCountObject <- synGet(synID, downloadFile=F)

# define annotations as a list
normalizedTranscriptCountObjectAnnotations <- list(consortium = "AMP-AD",
                                                   study = list("ROSMAP", "SampleSwap"),
                                                   center = list("Broad-Rush", "UFL_Mayo-ISB"),
                                                   disease = "Alzheimers Disease",
                                                   platform = "IlluminaHiSeq2000",
                                                   tissueTypeAbrv = "PFC",
                                                   tissueType = "Dorsolateral Prefrontal Cortex",
                                                   organism = "  Homo sapiens",
                                                   dataContact = "Ben Heavner ben.heavner@systemsbiology.org",
                                                   dataType = "mRNA",
                                                   fileType = "count")
synSetAnnotations(normalizedTranscriptCountObject) <- normalizedTranscriptCountObjectAnnotations

# define activity
# define activity
codeFile <- ("https://github.com/PriceLab/AMP-ad/blob/v0.3/SS/RushBroad/normalize_readcounts.R")   

act <- Activity(name='Normalize readcounts',
                used = as.list('syn3801548'),
                executed = as.list(codeFile))
generatedBy(normalizedTranscriptCountObject) = act

normalizedTranscriptCountObject <- synStore(normalizedTranscriptCountObject, forceVersion = T)

