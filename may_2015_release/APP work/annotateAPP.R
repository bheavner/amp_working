# goal: add required metadata, filenames, and provenance to: 

# syn3435792    			folder
# syn3483880    			AMP-AD_TAUAPPms_UFL-Mayo_ISB_IlluminaHiSeq2000_App_Covariates.csv

# syn3483519  				AMP-AD_TAUAPPms_UFL-Mayo_ISB_IlluminaHiSeq2000_App_GeneCounts.txt.gz
# syn3505856					AMP-AD_TAUAPPms_UFL-Mayo_ISB_IlluminaHiSeq2000_App_GeneCounts_Normalized.txt
# syn3483625					AMP-AD_TAUAPPms_UFL-Mayo_ISB_IlluminaHiSeq2000_App_TranscriptCounts.txt.gz
# syn3505860					AMP-AD_TAUAPPms_UFL-Mayo_ISB_IlluminaHiSeq2000_App_TranscriptCounts_Normalized.txt
# syn3439202    			tarred Snapr output (source, not for public release)

library(synapseClient) # for synapse data exchange
library(tools) # for file path changes

#Login to Synapse using credentials saved in .synapseConfig file
synapseLogin()

### folder ###
synID <- "syn3435792"

# get info about current synapse object to modify
folderObject <- synGet(synID, downloadFile=F)

# define annotations as a list
folderAnnotations <- list(consortium = "AMP-AD",
                          study = "TAUAPPms",
                          center = "UFL-Mayo-ISB",
                          disease = "Alzheimers Disease",
                          platform = "IlluminaHiSeq2000",
                          mouseModel = "APP",
                          tissueTypeAbrv = "FB",
                          tissueType = "Forebrain",
                          organism = "Mus musculus",
                          dataContact = "Ben Heavner ben.heavner@systemsbiology.org")
synSetAnnotations(folderObject) <- folderAnnotations

folderObject <- synStore(folderObject, forceVersion = T)

### covariates ###
file <- "AMP-AD_TAUAPPms_UFL-Mayo_ISB_IlluminaHiSeq2000_App_Covariates.csv"
synID <- "syn3483880"

# get info about current synapse object to modify
covariatesObject <- synGet(synID, downloadFile=F)

# define annotations as a list
covariatesObjectAnnotations <- list(consortium = "AMP-AD",
                                    study = "TAUAPPms",
                                    center = "UFL-Mayo-ISB",
                                    dataType = "Covariates",
                                    disease = "Alzheimers Disease",
                                    platform = "IlluminaHiSeq2000",
                                    mouseModel = "APP",
                                    tissueTypeAbrv = "FB",
                                    tissueType = "Forebrain",
                                    fileType = "csv",
                                    organism = "Mus musculus",
                                    dataContact = "Ben Heavner ben.heavner@systemsbiology.org")
synSetAnnotations(covariatesObject) <- covariatesObjectAnnotations

# define activity
codeFile <- ("https://github.com/PriceLab/AMP-ad/blob/v0.2/UFL/APP/makeAPPCovariates.R")   

act <- Activity(name='Covariate file generation',
                executed = as.list(codeFile))

generatedBy(covariatesObject) = act

covariatesObject <- synStore(covariatesObject, forceVersion = T)

##############
file <- "AMP-AD_TAUAPPms_UFL-Mayo_ISB_IlluminaHiSeq2000_App_GeneCounts.txt.gz"
synID <- "syn3483519"

# get info about current synapse object to modify
geneCountObject <- synGet(synID, downloadFile=F)

# define annotations as a list
geneCountObjectAnnotations <- list(consortium = "AMP-AD",
                            study = "TAUAPPms",
                            center = "UFL-Mayo-ISB",
                            dataType = "mRNA",
                            disease = "Alzheimers Disease",
                            platform = "IlluminaHiSeq2000",
                            mouseModel = "APP",
                            tissueTypeAbrv = "FB",
                            tissueType = "Forebrain",
                            fileType = "count",
                            organism = "Mus musculus",
                            dataContact = "Ben Heavner ben.heavner@systemsbiology.org")
synSetAnnotations(geneCountObject) <- geneCountObjectAnnotations

# define activity
codeFile <- list("https://github.com/PriceLab/AMP-ad/blob/v0.2/UFL/APP/merge_count_files.R",
                 "https://github.com/PriceLab/AMP-ad/blob/v0.2/UFL/APP/merge_mouse_app_rnaseq_counts.R")   

act <- Activity(name='Merge output count files from SNAPR',
                used = as.list('syn3439202'),
                executed = as.list(codeFile))

generatedBy(geneCountObject) = act

geneCountObject <- synStore(geneCountObject, forceVersion = T)

###############
file <- "AMP-AD_TAUAPPms_UFL-Mayo_ISB_IlluminaHiSeq2000_App_GeneCounts_Normalized.txt"
synID <- "syn3505856"

# get info about current synapse object to modify
normalizedGeneCountObject <- synGet(synID, downloadFile=F)

# define annotations as a list
normalizedGeneCountObjectAnnotations <- list(consortium = "AMP-AD",
                                   study = "TAUAPPms",
                                   center = "UFL-Mayo-ISB",
                                   dataType = "mRNA",
                                   disease = "Alzheimers Disease",
                                   platform = "IlluminaHiSeq2000",
                                   mouseModel = "APP",
                                   tissueTypeAbrv = "FB",
                                   tissueType = "Forebrain",
                                   fileType = "count",
                                   organism = "Mus musculus",
                                   dataContact = "Ben Heavner ben.heavner@systemsbiology.org")
synSetAnnotations(normalizedGeneCountObject) <- normalizedGeneCountObjectAnnotations

# define activity
codeFile <- ("https://github.com/PriceLab/AMP-ad/blob/v0.2/UFL/APP/normalize_readcounts.R")   

act <- Activity(name='Merge output count files from SNAPR',
                used = as.list('syn3483519'),
                executed = as.list(codeFile))

generatedBy(normalizedGeneCountObject) = act

normalizedGeneCountObject <- synStore(normalizedGeneCountObject, forceVersion = T)




##############
file <- "AMP-AD_TAUAPPms_UFL-Mayo_ISB_IlluminaHiSeq2000_App_TranscriptCounts.txt.gz"
synID <- "syn3483625"

# get info about current synapse object to modify
transcriptCountObject <- synGet(synID, downloadFile=F)

# define annotations as a list
transcriptCountObjectAnnotations <- list(consortium = "AMP-AD",
                                   study = "TAUAPPms",
                                   center = "UFL-Mayo-ISB",
                                   dataType = "mRNA",
                                   disease = "Alzheimers Disease",
                                   platform = "IlluminaHiSeq2000",
                                   mouseModel = "APP",
                                   tissueTypeAbrv = "FB",
                                   tissueType = "Forebrain",
                                   fileType = "count",
                                   organism = "Mus musculus",
                                   dataContact = "Ben Heavner ben.heavner@systemsbiology.org")
synSetAnnotations(transcriptCountObject) <- transcriptCountObjectAnnotations

# define activity
codeFile <- list("https://github.com/PriceLab/AMP-ad/blob/v0.2/UFL/APP/merge_count_files.R",
                 "https://github.com/PriceLab/AMP-ad/blob/v0.2/UFL/APP/merge_mouse_app_rnaseq_counts.R")  

act <- Activity(name='Merge output count files from SNAPR',
                used = as.list('syn3439202'),
                executed = as.list(codeFile))

generatedBy(transcriptCountObject) = act

transcriptCountObject <- synStore(transcriptCountObject, forceVersion = T)

###############
file <- "AMP-AD_TAUAPPms_UFL-Mayo_ISB_IlluminaHiSeq2000_App_TranscriptCounts_Normalized.txt"
synID <- "syn3505860"

# get info about current synapse object to modify
normalizedTranscriptCountObject <- synGet(synID, downloadFile=F)

# define annotations as a list
normalizedTranscriptCountObjectAnnotations <- list(consortium = "AMP-AD",
                                             study = "TAUAPPms",
                                             center = "UFL-Mayo-ISB",
                                             dataType = "mRNA",
                                             disease = "Alzheimers Disease",
                                             platform = "IlluminaHiSeq2000",
                                             mouseModel = "APP",
                                             tissueTypeAbrv = "FB",
                                             tissueType = "Forebrain",
                                             fileType = "count",
                                             organism = "Mus musculus",
                                             dataContact = "Ben Heavner ben.heavner@systemsbiology.org")
synSetAnnotations(normalizedTranscriptCountObject) <- normalizedTranscriptCountObjectAnnotations

# define activity
codeFile <- ("https://github.com/PriceLab/AMP-ad/blob/v0.2/UFL/APP/normalize_readcounts.R")   

act <- Activity(name='Merge output count files from SNAPR',
                used = as.list('syn3483519'),
                executed = as.list(codeFile))

generatedBy(normalizedTranscriptCountObject) = act

normalizedTranscriptCountObject <- synStore(normalizedTranscriptCountObject, forceVersion = T)

