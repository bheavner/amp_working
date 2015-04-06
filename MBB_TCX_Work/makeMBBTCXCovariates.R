# Goal: build MBB TCX covariates table.
# This will go in the folder at syn3163039

library(synapseClient) # for synapse upload
library(RCurl) # to grab google doc covariates files
require(gdata) # to read .xlsx file - http://www.r-bloggers.com/read-excel-files-from-r/

source("uploadToSynapse.R") # to annotate and upload file

#Login to Synapse using credentials saved in .synapseConfig file
synapseLogin()

# Grab existing info for inspection and assembly (James organized this across subdirectories in syn2924444)

#RNASeq_fullSampleInformation_NET.xlsx syn3163262
#mayo_tcx_rnaseq_clinical_vars.txt syn3163736
#mayo_tcx_rnaseq_tech_vars.txt syn3163738 
#mayo_tcx_rnaseq_sample_groups.txt syn3163739

# and I'll need to work with the samples that were re-run for higher readcounts
# https://docs.google.com/spreadsheets/d/1IN2Sgv517UMBnjYqcnIcv63KL-K42QCXl4P88XjwecY/

url1 <- getURL("https://docs.google.com/spreadsheets/d/1IN2Sgv517UMBnjYqcnIcv63KL-K42QCXl4P88XjwecY/export?format=csv")
reRunSamples <- read.csv(textConnection(url1))

metaDataFileSynapseID <- c('syn3163262', 'syn3163736', 'syn3163738', 'syn3163739')

for (metaDataFile in metaDataFileSynapseID) {
  message("Grabbing ", metaDataFile)
  
  # Download file from Synapse
  originalCountFile <- synGet(mergedCountFile)
  
  # unzip file and load for processing
  localFilePath <- getFileLocation(originalCountFile)
  
  if(!file.exists(substr(localFilePath, 1, nchar(localFilePath) - 3))) {
    gunzip(localFilePath)
  }
  
  localFilePath <- sub('.gz', '', localFilePath) #trim the .gz suffix
  
  Counts <- read.table(localFilePath, header = TRUE)
  





# better to use readxls

# excel summary files on google drive - permissions are set to "anyone with link can view"

# file “20141209 Taner mouse-mRNA Summary.xlsx”  includes Mouse_ID, Seq.Run.ID, Lane.Number, Clusters, RLIMS.ID; part of the Raw_RNAseq_file_name (based on Seq Run ID); it's at 
url1 <- getURL("https://docs.google.com/spreadsheets/d/1IQnJheILYLsUgWbwPT7Qmrmi-JngK0R3HTAJUfBYRUg/export?format=csv")

# “Copy of APP Randomization 9-16-14.xlsx” includes mouse_ID, Line == Experiment, RIN, Genotype (to transform), Sex, age in months; it's at
url2 <- getURL("https://docs.google.com/spreadsheets/d/1X_QRh-xw8q3lZ8IOosHdydxUsoP2Ucgxbfbk0FzMJx4/export?format=csv")

rawCovariates1 <- read.csv(textConnection(url1))
rawCovariates2 <- read.csv(textConnection(url2))

#remove first row ov rawCovariates2 (all blanks)
rawCovariates2 <- rawCovariates2[2:129,]

#rawCovariates1 has 188 observations. Those with rawCovariates1$Request.ID == "REQ-000000002442" overlaps the 128 observations in rawCovariates2. I'm not sure what the others are.

## Assemble APP covariates data frame
APPCovariates = data.frame("Mouse_ID" = rawCovariates2$Mouse.ID)  
#will need to check that these IDs agree with IDs in readcount data

APPCovariates <- cbind(APPCovariates, "Experiment" = rawCovariates2$Line)
APPCovariates <- cbind(APPCovariates, "RIN" = rawCovariates2$RIN)
APPCovariates <- cbind(APPCovariates, "Genotype" = rawCovariates2$Genotype)
APPCovariates <- cbind(APPCovariates, "Sex" = rawCovariates2$Sex)
APPCovariates <- cbind(APPCovariates, "Age_months" = rawCovariates2$Age)

#note: everything so far is from rawCovariates2, so when getting something from rawCovariates1, be sure it's in the right order
APPCovariates <- cbind(APPCovariates, "RLIMS.ID" = 
                         rawCovariates1$RLIMS.ID[match(APPCovariates$Mouse_ID, 
                                                       rawCovariates1$Sample.Name)])

APPCovariates <- cbind(APPCovariates, "Seq.Run.ID" = 
                         rawCovariates1$Seq.Run.ID[match(APPCovariates$Mouse_ID, 
                                                         rawCovariates1$Sample.Name)])
APPCovariates <- cbind(APPCovariates, "Lane.Number" = 
                         rawCovariates1$Lane.Number[match(APPCovariates$Mouse_ID, 
                                                          rawCovariates1$Sample.Name)])
APPCovariates <- cbind(APPCovariates, "Clusters" = 
                         rawCovariates1$Clusters[match(APPCovariates$Mouse_ID, 
                                                       rawCovariates1$Sample.Name)])

#TODO: filenames - do I need them? These are fastq files, so theree are two per sample - R1 and R2...
#APPCovariates <- cbind(APPCovariates, "Raw_RNAseq_file_name")


## Reformat values for consistency Tau covariates
APPCovariates$Genotype <- sub("NTG", "-", APPCovariates$Genotype)
APPCovariates$Genotype <- sub("NonTg", "-", APPCovariates$Genotype)
APPCovariates$Genotype <- sub("TG", "+", APPCovariates$Genotype)
APPCovariates$Genotype <- sub("Tg", "+", APPCovariates$Genotype)

APPCovariates$Sex <- sub("Female", "F", APPCovariates$Sex)
APPCovariates$Sex <- sub("Male", "M", APPCovariates$Sex)

APPCovariates$Age_months <- sub("mo", "", APPCovariates$Age_months)

# Reformat/rename IDs for consistency with readcount data
# First, replace space and # with _
APPCovariates$Mouse_ID <- sub("[#]", "_", APPCovariates$Mouse_ID)
APPCovariates$Mouse_ID <- sub("[ ]", "_", APPCovariates$Mouse_ID)

# next, replace - with . in all cases
APPCovariates$Mouse_ID <- sub('-', '.', APPCovariates$Mouse_ID)

# finally, I want to prepend all IDs that start with a number with "X". My regex-fu is weak:
APPCovariates$Mouse_ID <- sub('^1', 'X1', APPCovariates$Mouse_ID)
APPCovariates$Mouse_ID <- sub('^2', 'X2', APPCovariates$Mouse_ID)
APPCovariates$Mouse_ID <- sub('^3', 'X3', APPCovariates$Mouse_ID)
APPCovariates$Mouse_ID <- sub('^4', 'X4', APPCovariates$Mouse_ID)

# Write covariates table to local file
fileName <- "AMP-AD_TAUAPPms_UFL-Mayo-ISB_IlluminaHiSeq2000_App-Covariates.csv"
write.table(APPCovariates, file = fileName, quote = FALSE, sep = ",")

# upload to synapse with annotation, provenance, etc.
APPCovariatesId <- "syn3483880"
uploadToSynapse(fileName, APPCovariatesId)