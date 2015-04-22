# Goal: build sample swap covariates table. Per Dr. Younkin, "We ran RNASeq on ... 10 samples from rush/broad

# This will go in the folder at syn3388564

library(synapseClient) # for synapse upload
library(RCurl) # to grab google doc covariates files
library(gdata) # to read .xlsx file -  install.packages("gdata")

# source("uploadToSynapse.R") # to annotate and upload file

#Login to Synapse using credentials saved in .synapseConfig file
synapseLogin()

sampleSwapIds <- c("10249336", "50100518", "50104008", "20177982", "36492755",
                   "05689621", "11444465", "10315029", "20270920", "11615242")

# Grab existing info for inspection and assembly (James organized this across subdirectories in syn2924444)

#RNASeq_fullSampleInformation_NET.xlsx syn3163262
# has RNASubjectID, RNAId, Source, tissue, final diagnosis, RIN, unfiltered age, sex, and Braak
fullSampleInformation <- synGet("syn3163262")
fullSampleInformationFilePath <- getFileLocation(fullSampleInformation)
fullSample <- read.xls(fullSampleInformationFilePath, sheet = 1, header = TRUE, stringsAsFactors = FALSE)

#mayo_tcx_rnaseq_clinical_vars.txt syn3163736
# censored: has participant_id, age_at_last_assessment for some samples, sex, Braak
clinicalVars <- synGet("syn3163736")
clinicalVarsFilePath <- getFileLocation(clinicalVars)
clinical <- read.table(clinicalVarsFilePath, header = TRUE, stringsAsFactors = FALSE)

#mayo_tcx_rnaseq_tech_vars.txt syn3163738 
# has RNASubjectID, RNAId, Source, tissue, RIN
techVars <- synGet("syn3163738")
techVarsFilePath <- getFileLocation(techVars)
tech <- read.table(techVarsFilePath, header = TRUE, stringsAsFactors = FALSE)

#mayo_tcx_rnaseq_sample_groups.txt syn3163739
# may not need
sampleGroups <- synGet("syn3163739")
sampleGroupsFilePath <- getFileLocation(sampleGroups)
groups <- read.table(sampleGroupsFilePath, header = TRUE, stringsAsFactors = FALSE)
#Error in scan(file, what, nmax, sep, dec, quote, skip, nlines, na.strings,  : line 19 did not have 6 elements

#U01_288_AUT_TCx_RNAseq_Covars-Drives_02-06-2015_1447.xlsx syn3523879
# THIS IS WHAT I WANT TO USE FOR EVERYTHING.
# includes SubjectID, SampleID, seq Run ID, DriveA_FileName1, DriveA_FileName2, DriveB_FileName1, and DriveB_Filename2
rerunSamples  <- synGet("syn3523879")
rerunSamplesFilePath <- getFileLocation(rerunSamples)
rerun <- read.xls(rerunSamplesFilePath, sheet = 2, header = TRUE, stringsAsFactors = FALSE)

# the fullSampleInformation_NET.xlsx includes data from 3 sources: MayoBrainBank_Dickson, BannerSunHealth_TomBeach, and RUSH_BROAD. I need to generate different covariates files for the RUSH_BROAD and the other samples. RUSH_BROAD under syn3388564, others under syn3163039.

# Begin by building lists of sample IDs to include in each of the 2 covariates files (RUSH_BROAD and MBBTCX)

rushBroadCovariates  <- data.frame()
rushBroadCovariates$Samples <- fullSample[fullSample$Source == "RUSH-BROAD", "RNASubjectId"]



#remove first row of rawCovariates2 (all blanks)
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