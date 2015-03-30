# Goal: build APP covariates table with same info as Tau data. 
# NOTE: I expect 128 samples in this group
# This will go in the folder at syn3435792

# Log in to synapse and grab Tau covariates to compare
library(synapseClient)
library(RCurl) # to grab google doc covariates files

#Login to Synapse using credentials saved in .synapseConfig file
synapseLogin()

# get the normalized readcount and covariates files from synapse
covariatesFile <- synGet('syn2875343')

# load Tau covariates file to have handy for comparison
covariates <- read.table(getFileLocation(covariatesFile), 
                         header = TRUE, 
                         stringsAsFactors = FALSE)

colnames(covariates)
#  "Mouse_ID" "Experiment" "RIN" "Genotype" "Sex" "Age_months" "RLIMS.ID" "Seq.Run.ID" "Lane.Number" "Clusters" "Raw_RNAseq_file_name"

# excel summary files on google drive - permissions are set to "anyone with link can view"

# file “20141209 Taner mouse-mRNA Summary.xlsx”  includes Mouse_ID, Seq.Run.ID, Lane.Number, Clusters, RLIMS.ID; part of the Raw_RNAseq_file_name (based on Seq Run ID); it's at 
url1 <- getURL("https://docs.google.com/spreadsheets/d/1IQnJheILYLsUgWbwPT7Qmrmi-JngK0R3HTAJUfBYRUg/export?format=csv")

# “Copy of APP Randomization 9-16-14.xlsx” includes mouse_ID, Line == Experiment, RIN, Genotype (to transform), Sex, age in months; it's at
url2 <- getURL("https://docs.google.com/spreadsheets/d/1X_QRh-xw8q3lZ8IOosHdydxUsoP2Ucgxbfbk0FzMJx4/export?format=csv")

rawCovariates1 <- read.csv(textConnection(url1))
rawCovariates2 <- read.csv(textConnection(url2))
#remove first row (all blanks)
rawCovariates2 <- rawCovariates2[2:129,]

#rawCovariates1 has 188 observations. I think I want rawCovariates1$Request.ID == "REQ-000000002442" (that subset overlaps rawCovariates2)

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

#TODO:
## Reformat IDs for consistency with readcount data

# upload to synapse with correct file name, annotation, and provenance
# confirm wiki descriptions needed