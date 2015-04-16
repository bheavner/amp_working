# goal: build list of files to include in data release Per Dr. Younkin, "We ran RNASeq on 278 mayo temporal cortex samples plus 10 samples from rush/broad = 288 samples." NOTE: of the 288 RNAseq samples, 10 are “sample-swap” samples from RUSH/Broad and 278 are Mayo temporal cortex samples.

# so I want a list of the 10 samples (20 file names?) for RUSH-BROAD sample-swap and the 278 from MayoBrainBank_Dickson, BannerSunHealth_TomBeach

# I'll use this for generating zip file for data processing and for making covariates files.

library(synapseClient) # for synapse upload
library(RCurl) # to grab google doc covariates files
library(gdata) # to read .xlsx file
library(dplyr) # for subsetting data

#Login to Synapse using credentials saved in .synapseConfig file
synapseLogin()

# Grab existing info for inspection and assembly (James organized this across subdirectories in syn2924444)

#RNASeq_fullSampleInformation_NET.xlsx syn3163262
# has RNASubjectID, RNAId, Source, tissue, final diagnosis, RIN, unfiltered age, sex, and Braak
fullSampleInformation <- synGet("syn3163262")
fullSampleInformationFilePath <- getFileLocation(fullSampleInformation)
fullSample <- read.xls(fullSampleInformationFilePath, sheet = 1, header = TRUE, stringsAsFactors = FALSE)

#U01_288_AUT_TCx_RNAseq_Covars-Drives_02-06-2015_1447.xlsx syn3523879
# includes SubjectID, SampleID, seq Run ID, DriveA_FileName1, DriveA_FileName2, DriveB_FileName1, and DriveB_Filename2
rerunSamples  <- synGet("syn3523879")
rerunSamplesFilePath <- getFileLocation(rerunSamples)
rerun <- read.xls(rerunSamplesFilePath, sheet = 1, header = TRUE, stringsAsFactors = FALSE)

# get file names for RUSH-BROAD samples. Some have been rerun. Which ones?

rushBroadSamples <- fullSample[fullSample$Source == "RUSH-BROAD", "RNASubjectId"]
length(rushBroadSamples) #10
length(intersect(rushBroadSamples, rerun$SubjectID)) #10 - so all of them!

rushBroadFile1 <- rerun[rerun$Source == "RUSH-BROAD", "DriveA_FileName1"]
rushBroadFile2 <- rerun[rerun$Source == "RUSH-BROAD", "DriveA_FileName2"]
rushBroadFile3 <- rerun[rerun$Source == "RUSH-BROAD", "DriveB_FileName1"]
rushBroadFile4 <- rerun[rerun$Source == "RUSH-BROAD", "DriveB_FileName2"]

rushBroadFiles <- c(rushBroadFile1[rushBroadFile1 != "NULL"], 
                    rushBroadFile2[rushBroadFile2 != "NULL"], 
                    rushBroadFile3[rushBroadFile3 != "NULL"], 
                    rushBroadFile4[rushBroadFile4 != "NULL"])

# next, build list of files from MayoBrainBank_Dickson and BannerSunHealth_TomBeach, including the most recent if rerun.
MBBSamples <- (rerun$Source == "BannerSunHealth_TomBeach" | rerun$Source == "MayoBrainBank_Dickson")
MBBSampleIds <- rerun[MBBSamples, "SubjectID"]
length(MBBSampleIds) #278. Yay!

MBBFile1 <- rerun[MBBSamples, "DriveA_FileName1"]
MBBFile2 <- rerun[MBBSamples, "DriveA_FileName2"]
MBBFile3 <- rerun[MBBSamples, "DriveB_FileName1"]
MBBFile4 <- rerun[MBBSamples, "DriveB_FileName2"]

MBBFiles <- c(MBBFile1[MBBFile1 != "NULL"],
              MBBFile2[MBBFile2 != "NULL"], 
              MBBFile3[MBBFile3 != "NULL"], 
              MBBFile4[MBBFile4 != "NULL"])

length(MBBFiles) #568... should be 278 * 2 = 556
#it turns out there are 6 files on both drive A and drive B, for the following subject IDs:
# NA08-370
# NA05-054
# NA01-233
# NA03-163
# 00-32
# 10-20
length(unique(MBBFiles)) #556 - that's better. So that's the set of files I want to include in my .tar file from snapr output and samples for the MBB TCX covariates file.
