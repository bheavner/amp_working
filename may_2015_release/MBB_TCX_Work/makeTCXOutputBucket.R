# goal: use 01_288_AUT_TCx_RNAseq_Covars-Drives_02-06-2015_1447.xlsx syn3523879 to generate list of files that should be included for readcounts and data release; copy SNAPR output files from s3://mayo-u01-rnaseq/snapr to s3://mayo-u01-rnaseq/TCXOutput for subsequent processing steps.

library(synapseClient) # for synapse upload
library(gdata) # to read .xlsx file -  install.packages("gdata")

#Login to Synapse using credentials saved in .synapseConfig file
synapseLogin()

#U01_288_AUT_TCx_RNAseq_Covars-Drives_02-06-2015_1447.xlsx syn3523879
# includes SubjectID, SampleID, seq Run ID, DriveA_FileName1, DriveA_FileName2, DriveB_FileName1, and DriveB_Filename2
rerunSamples  <- synGet("syn3523879")
rerunSamplesFilePath <- getFileLocation(rerunSamples)
rerun <- read.xls(rerunSamplesFilePath, sheet = 1, header = TRUE, stringsAsFactors = FALSE)

# I only want the MBB TCX samples
MBBSamples <- (rerun$Source == "BannerSunHealth_TomBeach" | rerun$Source == "MayoBrainBank_Dickson")

# generate list of unique strings that should be in the SNAPR output from fastq file names (sample IDs not sufficient b/c some IDs rerun)
MBBFile1 <- rerun[MBBSamples, "DriveA_FileName1"]
MBBFile2 <- rerun[MBBSamples, "DriveA_FileName2"]
MBBFile3 <- rerun[MBBSamples, "DriveB_FileName1"]
MBBFile4 <- rerun[MBBSamples, "DriveB_FileName2"]

MBBFiles <- c(MBBFile1[MBBFile1 != "NULL"],
              MBBFile2[MBBFile2 != "NULL"], 
              MBBFile3[MBBFile3 != "NULL"], 
              MBBFile4[MBBFile4 != "NULL"])

MBBFiles  <- unique(MBBFiles) # there are 6 files on both drive A and drive B

# strip off path info
MBBString <- sub(".+/", '', MBBFiles)

# strip off things after L#
MBBString <- sub("_L[0-9]_.+", '', MBBString)

# I expect 278 are Mayo temporal cortex samples
MBBString <- unique(MBBString)
length(MBBString) # 278 - yay!


# list of files currently in s3://mayo-u01-rnaseq/snapr
system("aws s3 ls s3://mayo-u01-rnaseq/snapr/ > all_snapr_out.txt")
inBucket <- read.table("all_snapr_out.txt", stringsAsFactors = FALSE, skip=1)
inBucket <- inBucket$V4

# we don't care about .bam or .bam.bai files at the moment
notBam <- unique(inBucket[!grepl(".bam", inBucket)])

# if all are present, should be length(MBBString) * 12 files.
length(notBam)/12 #252.5 - so some are missing.


# DIVERSION: what is missing from s3://mayo-u01-rnaseq/snapr that is in the MBBString list?
processedList <- unique(sub(".snap.+", '', notBam))
processedList <- unique(sub("_L[0-9].+", '', processedList)) # N = 303

notProcessed  <- setdiff(MBBString, processedList) # N = 7
notNeeded  <- setdiff(processedList, MBBString)

# I need to run snapr on the NotProcessed List. Only 7 samples, though!
# Need to make file list for these samples:
toReRun <- c("NA04-258", "NA05-327", "05-18", "06-05", "06-15", "09-34", "09-50")

#source_files  <- read.table("snapr\ running/source_files", stringsAsFactors = FALSE)

# need to learn about data frames
#source_files <- source_files$V1

system("aws s3 ls s3://mayo-u01-rnaseq/ --recursive > inBucket.txt")
source_files <- read.table("inBucket.txt", stringsAsFactors = FALSE)
source_files <- source_files$V4

# get sample ID from source files:
# strip everything before last /
source_samples <- gsub(".+/", "", source_files)

#strip .fastq.gz from the end
source_samples <- gsub(".fastq.gz", "", source_samples)

#strip ._ from beginning of some sample names
source_samples <- gsub("\\._", "", source_samples)

# strip everything after a period
source_samples <- gsub("\\..+", "", source_samples)

# now, find samples in source_files that haven't been processed
toReRun %in% source_samples

# which files _should_ be there? (look at .xls)

# missing files:
# sample NA04-258 (from Drive A)
# 141212_SN414_0414_BC5RNRACXX/primary/NA04-258.FCC5RNRACXX_L5_R1_ITTAGGC.fastq.gz
# 141212_SN414_0414_BC5RNRACXX/primary/NA04-258.FCC5RNRACXX_L5_R2_ITTAGGC.fastq.gz

# sample NA05-327 (from Drive A; other is in s3://mayo-u01-rnaseq/tcx-rnaseq )
# 141212_SN414_0414_BC5RNRACXX/primary/NA05-327.FCC5RNRACXX_L1_R1_IACTTGA.fastq.gz

# sample 05-18 (from Drive B)
# 150102_SN616_0383_AC6366ACXX/primary/05-18.FCC6366ACXX_L6_R1_IATCACG.fastq.gz
# 150102_SN616_0383_AC6366ACXX/primary/05-18.FCC6366ACXX_L6_R2_IATCACG.fastq.gz

# sample 06-05 (from Drive B)
# 150102_SN616_0383_AC6366ACXX/primary/06-05.FCC6366ACXX_L4_R1_IGTGGCC.fastq.gz
# 150102_SN616_0383_AC6366ACXX/primary/06-05.FCC6366ACXX_L4_R2_IGTGGCC.fastq.gz

# sample 06-15 (from Drive B)
# 150102_SN616_0383_AC6366ACXX/primary/06-15.FCC6366ACXX_L2_R1_IATCACG.fastq.gz
# 150102_SN616_0383_AC6366ACXX/primary/06-15.FCC6366ACXX_L2_R2_IATCACG.fastq.gz

# sample 09-34 (from Drive B)
# 150102_SN616_0383_AC6366ACXX/primary/09-34.FCC6366ACXX_L3_R1_IGTGGCC.fastq.gz
# 150102_SN616_0383_AC6366ACXX/primary/09-34.FCC6366ACXX_L3_R2_IGTGGCC.fastq.gz

# sample 09-50 (from Drive B)
# 150102_SN616_0383_AC6366ACXX/primary/09-50.FCC6366ACXX_L2_R1_IGAGTGG.fastq.gz
# 150102_SN616_0383_AC6366ACXX/primary/09-50.FCC6366ACXX_L2_R2_IGAGTGG.fastq.gz

# which _are_ in s3?
# grab list of files not processed to try to process them again
source_files[!(source_samples %in% processed_samples)]

# at the moment, I can't find the source fastq files in our S3 buckets

# finally, grab list of files not processed to try to process them again
not_run <- source_files[!(source_samples %in% processed_samples)]

source_files[source_samples %in% toReRun]

## BACK TO WORK - COPY FILES WE HAVE AND WANT from s3://mayo-u01-rnaseq/snapr to s3://mayo-u01-rnaseq/TCXOutput

# select filenames from notBam that are on the MBBString list, and copy them to new bucket

toCopy <- unique(grep(paste(MBBString, collapse="|"), notBam, value=TRUE))

for (index in 1:length(toCopy)) {
  system(paste("aws s3 cp s3://mayo-u01-rnaseq/snapr/", 
               toCopy[index], 
               " s3://mayo-u01-rnaseq/TCXOutput/", 
               #" --dryrun", # comment out to run!
               sep=""))
}

# confirm copy worked;

# if worked, notBam should have the same things as s3://mayo-u01-rnaseq/snapr/ (more complicated in check_snapr.R)

# get dir
system(paste("aws s3 ls s3://mayo-u01-rnaseq/TCXOutput/", "> TCXCopied.txt"))
copiedFiles <- read.table("TCXCopied.txt", stringsAsFactors = FALSE, skip=1)
copiedFiles <- copiedFiles$V4
setdiff(copiedFiles, notBam) #05689621.FCC5RK6ACXX_L4IACTGAT.snap.fusions.reads.fa not copied. Do it now.

# yay!