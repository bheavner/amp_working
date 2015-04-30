# goal: find source files that need to be reprocessed, copy them to different s3 folder for reprocessing.

library(gdata) # to read .xlsx file -  install.packages("gdata")

# samples found to have 0 readcounts when attempting to normalize genecount files; 
# e.g. rownames(expr$samples[expr$samples$lib.size == 0, ])
geneMissedSamples  <- c("11344_TCX", "11480_TCX", "11424_TCX", "11313_TCX", 
                        "11328_TCX", "11476_TCX", "11485_TCX", "11429_TCX", 
                        "11321_TCX", "11291_TCX", "11330_TCX", "11428_TCX", 
                        "11399_TCX", "11264_TCX", "1919_TCX",  "6819_TCX", 
                        "6929_TCX",  "6913_TCX", "6872_TCX", "6880_TCX", 
                        "11412_TCX", "11267_TCX", "11448_TCX", "11279_TCX", 
                        "11358_TCX", "1104_TCX", "11377_TCX", "11287_TCX", 
                        "11275_TCX", "11438_TCX", "132_TCX", "736_TCX", 
                        "896_TCX", "142_TCX", "1098_TCX","731_TCX", 
                        "744_TCX", "1302_TCX", "1951_TCX")

transcriptMissedSamples <- c("00-03", "02-24", "03-41", "04-26", "05-18", 
                               "05-36", "08-23", "10-63", "11-71",  "12-21",  
                               "12-44",  "13-40", "97-45", "99-25", 
                               "NA00-082", "NA00-187", "NA02-081", "NA02-134", 
                               "NA03-138", "NA03-295", "NA04-303", "NA05-054", 
                               "NA05-224", "NA06-082", "NA06-221", "NA07-175", 
                               "NA07-298", "NA08-091", "NA08-166", "NA08-219", 
                               "NA09-245", "NA09-278", "NA10-112", "NA10-117", 
                               "NA10-290", "NA10-381", "NA10-412", "NA11-215", 
                               "NA98-222")

# translate geneMissed Sample IDs to filename IDs
keyPath  <- "/Users/bheavner/Desktop/TCXIDsKey.xlsx"

keyHash <- read.xls(keyPath, sheet = 1, header = TRUE, stringsAsFactors = FALSE)

# replace colnames(geneCounts) with corresponding keyHash$IlluminaSampleID values
geneMissedSamples <- keyHash$Path_ID[match(geneMissedSamples, keyHash$IlluminaSampleID)]

# build list of samples to reprocess
toReRun  <- union(geneMissedSamples, transcriptMissedSamples)

# Get list of files in s3 bucket
system('aws s3 ls s3://mayo-u01-rnaseq/ --recursive | grep -e \".fastq\" | awk \'{print $4}\' > file_list')
# allFiles <- read.table('file_list')

for (sample in toReRun) {
  system(paste('grep ', '"', sample, '\\.', '" ', 'file_list ', '>> sourceFiles', sep = ""))
}

sourceFiles  <- read.table('sourceFiles')

# Since there are 39 samples to rerun, I need 78 paired fastq files. I currently have 110.
# some are missed files, some are duplicates. I'll mannually edit this list.

#find samples with multiple source files
index = 1

for (sample in toReRun) {
  message(
    paste(index, (
      system(
        paste('grep ', '"', sample, '\\.', '" ', 'file_list ', '| wc -l', sep = ""), 
        intern = TRUE)
    )))
    index <- index + 1
}

toReRun[c(5, 9, 16, 17, 20, 27, 28, 29, 32, 35)]

# manually remove duplicate source files from the source list as appropriate for the following samples:
# "05-18"    "11-71"    "NA00-187" "NA02-081" "NA03-295" "NA07-298" "NA08-091" "NA08-166" "NA09-278" "NA10-290"

# looks like NA08-166 was run twice (I bet empty flow cell again). I'm keeping the 2nd run.

sourceFilesEdited  <- read.table('sourceFiles_curated') # it is 78 entries long - good sign.

index = 1
for (sample in toReRun) {
  message(
    paste(index, (
      system(
        paste('grep ', '"', sample, '\\.', '" ', 'sourceFiles_curated ', '| wc -l', sep = ""), 
        intern = TRUE)
    )))
  index <- index + 1
}

# copy source files to TCX_rerun_2/ folder in s3://mayo-u01-rnaseq/ bucket
# test with dryrun
for (file in sourceFilesEdited$V1) {
  message(system(
    paste('aws s3 cp s3://mayo-u01-rnaseq/', file, ' s3://mayo-u01-rnaseq/TCX_rerun_2/', ' --dryrun', sep = "")))
}

# actually do it
for (file in sourceFilesEdited$V1) {
  message(system(
    paste('aws s3 cp s3://mayo-u01-rnaseq/', file, ' s3://mayo-u01-rnaseq/TCX_rerun_2/', sep = "")))
}

# confirm they were all copied
system('aws s3 ls s3://mayo-u01-rnaseq/TCX_rerun_2/ | awk \'{print $4}\' > copied_files')
copied  <- read.table('copied_files') # there are 78! So now run snapr on s3://mayo-u01-rnaseq/TCX_rerun_2/



## some still have zero readcounts (00-03, 11-71, NA00-187, NA05-224, NA09-245, NA10-381). Run them again in TCX_rerun_3/
sourceFilesEdited  <- read.table('secondRunSourceFiles_curated') # it is 12 entries long - good sign.

# copy source files to TCX_rerun_3/ folder in s3://mayo-u01-rnaseq/ bucket
# test with dryrun
for (file in sourceFilesEdited$V1) {
  message(system(
    paste('aws s3 cp s3://mayo-u01-rnaseq/', file, ' s3://mayo-u01-rnaseq/TCX_rerun_3/', ' --dryrun', sep = "")))
}

# actually do it
for (file in sourceFilesEdited$V1) {
  message(system(
    paste('aws s3 cp s3://mayo-u01-rnaseq/', file, ' s3://mayo-u01-rnaseq/TCX_rerun_3/', sep = "")))
}

# confirm they were all copied
system('aws s3 ls s3://mayo-u01-rnaseq/TCX_rerun_3/ | awk \'{print $4}\' > copied_files')
copied  <- read.table('copied_files') # there are 78! So now run snapr on s3://mayo-u01-rnaseq/TCX_rerun_2/