# goal: given list of sample IDs, check for snapr output files. In this case, I'll use the 10 samples from the Rush/Broad sample swap.

sampleSwapIds <- c("10249336", "50100518", "50104008", "20177982", "36492755",
                   "05689621", "11444465", "10315029", "20270920", "11615242")

for (index in 1:length(sampleSwapIds) ) {
  system(paste("aws s3 ls s3://mayo-u01-rnaseq/snapr/|grep", 
               sampleSwapIds[index], 
         ">> processOutput.txt"))
}

processedFiles <- read.table("processOutput.txt", stringsAsFactors = FALSE)

numMissed  <- length(sampleSwapIds) - (nrow(processedFiles) / 12)
numMissed

# which is missing?
outputs  <- rep(NA,length(sampleSwapIds))

if(numMissed > 0) {
  for (index in 1:length(sampleSwapIds) ) {
    outputs[index]  <- sum(grepl(sampleSwapIds[index], processedFiles$V4))
  }
}

sampleSwapIds[!as.logical(outputs)]
#11615242 is missing.

# copy sample swap output (not the .bam or .bai files for now) to s3://mayo-u01-rnaseq/rush_broad_ss/

toCopy <- processedFiles$V4[!grepl(".bam", processedFiles$V4)]
toCopy <- unique(toCopy)

for (index in 1:length(toCopy)) {
  system(paste("aws s3 cp s3://mayo-u01-rnaseq/snapr/", 
               toCopy[index], 
               " s3://mayo-u01-rnaseq/rush_broad_ss/snapr/", 
               " --dryrun", # comment out to run!
               sep=""))
}

# confirm copy worked; if so, remove output from s3://mayo-u01-rnaseq/snapr/

# if worked, processedFiles$V4[!grepl(".bam", processedFiles$V4)] should have the same things as s3://mayo-u01-rnaseq/rush_broad_ss/snapr/

# get dir
system(paste("aws s3 ls s3://mayo-u01-rnaseq/rush_broad_ss/snapr/", "> copied.txt"))
copiedFiles <- read.table("copied.txt", stringsAsFactors = FALSE, skip=1)

nonBam <- processedFiles$V4[!grepl(".bam", processedFiles$V4)]

setdiff(nonBam, copiedFiles$V4) #05689621.FCC5RK6ACXX_L4IACTGAT.snap.fusions.reads.fa not copied. Do it now.


# remove output from snapr folder
#CAREFUL HERE! CHECK YOU'VE COPIED WHAT YOU WANT FIRST!!!
for (index in 1:nrow(copiedFiles)) {
  system(paste("aws s3 rm s3://mayo-u01-rnaseq/snapr/", copiedFiles$V4[index], 
               #" --dryrun", # comment out to run -- CAREFUL!
               sep=""))
}
