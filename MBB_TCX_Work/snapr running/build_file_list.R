snapr_files <- read.table("snapr_files", stringsAsFactors = FALSE)
source_files  <- read.table("source_files", stringsAsFactors = FALSE)

# need to learn about data frames
snapr_files <- snapr_files$V1
source_files <- source_files$V1


# get sample ID from snapr files:
# strip snapr/ from beginning
processed_samples <- gsub("snapr/", "", snapr_files)

#strip .snap.bam from end
processed_samples <- gsub(".snap.bam", "", processed_samples)

# strip everything after a period
processed_samples <- gsub("\\..+", "", processed_samples)


# get sample ID from source files:
# strip everything before last /
source_samples <- gsub(".+/", "", source_files)

#strip .fastq.gz from the end
source_samples <- gsub(".fastq.gz", "", source_samples)

#strip ._ from beginning of some sample names
source_samples <- gsub("\\._", "", source_samples)

# strip everything after a period
source_samples <- gsub("\\..+", "", source_samples)

# now, find files in source_files that haven't been processed
source_samples %in% processed_samples

# and grab list of files not processed
not_run <- source_files[!(source_samples %in% processed_samples)]

# don't do the following - it breaks downstream.
#not_run <- paste("s3://mayo-u01-rnaseq", not_run)
#not_run <- gsub(" ", "/", not_run)

write(not_run, "not_run")
