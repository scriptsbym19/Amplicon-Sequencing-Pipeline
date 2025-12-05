# Amplicon-Sequencing-Pipeline
Amplicon sequencing of the 16S rRNA gene using the DADA2 pipeline and R Studio 

This is an R markdown file for processing Illumina sequencing output fast.q files using the DADA2 pipeline
# Set Up 
You should always start by setting your working directory and loading your libraries 

```{}
setwd("your/working/directory/path")
```

These are the libraries needed for this pipeline:

```{}
library(tidyverse)
library(ggplot2)
library(scales)
library(RColorBrewer)
library(phyloseq)
library(Biostrings)
library(dada2)
```
Your raw sequencing data should be saved in your working directory or wherever you direct your path so that the files can be called on in the next steps. 
By creating an object (path) with the directory, the files can easily be called on later:

```{}
path <- "your/path/to/your/sequencing/files"
```

This step is optional but it verifies that your data was properly unzipped and uploaded:

```{}
list.files(path)
```

# Step One: Organizing and Preparing Your Data 
First, we need to organize the data. You will need to know what naming system is being used depending on where or how it is being done. In this case, R1 is for the forward reads and R2 is for the reverse reads
Start by differentiating between the forward and reverse reads. We will use fnFs for the forward reads:

```{}
fnFs <- sort(list.files("C:/Users/morga/Downloads/WorkingD_R/A3_ASV-2025", pattern="_R1_001.fastq", full.names = TRUE))
```

We will use fnRs for the reverse reads:

```{}
fnRs <- sort(list.files("C:/Users/morga/Downloads/WorkingD_R/A3_ASV-2025", pattern="_R2_001.fastq", full.names = TRUE))
```

Next, we extract the file names from the data:

```{}
sample.names<-sapply(strsplit(basename(fnFs),"_"),`[`,1)
```

After the files are organized and named, we need to assess the quality of the reads. This step will generate quality plots that we can use to determine where to truncate (trim) the reads to get the highest quality final product 
This visualizes the forward reads:

```{}
plotQualityProfile(fnFs[1:2])
```

This visualizes the reverse reads: 

```{}
plotQualityProfile(fnRs[1:2])
```

STOP HERE! This is where you need to assess your qaulity profile to determine where to truncate the forward and reverse reads. It is usually a good place to truncate when the quality score dips past 30

This step is creating/naming new folders for the filtered forward and reverse files:

```{}
filtFs<-file.path("your/path/to/your/sequencing/files","filtered",paste0(sample.names,"_F_filt.fastq.gz"))
filtRs<-file.path("your/path/to/your/sequencing/files","filtered",paste0(sample.names,"_R_filt.fastq.gz"))
```

The object sample.names is now the filtered files:

```{}
names(filtFs)<-sample.names
names(filtRs)<-sample.names
```

This is where we truncate, this number varies every time so always review the profile and change the ranges in the 'truncLen' part of the script
multithread will = FALSE on Windows but TRUE on Mac

```{}
out<-filterAndTrim(fnFs,filtFs,fnRs,filtRs,truncLen=c(240,190),maxN=0,maxEE=c(2,2),truncQ=2,rm.phix=TRUE,compress=TRUE,multithread=FALSE)
head(out)
```

This step can take a while to run. We will check the error rates based off the trimmed data:

```{}
errF <- learnErrors(filtFs, multithread=FALSE)
errR <- learnErrors(filtRs, multithread=FALSE)
```

This plot is an optional verification step to ensure the error rates are normal:

```{}
plotErrors(errF, nominalQ=TRUE)
```

The next step is to apply the core sample inference algorithm to the filtered data to determine the unique sequences before we try to merge paired reads
The same as the other steps, start with the forward reads and repeat for the reverse:

```{}
dadaFs<-dada(filtFs,err=errF,multithread = FALSE)
dadaRs <- dada(filtRs, err=errR, multithread=FALSE)
```

Then, as an optional step, you can stop again to inspect the data so far:

```{}
dadaFs[[1]]
```

# Step Two: Merging Paired Reads, Removing Chimeras and Assigning Taxonomy
This step will bring the forward and reverse reads together to determine the number of paired unique sequences. The number of reads out determined previously will be smaller since not all the reads will be able to be paired:

```{}
mergers <- mergePairs(dadaFs, filtFs, dadaRs, filtRs, verbose=TRUE)
```

After creating a new object (mergers), you inspect it from the first sample. If 'head' isn't working, you can copy and paste it outside of the chunk or into the console directly to try running it:

```{}
head(mergers[[1]])
```

This step can take a while to run. Next, we assemble the ASV (amplicon sequence variant) table:

```{}
asvtable <- makeSequenceTable(mergers) 
dim(asvtable)
```

This step is to inspect the distribution of the sequence lengths:

```{}
table(nchar(getSequences(asvtable)))
```

Next, we need to remove chimeras 

```{}
asvtable.nochim <- removeBimeraDenovo(asvtable, method="consensus", multithread=FALSE, verbose=TRUE)
```

```{}
dim(asvtable.nochim)
sum(asvtable.nochim)/sum(asvtable)
```





