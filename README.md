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

The second number of the results should be a decimal that can be converted to a percent. This tells you that 44% of the reads were quality reads and that the remaining 56% were chimeras and needed to be removed

Next, we will assess the number of sequences as they were throughout every step of the pipeline. This number should decrease through each of the trimming, pairing and filtering steps

```{}
getN <- function(x) sum(getUniques(x))
track <- cbind(out, sapply(dadaFs, getN), sapply(dadaRs, getN), sapply(mergers, getN), rowSums(asvtable.nochim))
```

This step assigns names to all of the columns and prints it so it can be viewed and analyzed:

```{}
colnames(track) <- c("input", "filtered", "denoisedF", "denoisedR", "merged", "nonchim")
rownames(track) <- sample.names
head(track)
```

Next, we need to assign taxonomy to the sequences:

```{}
phyla <- assignTaxonomy(asvtable.nochim, "path/to/your/database", multithread=FALSE)
```

Now we can view the assignments created in the last step:

```{}
phyla.print <- phyla
```

This extracts the row names of the sequences we want to view:

```{}
rownames(phyla.print) <- NULL
head(phyla.print)
```

# Step Two: Saving Data Up to This Point 
STOP HERE! This is an important spot to stop and save your work up until this point to prevent losing any of your files. This is a good stopping point. 

First, save the dataframe that assigns taxonomy to each ASV:
```{}
write.csv(phyla, file=“your/directory/path/phyla_A3.csv”)
```

Next, save the table that displayed the abundance of the chimeras:

```{}
write.csv(asvtable.nochim,file=“your/directory/path/asvtable.nochim_A3.csv”)
```

This saves the tracking table that displays the sequences from each step along the pipeline and how many were lost between each step:

```{}
write.csv(track,file=“your/directory/path/tracking_A3.csv”) 
```

# Step Three: Preparing the Data for Plotting
We will start from the saved CSV files. We are reading all the saved files into the environment as objects so that they can be used to create a plot.
Some of the following steps are for prepping data for your viewing convenience later and are not mandatory. Read through the steps carefully but skip ahead to Step Four if you are looking to skip to plotting with phyloseq 

```{}
phyla <- read.csv(file = 'your/path/to/phyla_A3.csv')
asvtable.nochim <- read.csv(file = 'your/path/to/asvtable.nochim_A3.csv', header = FALSE)
track <- read.csv(file = 'your/path/to/tracking_A3.csv')
```

Here, we will correct the formatting of the objects so that they can be used to make plots through a few steps of flipping and moving or deleting rows and columns. Each step is followed by a view to verify that the step worked which are optional, but recommended here

```{}
flipped_asvtable.nochim<-as.data.frame(t(asvtable.nochim))
View(flipped_asvtable.nochim)
```

In the flipped CSV, column V1 row V1 is blank, so next we're going to copy the first row into the header and then delete the first row in the next two steps
Step 1 is to copy the first row:

```{}
colnames(flipped_asvtable.nochim) <- flipped_asvtable.nochim[1,]
```

This is saying that the column names should be the names in the 1st row in the flipped ASV table and replacing the placeholders with the actual names. Now we can view to verify it worked.

```{}
View(flipped_asvtable.nochim)
```

Step 2 is to now delete the first row and inspect agan to ensure it is correct:

```{}
flipped_asvtable.nochim <- flipped_asvtable.nochim[-1,]
View(flipped_asvtable.nochim)
```

Next, we will change the names of the sequences to 'ASVs' rather than the sequences themselves. We will change all the sequences in the first row to 'ASV1, 2, etc.' so we can always go back and reference which ASV corresponded with which taxonomic name. This is a nice-to-do-step, not a mandatory step but it makes it look nice

```{}
rownames(flipped_asvtable.nochim) <- paste0("ASV", 1:nrow(flipped_asvtable.nochim))
```

Now remove the sequences column and save it:

```{}
flipped_asvtable.nochim_forself <- flipped_asvtable.nochim[,-1]
```

This transposed file can be saved in case it is useful later and now it is formatted nicely 

```{}
write.csv(flipped_asvtable.nochim, file = 'your/path/to/save/flipped_asvtable.nochim.csv')
write.csv(flipped_asvtable.nochim_forself, file ='your/path/to/save/flipped_asvtable.nochim_forself.csv')
```

Now, the flipped data can be saved in one data sheet for yourself using the cbind function. Ensure you have verified your directory before this step

```{}
ASVabund<-cbind(flipped_asvtable.nochim,phyla)
write.csv(ASVabund,file='your/path/to/save/ASVabund.csv')
```

# Step Four: Plotting with Phyloseq and ggplot2
After all the adjustments we made, phyloseq does not like the sequences in column 1 of the final product so we will remove it and view it one more time to verify before beginning the next step

```{}
phyla<-phyla[-1]
View(phyla)
```

STOP HERE! Confirm that you have loaded all your required libraries for the Phyloseq library:

```{}
library(phyloseq)
library(Biostrings)
library(ggplot2)
library(RColorBrewer)
library(tidyverse)
```

We are going to format an object using the phyla data and phyloseq's preferred formatting for analysis. This step will name and create a matrix using the phyla data

```{}
phylamatrix<-as.matrix(phyla)
```

Now, we will make an ASV table in phyloseq formatting using flipped_asvtable.nochim. Here, we are creating and naming the new object while removing the first column that does not fit the phyloseq formatting and viewing it to verify the step 

```{}
otumatrix <-flipped_asvtable.nochim[,-1]
view(otumatrix)
```

Next, we want to convert the object to a matrix format. Since a matrix can only store one type of data, we will have one OTU (ASV) table and one for the phyla:

```{}
otumatrix<-as.matrix(otumatrix)
phylamatrix<-as.matrix(phylamatrix)
```

Now, we need to ensure that the row names are 'ASV' for both files:

```{}
rownames(otumatrix) <- paste0("ASV", 1:nrow(otumatrix))
rownames(phylamatrix) <- paste0("ASV", 1:nrow(otumatrix))
```

Next, we need to ensure that the OTU data is recognized as numeric, not character data:

```{}
class(otumatrix)<-"numeric"
```

Now that we have completed all the formatting, we can use phyloseq to analyze the data. We will use the following commands to tell phyloseq where the OTU and phyla files are, where the data is located, combine the files and print to view:

```{}
OTU=otu_table(otumatrix,taxa_are_rows = TRUE)
PHYLA=tax_table(phylamatrix)
physeq=phyloseq(OTU,PHYLA)
physeq
```

Assign an object with the correct data for ease of use in the next plotting steps:

```{}
sample_names(physeq)
sample_names<-sample_names(physeq)
```

STOP HERE! This portion should be typed in the console if you want the plot to appear in the plots window for proper viewing. First, let's merge the ASVs of each phyla together so easier to interpret the plot. We'll use 'tax_glom' to glom together the taxa in the phyla column

```{}
ps_phylum <- tax_glom(physeq, "Phylum")
```

Next, this glommed data can be used to display the relative abundance of each phylum. Then, we use psmelt to melt away the phyloseq formatting and make it easier for plotting using ggplot2, and we factor the values of Phylum

```{}
ps_phylum_relabun <- transform_sample_counts(ps_phylum, function(ASV) ASV/sum(ASV))
taxa_abundance_table_phylum <- psmelt(ps_phylum_relabun)
taxa_abundance_table_phylum$Phylum<-factor(taxa_abundance_table_phylum$Phylum)
```

Now we can (finally) plot the data. First we create a title so we can use it in our plot line. Then we create a plot as an object (Or you can omit the name <- to just generate the plot and not save it to your environment)

```{}
title = "Relative Abundance of Phyla in Pumice Rock Samples in the South Pacific Ocean"
p_realabun<-plot_bar(ps_phylum_relabun, fill = "Phylum", title=title) + ylab("Relative Abundance (%)")
```
And you're done!

# Additional Steps
OR what if you wanted to use this same data to make a plot of relative abundance for order instead of phyla? We'll follow all the same steps as before but substituting for another type of data within the same data set to create a new plot, starting with formatting:


```{}
ordermatrix<-as.matrix(phyla$Order)
rownames(ordermatrix) <- paste0("ASV", 1:nrow(otumatrix))
ORDER=tax_table(ordermatrix)
ordseq=phyloseq(OTU,ORDER)
sample_names(ordseq)
```

```{}
order_names<-sample_names(ordseq)
ps_order <- tax_glom(physeq, "Order")
ps_order_relabun <- transform_sample_counts(ps_order, function(ASV) ASV/sum(ASV))
taxa_abundance_table_order <- psmelt(ps_order_relabun)
taxa_abundance_table_phylum$Phylum<-factor(taxa_abundance_table_phylum$Phylum)
title2 = "Relative Abundance of Orders in Pumice Rock Samples in South Pacific Ocean"
o_realabun <- plot_bar(ps_order_relabun, fill = "Order", title=title2) + ylab("Relative Abundance (%)")
o_realabun
```

NOW we're done :)
