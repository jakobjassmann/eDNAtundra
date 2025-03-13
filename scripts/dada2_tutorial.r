# Metabarcoding Bioinformatics tutorial by Katy Maher and Helen Hipperson 

# Install depdendencies
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install("dada2")

# Dependencies
library(dada2)
library(Biostrings)
library(ShortRead)

# Set working directory
setwd("C:/Users/jakob/repositories/eDNAtundra")

# Get data 
download.file("https://ftp.ebi.ac.uk/pub/training/2024/Exploring-environmental-DNA/12S_data.tar.gz", "12S_data.tar.gz")
untar("12S_data.tar.gz", exdir = "data/")

# Set paths 
input.path <- "data/12S_data/cutadapt"
output.path <- "data/test_metabarcoding"

# Check presence of files (R1 = forward, R2 = reverse reads)
list.files(input.path)

# Assign the two read groups to vectors
fnFs <- sort(list.files(input.path, pattern = "_L001_R1_001.fastq", full.names = TRUE))
fnRs <- sort(list.files(input.path, pattern = "_L001_R2_001.fastq", full.names = TRUE))

# Extract sample names (helper function)
get.sample.name <- function(fname) strsplit(basename(fname), "_")[[1]][2]
sample.names <- unname(sapply(fnFs, get.sample.name))

# check the quality for the first files using dada2
plotQualityProfile(fnFs[1:1])
plotQualityProfile(fnRs[1:1])
# The heatmaps show the frequency of the quality scores at each positoin
# the green line is the median, organge are quartilse
# the read line at the bottom gives the number of reads at a particular lenght
# the x-axis is the lenght (cycle) along the sequence. 

# Check the other files
plotQualityProfile(fnFs[2:3])
plotQualityProfile(fnRs[2:3])

## Automatic filtering of bad reas

# Set out file paths
filtFs <- file.path(output.path, "filtered", basename(fnFs))
filtRs <- file.path(output.path, "filtered", basename(fnRs))

# Automated filtering with filterAndTrim()
# maxN=0 After truncation, sequences with more than 0 Ns will be discarded. (DADA2 requires sequences contain no Ns)
# truncQ = 2 Truncate reads at the first instance of a quality score less than or equal to 2
# rm.phix = TRUE Discard reads that match against the phiX genome
# maxEE=c(2, 2) After truncation, reads with higher than 2 "expected errors" will be discarded
# minLen = 60 Remove reads with length less than 60
# multithread = TRUE input files are filtered in parallel

out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs, maxN = 0, maxEE = c(2, 2),
                     truncQ = 2, minLen = 60, rm.phix = TRUE, compress = TRUE,
                     multithread = TRUE)
out
# The low values are likely the controls. 

## Identification of ASVs with DADA

# 1. Train the error models
errF <- learnErrors(filtFs, multithread = TRUE)
errR <- learnErrors(filtRs, multithread = TRUE)

# Plot errors
plotErrors(errF, nominalQ = TRUE)
# You will see some warning messages, but do not worry about them. This is a message from the plotting function to let you know that there were some zero values in the data plotted (which turn into infinities on the log-scale). This is expected, it results from the fact that not every combination of error type (e.g. A->C) and quality score (e.g. 33) is observed in your data, which is normal.
# To interpret the plots:
# The error rates for each possible transition (e.g. A→C, A→G) are shown
# Red line - expected based on the quality score. (These are plotted when nominalQ = TRUE is included in the plot command)
# Black line - estimate
# Black dots - observed
# What we are expecting to see here is that the observed match up with estimates. Here we can see that the black dots track well with the black line. We can also see that the error rates drop with increasing quality score as we would expect. Sanity checks complete, we can proceed with the analysis.
# If you are worried about what your error plots look like when you fit your own dataset, one possible way to to improve the fit is to try increasing the number of bases the function is using (the default is 100 million).

# 2. Dereplication
exists <- file.exists(filtFs)
# check that all the samples are still present after filtering
derepFs <- derepFastq(filtFs[exists], verbose=TRUE)
derepRs <- derepFastq(filtRs[exists], verbose=TRUE)
# Name the derep-class objects by the sample names
names(derepFs) <- sample.names[exists]
names(derepRs) <- sample.names[exists]

# 3. Inference of ASVs
# Using dada and the error model generated
dadaFs <- dada(derepFs, err = errF, multithread = TRUE)
dadaRs <- dada(derepRs, err = errR, multithread = TRUE)

# 4. Merge paired-end reads
mergers <- mergePairs(dadaFs, derepFs, dadaRs, derepRs, verbose=TRUE)

# 5. Create the ASV matrix
seqtab <- makeSequenceTable(mergers)
dim(seqtab)

6. Chimera detection and removal
seqtab.nochim <- removeBimeraDenovo(seqtab, method="consensus",
                                    multithread=TRUE, verbose=TRUE)
# Number of non-chimeric sequnces
dim(seqtab.nochim)
# Proportion of non-chimeric sequence reads
sum(seqtab.nochim)/sum(seqtab)
# Range of sequence lenghts in ASVs
table(nchar(getSequences(seqtab.nochim)))

# 7. Sequence tracking sanity check
getN <- function(x) sum(getUniques(x))
track <- cbind(out, sapply(dadaFs, getN),
               sapply(dadaRs, getN),
               sapply(mergers, getN),
               rowSums(seqtab.nochim))
colnames(track) <- c("input", "filtered",
                     "denoisedF", "denoisedR",
                     "merged", "nonchim")
rownames(track) <- sample.names
track

## Assign taxonomy
taxa <- assignTaxonomy(seqtab.nochim,
                       "data/12S_data/MiFish_Reference_Database_taxonomy.fasta",
                       multithread=TRUE, verbose = T)

taxa.print <- taxa
rownames(taxa.print) <- NULL
head(taxa.print)

## Futher analysis

# 1. Add meta data
meta<-read.csv("data/12S_data/sample_info.csv",
               row.names = 1)
# Remove controls
seqtab.rmcontrol<-seqtab.nochim[-c(25,31),]
meta.rmcontrol<-meta[-c(25,31),]
# Check number of sequences per sample
rowSums(seqtab.rmcontrol)
# Remove low sequence samples
seqtab.rmlow<-seqtab.rmcontrol[-c(13,14,23),]
meta.rmlow<-meta.rmcontrol[-c(13,14,23),]
# Print the minimum sequence number in one sample.
min(rowSums(seqtab.rmlow))

# 2. Install R packages for further analysis
install.packages(c('vegan', 'ggplot2'))
BiocManager::install("phyloseq")
#Load the packages
library(phyloseq)
library(vegan)
library(ggplot2)

# 3. Rarefraction curves
rarecurve(seqtab.rmlow, step=100, col=meta.rmlow$COLOUR, lwd=2, ylab="ASVs", label=F)
# add a vertical line to represent the fewest sequences in any sample
abline(v=(min(rowSums(seqtab.rmlow))))
# You can see that each sample levels off at a different sequencing depth and each sample has been sequenced at different depths. This also gives us an indication of how many ASVs are in each sample, with the aquarium samples here (dark blue lines) amongst the least diverse, as we'd expect.
# Rarefaction is one of the common approaches people use to correct for differences in library sequencing depth between samples and was originally proposed for traditional ecological studies (Sanders 1968). It involves randomly removing reads until you reach a number often equal to or less than the number of reads in the smallest sample. Random subsampling can result in a loss of data and generate artificial variation. Due to these concerns other methods of transformation have been suggested.
# Other methods of normalisation include normalising counts as a proportion of the total library size (McMurdie and Holmes 2014), centered log-ratio transformations (Gloor et al. 2017), geometric mean pairwise ratios (Chen et al. 2018), and also variance stabilising transformations and relative log expressions (Badri et al. 2018).
# Choosing the type of normalisation method to use for your own dataset is not trivial and the most appropriate method can vary between studies and datasets. Further reading is recommended. A couple of papers to get you started are listed below:
# McMurdie and Holmes (2014)
# Willis (2019)
# Cameron et al. (2021)

# 4. Alpha diversity

# Make phyloseq object
seqtab.rmlow2<-t(as.data.frame(seqtab.rmlow))
phylo_asv <- otu_table(seqtab.rmlow2, taxa_are_rows=TRUE)
phylo_tax <- tax_table(taxa)
phylo_samples <- sample_data(meta.rmlow)

phylo <- phyloseq(phylo_asv, phylo_tax, phylo_samples)

sample_names(phylo)
rank_names(phylo)
sample_variables(phylo)

# Plot alpha diversity
plot_richness(phylo,
              measures=c("Shannon", "Simpson"),
              color = "SITE")
plot_richness(phylo, x="SITE", measures=c("Shannon", "Simpson"),
              color = "SITE") + geom_boxplot()

