#BiocManager::install("dada2", version = "3.16")
library(dada2)
library(ggplot2)
library(phyloseq)
library(vegan)
library(tidyverse)

#Path below is for raw files of kept runs only
path <- "data/sequences"
list.files(path)
# Forward and reverse fastq filenames have format: SAMPLENAME_R1_cut.fastq.gz and SAMPLENAME_R2_cut.fastq.gz
# Samplename is everything before the first underscore
fnFs <- sort(list.files(path, pattern="_L001_R1_cut.fastq", full.names = TRUE))
fnRs <- sort(list.files(path, pattern="_L001_R2_cut.fastq", full.names = TRUE))
# Extract sample names, assuming filenames have format: SAMPLENAME_XXX.fastq
sample.names <- sapply(strsplit(basename(fnFs), "_"), `[`, 1) #only needs to be done for Forwards

plotQualityProfile(fnFs[1:8], aggregate = TRUE)
plotQualityProfile(fnRs[9:12], aggregate = TRUE) 

# Perform filtering and trimming
# Assign the filenames for the filtered fastq.gz files.
# Make directory and filenames for the filtered fastqs
filt_path <- file.path(path, "filtered") # Create a "filtered" subdirectory
filtFs <- file.path(filt_path, paste0(sample.names, "_F_filt.fastq.gz")) # Put forwards in filtered directory
filtRs <- file.path(filt_path, paste0(sample.names, "_R_filt.fastq.gz")) # Put reverses in filtered directory

# Filter the forward and reverse reads
out <- filterAndTrim(fwd=fnFs, filt=filtFs, #fastq and corresponding filtered files (for forward reads)
                     rev=fnRs, filt.rev=filtRs, #fastq and corresponding filtered files (for reverse reads)
                     truncLen=c(150,150), #truncate reads after 150 bases (this is based on quality profiles)
                     maxN=0, # 0 is default - sequences with ANY Ns are discarded (DADA does not allow Ns, so must remove sequences with them)
                     maxEE=c(2,2), #After truncation, reads with higher than maxEE "expected errors" will be discarded
                     truncQ=2, # 2 is default - truncate reads at the first instance of a quality score less than or equal to truncQ
                     rm.phix=TRUE, # TRUE is default - discard reads that match against the phiX genome
                     compress=TRUE, # if TRUE, output files are zipped
                     multithread=TRUE) #On Windows set multithread=FALSE

# Learn the Error Rates, it TAKES TIME! do forward first, and then reverse
errF <- learnErrors(filtFs, multithread=TRUE) #calculate forward read error rates
errR <- learnErrors(filtRs, multithread=TRUE) #calculate reverse read error rates

# Dereplicate the filtered fastq files
derepFs <- derepFastq(filtFs, verbose=TRUE)
derepRs <- derepFastq(filtRs, verbose=TRUE)

# Name the derep-class objects by the sample names
names(derepFs) <- sample.names
names(derepRs) <- sample.names

# Infer the sequence variants in each sample
dadaFs <- dada(derepFs, err=errF, multithread=TRUE, pool = TRUE) #75 samples were pooled: 4200598 reads in 371863 unique sequences.
dadaRs <- dada(derepRs, err=errR, multithread=TRUE, pool = TRUE) # 75 samples were pooled: 4200598 reads in 446528 unique sequences.

# Inspecting the dada-class object returned by dada:
dadaFs[[1]]  #dada-class: object describing DADA2 denoising results
#1019 sequence variants were inferred from 11096 input unique sequences.
#Key parameters: OMEGA_A = 1e-40, OMEGA_C = 1e-40, BAND_SIZE = 16

# Merge the denoised forward and reverse reads- by default, will only merge if the forward and reverse reads overlap identically by at least 12 bps
# this can be changed with function arguments
mergers <- mergePairs(dadaFs, derepFs, dadaRs, derepRs, verbose=TRUE)
# Inspect the merger data.frame from the first sample
head(mergers[[2]])

# Construct sequence table (analogous to an OTU table)
seqtab <- makeSequenceTable(mergers) # The sequences being tabled vary in length.

dim(seqtab) #Dimensions of sequence table = 75 (samples) x 9388 (sequence variants)

# Inspect distribution of sequence lengths
table(nchar(getSequences(seqtab))) #anywhere from 156 to 288 bps, but VAST majority @ 250-258 bps
#seqtab2 <- seqtab[,nchar(colnames(seqtab)) %in% 249:261] #Run this line remove non-target-length sequences from your sequence table


#Remove chimeric sequences:
seqtab.nochim <- removeBimeraDenovo(seqtab,
                                    method="consensus", #Samples are independently checked for bimeras, and a consensus decision on each sequence variant is made (as opposed to pooled)
                                    multithread=TRUE, verbose=TRUE)
dim(seqtab.nochim) #dimensions of sequence table with chimeras removed 75 x 5461
1-sum(seqtab.nochim)/sum(seqtab) #0.03 - this tells us what percentage of seqtab was chimeric

# Track reads through the pipeline
# As a final check of our progress, we’ll look at the number of reads that made it through each step in the pipeline
getN <- function(x) sum(getUniques(x))
track <- cbind(out, #input and filtered
               sapply(dadaFs, getN), #denoised
               sapply(mergers, getN), #merged
               rowSums(seqtab), #tabled
               rowSums(seqtab.nochim)) #nonchim
# Note: If processing a single sample, remove the sapply calls: e.g. replace sapply(dadaFs, getN) with getN(dadaFs)
colnames(track) <- c("input", "filtered", "denoised", "merged", "tabled", "nonchim") #name columns appropriately
rownames(track) <- sample.names #label rows with Sample IDs
head(track) #inspect
write.table(track, "data/readStats.txt",sep="\t",col.names=NA) #uncomment this to write to .txt file

# SAVE THIS FILE SO YOU DON'T HAVE TO REPEAT ALL OF THE ABOVE STEPS, adjust name 
#uncomment the below to write to .txt file
saveRDS(seqtab.nochim, file="data/seqtab.nochim.rds")

# Assign taxonomy
# Make sure the appropriate database is available in the DADA2 directory
taxa <- assignTaxonomy(seqtab.nochim, "data/silva_nr_v132_train_set.fa.gz", multithread=TRUE) #Classifies sequences against reference training dataset (silva)

taxon <- as.data.frame(taxa,stringsAsFactors=FALSE) #convert to dataframe 
#"taxon" has each unique ASV in a row, with columns for Kingdom thru Genus

# FIX the NAs in the taxa table
taxon$Phylum[is.na(taxon$Phylum)] <- taxon$Kingdom[is.na(taxon$Phylum)]
taxon$Class[is.na(taxon$Class)] <- taxon$Phylum[is.na(taxon$Class)]
taxon$Order[is.na(taxon$Order)] <- taxon$Class[is.na(taxon$Order)]
taxon$Family[is.na(taxon$Family)] <- taxon$Order[is.na(taxon$Family)]
taxon$Genus[is.na(taxon$Genus)] <- taxon$Family[is.na(taxon$Genus)]

write.table(taxon,"data/dada2_silva_taxa_table.txt",sep="\t",col.names=NA)
write.table(seqtab.nochim, "data/dada2_silva_otu_table.txt",sep="\t",col.names=NA)

# check total number of reads and number that made it through DADA2
track%>%
  as.data.frame()%>%
  rownames_to_column("ID")%>%
  filter(!ID == "BLANKRun1")%>%
  summarize(totalReads = sum(input),
            totalNonchim = sum(nonchim),
            mean = mean(nonchim),
         min = min(nonchim), 
         max = max(nonchim))
