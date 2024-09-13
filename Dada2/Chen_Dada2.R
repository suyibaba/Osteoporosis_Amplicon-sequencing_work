######Loading required packages and set the path
library(knitr)
library(rmarkdown)
library(dada2)
packageVersion("dada2")
library(phyloseq)
library(Rcpp)
library(ShortRead)
BiocManager::install("BiocParallel")
library(BiocParallel)
path = "/Users/oluwamayowaakinsuyi/Desktop/archive/all" 
list.files(path)

# Forward and reverse fastq filenames have format: SAMPLENAME = _1.fastq and _2.fastq respectively
# Forward and reverse fastq filenames have format: SAMPLENAME = _1.fastq and _2.fastq respectively
fnFs <- sort(list.files(path, pattern="_1.fastq", full.names = TRUE))
fnFs 

# Extracting sample names, assuming filenames have format: SAMPLENAME_1.fastq
sample.names <- sapply(strsplit(basename(fnFs), "_"), `[`, 1)
sample.names 

## plot quality profiles for both the forward and reverse sequences, this helps to access the quality of the reads 
plotQualityProfile(fnFs[1:6])
plotQualityProfile(fnRs[1:6])

####
filtFs <- file.path(path, "filtered", paste0(sample.names, "_F_filt.fastq.gz"))
filtRs <- file.path(path, "filtered", paste0(sample.names, "_R_filt.fastq.gz"))
names(filtFs) <- sample.names
names(filtRs) <- sample.names
out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs, maxN=0, maxEE=c(2,2), truncLen= c(270,210), truncQ=2, rm.phix=TRUE, trimLeft = c(19,19), compress=TRUE, multithread=TRUE)
head(out)

##Learn The Error Rates
## The plotErrors functions is used for sanity checks, as it shows the error rates for each possible transition from A -C or A-G

errF <- learnErrors(filtFs, multithread=TRUE)
errR <- learnErrors(filtRs, multithread=TRUE)
plotErrors(errF, nominalQ=TRUE)

#Apply the core sample inference algorithm to the filtered and trimmed sequence data.
dadaFs <- dada(filtFs, err=errF, multithread=TRUE)
dadaRs <- dada(filtRs, err=errR, multithread=TRUE)

#Inspecting the returned dada-class object
dadaFs
dadaRs

#######
mergers <- mergePairs(dadaFs, filtFs, dadaRs, filtRs, verbose=TRUE)
View(mergers)
# Inspect the merger data.frame from the first sample
head(mergers[[1]])

#Construction of the Amplicon Sequence Variant Table(ASV)
seqtab4 <- makeSequenceTable(mergers)
View(seqtab4)
#Check the dimension of the seqtab
dim(seqtab4)
#Inspect the distribution of sequence length IN seqtab
table(nchar(getSequences(seqtab4)))
seqtab4


# Remove Chimeras from the seqtab.nochim matrix
seqtab.nochim4 <- removeBimeraDenovo(seqtab4, method="consensus", multithread=TRUE, verbose=TRUE)
class(seqtab.nochim4)

#Check the dimension of the new seqtab.nochim matrix
dim(seqtab.nochim4)
sum(seqtab.nochim4)/sum(seqtab4)
seqtab.nochim4
View(seqtab.nochim4)

# Assign Taxonomy

taxa4_rdp <- assignTaxonomy(seqtab.nochim4,"/Users/if-mcs-1247c/Desktop/line/rdp_train_set_18.fa.gz", multithread=TRUE)
View(taxa4_rdp)
dim(taxa4_rdp)

#Save the DADA2 objects as RDS
saveRDS(seqtab.nochim4, "/Users/if-mcs-1247c/Desktop/untitled\ folder\ 3/Osteoporosis/Chein_2021/seqtab.rds")
saveRDS(taxa4_rdp,"/Users/if-mcs-1247c/Desktop/untitled\ folder\ 3/Osteoporosis/Chein_2021/taxa4_rdp.rds")
