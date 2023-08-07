Sys.time()
# Create the paths needed to install libraries
.libPaths(c("/dartfs/rc/lab/C/ChaudharyB/shared/R/4.2.3", .libPaths()))
.libPaths(c("/optnfs/el7/Rlibs/4.2",.libPaths()))

library(BiocManager)

# if (!requireNamespace("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# BiocManager::install("dada2", version = "3.16")

library(dada2)
packageVersion("dada2")
library(ShortRead)
packageVersion("ShortRead")
library(Biostrings)
packageVersion("Biostrings")

# path <- "~/Documents/ChaudharyLab/Dust/minimalistSubset" 
path <- "/dartfs/rc/lab/C/ChaudharyB/SequencingRuns/ITS_SSU_5-23-23/ITS"

list.files(path)
fnFs <- sort(list.files(path, pattern = "_R1_001.fastq", full.names = TRUE))
fnRs <- sort(list.files(path, pattern = "_R2_001.fastq", full.names = TRUE))
sample.names <- sapply(strsplit(basename(fnFs), "_"), `[`, 1)

# 
ITS1f <- "CTTGGTCATTTAGAGGAAGTAA"  ## CHANGE ME to your forward primer sequence
ITS2r <- "GCTGCGTTCTTCATCGATGC"  ## CHANGE ME...
# 
# 
allOrients <- function(primer) {
  # Create all orientations of the input sequence
  require(Biostrings)
  dna <- DNAString(primer)  # The Biostrings works w/ DNAString objects rather than character vectors
  orients <- c(Forward = dna, Complement = complement(dna), Reverse = reverse(dna), 
               RevComp = reverseComplement(dna))
  return(sapply(orients, toString))  # Convert back to character vector
}
ITS1f.orients <- allOrients(ITS1f)
ITS2r.orients <- allOrients(ITS2r)
# 
fnFs.filtN <- file.path(path, "filtN", basename(fnFs)) # Put N-filterd files in filtN/ subdirectory
fnRs.filtN <- file.path(path, "filtN", basename(fnRs))

Sys.time()
print("now running filterAndTrim(). This may take a while, go touch some grass.")

filterAndTrim(fnFs, fnFs.filtN, fnRs, fnRs.filtN, maxN = 0, multithread = 6) # Here the multithread is set to six cores but we can set however many we wont. I had an error when I just put TRUE
#
print('primer hits')
primerHits <- function(primer, fn) {
  # Counts number of reads in which the primer is found
  nhits <- vcountPattern(primer, sread(readFastq(fn)), fixed = FALSE)
  return(sum(nhits > 0))
}

print('now running rbind')

rbind(ITS1f.ForwardReads = sapply(ITS1f.orients, primerHits, fn = fnFs.filtN[[1]]), 
      ITS1f.ReverseReads = sapply(ITS1f.orients, primerHits, fn = fnRs.filtN[[1]]), 
      ITS2r.ForwardReads = sapply(ITS2r.orients, primerHits, fn = fnFs.filtN[[1]]), 
      ITS2r.ReverseReads = sapply(ITS2r.orients, primerHits, fn = fnRs.filtN[[1]])
)

Sys.time()
print("installing cutadapt") 

# Install cutadapt via instructions online
cutadapt <- "/dartfs-hpc/rc/home/6/f0066g6/.conda/envs/cutadaptenv/bin/cutadapt" # CHANGE ME to the cutadapt path on your machine
system2(cutadapt, args = "--version") # Run shell commands from R
#
path.cut <- file.path(path, "cutadapt")
if(!dir.exists(path.cut)) dir.create(path.cut)
fnFs.cut <- file.path(path.cut, basename(fnFs))
fnRs.cut <- file.path(path.cut, basename(fnRs))

print("now we are doing some trimming")

ITS1f.RC <- dada2:::rc(ITS1f)
ITS2r.RC <- dada2:::rc(ITS2r)
# Trim FWD and the reverse-complement of REV off of R1 (forward reads)
R1.ITS.flags <- paste("-g", ITS1f, "-a", ITS2r.RC)

# Trim REV and the reverse-complement of FWD off of R2 (reverse reads)
R2.ITS.flags <- paste("-G", ITS2r, "-A", ITS1f.RC)
#
print("now cutadapt is running")
# Run Cutadapt
for(i in seq_along(fnFs)) {
  system2(cutadapt, args = c(R1.ITS.flags, R2.ITS.flags, "-n", 2, # -n 2 required to remove FWD and REV from reads
                             "-o", fnFs.cut[i], "-p", fnRs.cut[i], # output files
                             fnFs.filtN[i], fnRs.filtN[i])) # input files
}
#
print("running rbind() again")
rbind(ITS1f.ForwardReads = sapply(ITS1f.orients, primerHits, fn = fnFs.cut[[1]]),
      ITS1f.ReverseReads = sapply(ITS1f.orients, primerHits, fn = fnRs.cut[[1]]),
      ITS2r.ForwardReads = sapply(ITS2r.orients, primerHits, fn = fnFs.cut[[1]]),
      ITS2r.ReverseReads = sapply(ITS2r.orients, primerHits, fn = fnRs.cut[[1]])
)
#
# # Forward and reverse fastq filenames have the format:
cutFs <- sort(list.files(path.cut, pattern = "_R1_001.fastq.gz", full.names = TRUE))
cutRs <- sort(list.files(path.cut, pattern = "_R2_001.fastq.gz", full.names = TRUE))
#
# Extract sample names, assuming filenames have format:
get.sample.name <- function(fname) strsplit(basename(fname), "_")[[1]][1]
sample.names <- unname(sapply(cutFs, get.sample.name))
head(sample.names)
#
#
#
filtFs <- file.path(path.cut, "filtered", basename(cutFs))
filtRs <- file.path(path.cut, "filtered", basename(cutRs))

out <- filterAndTrim(cutFs, filtFs, cutRs, filtRs, maxN = 0, maxEE = c(2, 2), truncQ = 2,
                     minLen = 50, rm.phix = TRUE, compress = TRUE, multithread = 6)  # on windows, set multithread = FALSE
head(out)

errF <- learnErrors(filtFs, multithread = 40)
errR <- learnErrors(filtRs, multithread = 40)

# plotErrors(errF, nominalQ = TRUE)

dadaFs <- dada(filtFs, err = errF, multithread = 40)
dadaRs <- dada(filtRs, err = errR, multithread = 40)

mergers <- mergePairs(dadaFs, filtFs, dadaRs, filtRs, verbose=TRUE)

seqtab <- makeSequenceTable(mergers)
dim(seqtab)
print('remove chimeras')
seqtab.nochim <- removeBimeraDenovo(seqtab, method="consensus", multithread=TRUE, verbose=TRUE)
write.csv(seqtab, file = 'ASV_Table_ITS.csv')

print('Inspect distribution of sequence lengths:')
table(nchar(getSequences(seqtab.nochim)))

print('Tracking reads through the pipeline')
getN <- function(x) sum(getUniques(x))
track <- cbind(out, sapply(dadaFs, getN), sapply(dadaRs, getN), sapply(mergers, getN),
               rowSums(seqtab.nochim))
# If processing a single sample, remove the sapply calls: e.g. replace
# sapply(dadaFs, getN) with getN(dadaFs)
colnames(track) <- c("input", "filtered", "denoisedF", "denoisedR", "merged", "nonchim")
rownames(track) <- sample.names
head(track)


# Assigning taxonomy with Unite
print('Assigning taxonomy with Unite')

unite.ref <- "/dartfs/rc/lab/C/ChaudharyB/databases/UNITE_general_release_18.07.2023.fasta"  # CHANGE ME to location on your machine
taxa <- assignTaxonomy(seqtab.nochim, unite.ref, multithread = 40, tryRC = TRUE)

taxa.print <- taxa  # Removing sequence rownames for display only
rownames(taxa.print) <- NULL
head(taxa.print)
write.csv(taxa.print, file = "taxaUNITE.csv")
# ### Running things in parralell

Sys.time()
time <- start - Sys.time()
print(time)



# library(plyr)
# # install.packages("doParallel")
# library(doParallel)
# cores <- detectCores()
# cores
# registerDoParallel(cores=cores)
# ?registerDoParallel


