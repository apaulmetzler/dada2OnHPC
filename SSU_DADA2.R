start <- Sys.time()
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
path <- "/dartfs/rc/lab/C/ChaudharyB/SequencingRuns/ITS_SSU_5-23-23/SSU"

list.files(path)
fnFs <- sort(list.files(path, pattern = "_R1_001.fastq", full.names = TRUE))
fnRs <- sort(list.files(path, pattern = "_R2_001.fastq", full.names = TRUE))
sample.names <- sapply(strsplit(basename(fnFs), "_"), `[`, 1)

# 
WANDAf <- "CAGCCGCGGTAATTCCAGCT"  ## CHANGE ME to your forward primer sequence
AML2r <- "GAACCCAAACACTTTGGTTTCC"  ## CHANGE ME...
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
WANDAf.orients <- allOrients(WANDAf)
AML2r.orients <- allOrients(AML2r)
# 
fnFs.filtN <- file.path(path, "filtN", basename(fnFs)) # Put N-filterd files in filtN/ subdirectory
fnRs.filtN <- file.path(path, "filtN", basename(fnRs))

Sys.time()
print("now running filterAndTrim(). This may take a while, go touch some grass.")

filterAndTrim(fnFs, fnFs.filtN, fnRs, fnRs.filtN, maxN = 0, multithread = 12) # Here the multithread is set to six cores but we can set however many we wont. I had an error when I just put TRUE
#
print('primer hits')
primerHits <- function(primer, fn) {
  # Counts number of reads in which the primer is found
  nhits <- vcountPattern(primer, sread(readFastq(fn)), fixed = FALSE)
  return(sum(nhits > 0))
}

print('now running rbind')

rbind(WANDAf.ForwardReads = sapply(WANDAf.orients, primerHits, fn = fnFs.filtN[[1]]), 
      WANDAf.ReverseReads = sapply(WANDAf.orients, primerHits, fn = fnRs.filtN[[1]]), 
      AML2r.ForwardReads = sapply(AML2r.orients, primerHits, fn = fnFs.filtN[[1]]), 
      AML2r.ReverseReads = sapply(AML2r.orients, primerHits, fn = fnRs.filtN[[1]])
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

WANDAf.RC <- dada2:::rc(WANDAf)
AML2r.RC <- dada2:::rc(AML2r)
# Trim FWD and the reverse-complement of REV off of R1 (forward reads)
R1.SSU.flags <- paste("-g", WANDAf, "-a", AML2r.RC)

# Trim REV and the reverse-complement of FWD off of R2 (reverse reads)
R2.SSU.flags <- paste("-G", AML2r, "-A", WANDAf.RC)
#
print("now cutadapt is running")
# Run Cutadapt
for(i in seq_along(fnFs)) {
  system2(cutadapt, args = c(R1.SSU.flags, R2.SSU.flags, "-n", 2, # -n 2 required to remove FWD and REV from reads
                             "-o", fnFs.cut[i], "-p", fnRs.cut[i], # output files
                             fnFs.filtN[i], fnRs.filtN[i])) # input files
}
#
print("running rbind() again")
rbind(WANDAf.ForwardReads = sapply(WANDAf.orients, primerHits, fn = fnFs.cut[[1]]),
      WANDAf.ReverseReads = sapply(WANDAf.orients, primerHits, fn = fnRs.cut[[1]]),
      AML2r.ForwardReads = sapply(AML2r.orients, primerHits, fn = fnFs.cut[[1]]),
      AML2r.ReverseReads = sapply(AML2r.orients, primerHits, fn = fnRs.cut[[1]])
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
print('running filter and trim again')
out <- filterAndTrim(cutFs, filtFs, cutRs, filtRs, maxN = 0, maxEE = c(3, 3), truncQ = 3,
                     minLen = 50, rm.phix = TRUE, compress = TRUE, multithread = 12)  # on windows, set multithread = FALSE
head(out)
print('finished running filter and trim')
errF <- learnErrors(filtFs, multithread = 50)
errR <- learnErrors(filtRs, multithread = 50)
print('error objects created')
# plotErrors(errF, nominalQ = TRUE)

dadaFs <- dada(filtFs, err = errF, multithread = 50)
dadaRs <- dada(filtRs, err = errR, multithread = 50)
print('dada2 objects created')
print('running mergePairs()')
mergers <- mergePairs(dadaFs, filtFs, dadaRs, filtRs, verbose=TRUE)
print('running makeSequenceTable')
seqtab <- makeSequenceTable(mergers)
dim(seqtab)
print('remove chimeras')
seqtab.nochim <- removeBimeraDenovo(seqtab, method="consensus", multithread=TRUE, verbose=TRUE)
write.csv(seqtab, file = 'ASV_Table_SSU.csv')

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

# 
# Assigning taxonomy with Unite
print('Assigning taxonomy with Maarjam')

maarjam.ref <- "/dartfs/rc/lab/C/ChaudharyB/databases/MaarjAMDADA2FormattedDatabase.fasta"  # CHANGE ME to location on your machine
taxa <- assignTaxonomy(seqtab.nochim, maarjam.ref, multithread = 40, tryRC = TRUE)

taxa.print <- taxa  # Removing sequence rownames for display only
rownames(taxa.print) <- NULL
head(taxa.print)
write.csv(taxa.print, file = "taxaMAARJAM.csv")
#
Sys.time()
end <- Sys.time()
time <- start-end
print(time)