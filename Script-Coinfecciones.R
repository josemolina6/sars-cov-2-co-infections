#!/usr/bin/Rscript

#if (!requireNamespace("BiocManager", quietly = TRUE))
#  install.packages("BiocManager")
#BiocManager::install("dada2", version = "3.13")

library(dada2); packageVersion("dada2")

path <- getwd() # CHANGE ME to the directory containing the fastq files after unzipping.
list.files(path) # verificar que estan los archivos fastq (2 archivos), el script y la Base de datos BD-Sars-cov-2


# Start the clock!
ptm <- proc.time()

# Forward and reverse fastq filenames have format: SAMPLENAME_R1_001.fastq and SAMPLENAME_R2_001.fastq

fnFs <- sort(list.files(path, pattern="_L001_R1_001.fastq.gz", full.names = TRUE))
fnRs <- sort(list.files(path, pattern="_L001_R2_001.fastq.gz", full.names = TRUE))
# Extract sample names, assuming filenames have format: SAMPLENAME_XXX.fastq
sample.names <- sapply(strsplit(basename(fnFs), "_"), `[`, 1)

sample.names

#plotQualityProfile(fnFs[1:2])
#plotQualityProfile(fnRs[1:2])


# Place filtered files in filtered/ subdirectory
filtFs <- file.path(path, "filtered", paste0(sample.names, "_F_filt.fastq.gz"))
filtRs <- file.path(path, "filtered", paste0(sample.names, "_R_filt.fastq.gz"))
names(filtFs) <- sample.names
names(filtRs) <- sample.names

out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs, truncLen=c(50,50),
                     maxN=0, maxEE=c(2,2), truncQ=2, rm.phix=TRUE,
                     compress=TRUE, multithread=FALSE) # On Windows set multithread=FALSE
head(out)

errF <- learnErrors(filtFs, multithread=TRUE)
errR <- learnErrors(filtRs, multithread=TRUE)
plotErrors(errF, nominalQ=TRUE)


dadaFs <- dada(filtFs, err=errF, multithread=TRUE)
dadaRs <- dada(filtRs, err=errR, multithread=TRUE)
dadaFs[[1]]

mergers <- mergePairs(dadaFs, filtFs, dadaRs, filtRs, verbose=TRUE)
# Inspect the merger data.frame from the first sample
head(mergers[[1]])


seqtab <- makeSequenceTable(mergers)
dim(seqtab)


# Inspect distribution of sequence lengths
table(nchar(getSequences(seqtab)))


seqtab.nochim <- removeBimeraDenovo(seqtab, method="consensus", multithread=TRUE, verbose=TRUE)
dim(seqtab.nochim)

sum(seqtab.nochim)/sum(seqtab)


taxa <- assignTaxonomy(seqtab.nochim, "BD_sars-cov-30-nov-2.fasta", multithread=TRUE)
#taxa <- addSpecies(taxa, "silva_species_assignment_v138.1.fa.gz")

taxa.print <- taxa # Removing sequence rownames for display only
rownames(taxa.print) <- NULL

colnames(taxa)<-c("Virus","Linaje")
head(taxa)

tablataxa<-as.data.frame(taxa)
tablataxa[,"Secuencia"]<-(rownames(taxa))


tablafinal<-as.data.frame(t(as.data.frame(seqtab.nochim)))
tablafinal[,"Secuencia"]<-(rownames(tablafinal))

tabla<-merge(tablafinal, tablataxa, by = "Secuencia")
colnames(tabla)[2]<-sample.names[1]
table(tabla$Linaje)

Paexportar<-as.data.frame(table(tabla$Linaje))
colnames(Paexportar)[1]<-sample.names[1]
Paexportar



# Stop the clock
proc.time() - ptm

write.csv(Paexportar, paste("Analisis_de_",sample.names[1],".csv"))



