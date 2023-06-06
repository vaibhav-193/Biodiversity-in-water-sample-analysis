path <- "D:/BI Project"

fnFs <- sort(list.files(path,pattern = "_1_trim.fastq"))
fnRs <- sort(list.files(path,pattern = "_2_trim.fastq"))


sample.names <- sapply(strsplit(fnFs, "_"),'[',1)

fnFs <- file.path(path, fnFs)
fnRs <- file.path(path, fnRs)

plotQualityProfile(fnFs)
plotQualityProfile(fnRs)

filt_path <- file.path(path, "filtered")

filtFs <- file.path(filt_path, paste0(sample.names, "_F_filt.fastq"))
filtRs <- file.path(filt_path, paste0(sample.names, "_R_filt.fastq"))

out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs, truncLen = c(100, 100),
                     maxN=0, maxEE=c(1,1), truncQ=2, rm.phix=TRUE,
                     compress=TRUE, multithread = FALSE)
head(out)
#Estimate the error model for DADA2 algorithm using reverse reads
errF <- learnErrors(filtFs)

errR <- learnErrors(filtRs)


plotErrors(errF, nominalQ = TRUE)
plotErrors(errR, nominalQ = TRUE)

dadaFs <- dada(filtFs, err=errF, multithread = FALSE)
dadaRS <- dada(filtRs, err=errR)

mergers <- mergePairs(dadaFs, filtFs, dadaRS, filtRs, verbose=TRUE)
# Inspect the merger data.frame from the first sample
head(mergers[[1]])

seqtab <- makeSequenceTable(mergers)
dim(seqtab)

table(nchar(getSequences(seqtab)))

seqtab.nochim <- removeBimeraDenovo(seqtab, method="consensus", multithread=FALSE, verbose=TRUE)
dim(seqtab.nochim)

sum(seqtab.nochim)/sum(seqtab)

getN <- function(x) sum(getUniques(x))
track <- cbind(out, sapply(dadaFs, getN), sapply(dadaRS, getN), sapply(mergers, getN), rowSums(seqtab.nochim))
# If processing a single sample, remove the sapply calls: e.g. replace sapply(dadaFs, getN) with getN(dadaFs)
colnames(track) <- c("input", "filtered", "denoisedF", "denoisedR", "merged", "nonchim")
rownames(track) <- sample.names
head(track)

taxa <- assignTaxonomy(seqtab.nochim, file.path(path, "silva_nr99_v138.1_train_set.fa.gz"), multithread=FALSE)

taxa.print <- taxa # Removing sequence rownames for display only
rownames(taxa.print) <- NULL
head(taxa.print)

samples.out <- rownames(seqtab.nochim)

samdf <- data.frame(Country= c("Germany", "India"), Region=c('Arctic', 'Equator'))
rownames(samdf) <- samples.out

#Construct phyloseq object directly from dada2 
ps <- phyloseq(otu_table(seqtab.nochim, taxa_are_rows=FALSE), 
               sample_data(samdf), 
               tax_table(taxa))
ps <- prune_samples(sample_names(ps) != "Mock", ps) # Remove mock sample



dna <- Biostrings::DNAStringSet(taxa_names(ps))
names(dna) <- taxa_names(ps)
ps <- merge_phyloseq(ps, dna)
taxa_names(ps) <- paste0("ASV", seq(ntaxa(ps)))
ps

#visualize alpha_diversity
plot_richness(ps, x="Country", measures=c("Shannon", "Simpson"), color="Region")

#ordinate
# Transform data to proportions as appropriate for Bray-Curtis distances
ps.prop <- transform_sample_counts(ps, function(otu) otu/sum(otu))
ord.nmds.bray <- ordinate(ps.prop, method="PCoA", distance="bray")
ord.nmds.bray <- ordinate(ps.prop, method="DCA", distance="euclidean")

top20 <- names(sort(taxa_sums(ps), decreasing=TRUE))[1:20]
ps.top20 <- transform_sample_counts(ps, function(OTU) OTU/sum(OTU))
ps.top20 <- prune_taxa(top20, ps.top20)
plot_bar(ps.top20, x="Country", fill="Family") + facet_wrap(~Region, scales="free_x")
