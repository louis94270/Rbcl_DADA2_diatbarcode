##### Keck's rbcL DADA2 pipeline
##### https://github.com/fkeck/DADA2_diatoms_pipeline
##### modified by L. Astorg on Apr 22 2022

library(dada2)
library(tidyverse)

# Your path to your sequences here
path <- "/path/to/demultiplexed/sequences"

path_results <- file.path(path, "results")
if (!dir.exists(path_results)) dir.create(path_results)

# Set patterns to discriminate your forward and reverse read files
file_pattern <- c("F" = "_R1_001.fastq", "R" = "_R2_001.fastq")

#### REMOVAL OF PRIMERS ####
# This part is optional. If the primers have already been removed, skip it.
# This script skips all the verifications and checkup of the original tutorial
# If you are not sure about what you are doing, it is strongly advised to read
# the official DADA2 ITS Pipeline Workflow.
# https://benjjneb.github.io/dada2/ITS_workflow.html

fas_Fs_raw <- sort(list.files(path,
                              pattern = file_pattern["F"],
                              full.names = TRUE))
fas_Rs_raw <- sort(list.files(path,
                              pattern = file_pattern["R"],
                              full.names = TRUE))

# This is our set of primers (from Vasselon et al. 2017 )
FWD <- c("AGGTGAAGTAAAAGGTTCWTACTTAAA",
         "AGGTGAAGTTAAAGGTTCWTAYTTAAA",
         "AGGTGAAACTAAAGGTTCWTACTTAAA")

REV <- c("CCTTCTAATTTACCWACWACTG",
         "CCTTCTAATTTACCWACAACAG")

FWD_RC <- dada2:::rc(FWD)
REV_RC <- dada2:::rc(REV)

path_cut <- file.path(path, "cutadapt")
if (!dir.exists(path_cut)) dir.create(path_cut)
fas_Fs_cut <- file.path(path_cut, basename(fas_Fs_raw))
fas_Rs_cut <- file.path(path_cut, basename(fas_Rs_raw))

R1_flags <- paste(paste("-g", FWD, collapse = " "),
                  paste("-a", REV_RC, collapse = " "))
R2_flags <- paste(paste("-G", REV, collapse = " "),
                  paste("-A", FWD_RC, collapse = " "))

cutadapt <- "/path/to/cutadapt" # Path to the executable
for (i in seq_along(fas_Fs_raw)) {
  cat("Processing",
      "-----------",
      i,
      "/",
      length(fas_Fs_raw),
      "-----------\n")
  system2(cutadapt, args = c(R1_flags, R2_flags,
                             "--discard-untrimmed",
                             "--max-n 0",
                             # Optional strong constraint on expected length
                             #paste0("-m ", 250-nchar(FWD)[1], ":", 250-nchar(REV)[1]),
                             #paste0("-M ", 250-nchar(FWD)[1], ":", 250-nchar(REV)[1]),
                             "-o", fas_Fs_cut[i], "-p", fas_Rs_cut[i],
                             fas_Fs_raw[i], fas_Rs_raw[i]))
}

out_1 <- cbind(ShortRead::qa(fas_Fs_raw)[["readCounts"]][, "read", drop = FALSE],
               ShortRead::qa(fas_Fs_cut)[["readCounts"]][, "read", drop = FALSE])

head(out_1)

###################################################

#### DADA2 Pipeline starts here ####

path_process <- path_cut
fas_Fs_process <- sort(list.files(path_process,
                                  pattern = file_pattern["F"],
                                  full.names = TRUE))
fas_Rs_process <- sort(list.files(path_process,
                                  pattern = file_pattern["R"],
                                  full.names = TRUE))

sample_names <- sapply(strsplit(basename(fas_Fs_process), "_"),
                                function(x) {
                                  paste0(x[c(1)], collapse = "_")})

#### Inspect read quality profiles ####
plotQualityProfile(fas_Fs_process[2])
plotQualityProfile(fas_Rs_process[2])

pdf(file.path(path_results, "Read_quality_profile_aggregated.pdf"))
p <- plotQualityProfile(sample(fas_Fs_process,
                                replace = FALSE,
                                size = ifelse(length(fas_Fs_process) < 100,
                                length(fas_Fs_process),
                                100)),
                                aggregate = TRUE)
p + ggplot2::labs(title = "Forward")
p <- plotQualityProfile(sample(fas_Rs_process,
                                replace = FALSE,
                                size = ifelse(length(fas_Rs_process) < 100,
                                length(fas_Rs_process),
                                100)),
                                aggregate = TRUE)
p + ggplot2::labs(title = "Reverse")
dev.off()


#### FILTER AND TRIM ####
fas_Fs_filtered <- file.path(path, "filtered", basename(fas_Fs_process))
fas_Rs_filtered <- file.path(path, "filtered", basename(fas_Rs_process))
all.equal(basename(fas_Fs_raw), basename(fas_Fs_filtered))

names(fas_Fs_filtered) <- sample_names
names(fas_Rs_filtered) <- sample_names

# Think twice before copying next command
# Check DADA2 official tutorial and the help of filterAndTrim
out_2 <- filterAndTrim(fas_Fs_process,
                        fas_Fs_filtered,
                        fas_Rs_process,
                        fas_Rs_filtered,
                        truncLen = c(200, 170),
                        maxN = 0,
                        maxEE = c(2, 2),
                        truncQ = 2,
                        rm.phix = TRUE,
                        compress = TRUE,
                        multithread = TRUE)
head(out_2)


#### LEARN THE ERROR RATES ####
error_F <- learnErrors(fas_Fs_filtered, multithread = TRUE, randomize = TRUE)
error_R <- learnErrors(fas_Rs_filtered, multithread = TRUE, randomize = TRUE)

pdf(file.path(path_results, "Error_rates_learning.pdf"))
p <- plotErrors(error_F, nominalQ = TRUE)
p + ggplot2::labs(title = "Error Forward")
p <- plotErrors(error_R, nominalQ = TRUE)
p + ggplot2::labs(title = "Error Reverse")
dev.off()


#### DEREPLICATION, SAMPLE INFERENCE & MERGE PAIRED READS ####

merged_list <- vector("list", length(sample_names))
names(merged_list) <- sample_names

for (i in sample_names) {
  cat("Processing -------",
      which(sample_names == i),
      "/",
      length(sample_names),
      "-------",
      i,
      "\n")
  derep_Fs <- derepFastq(fas_Fs_filtered[[i]], verbose = TRUE)
  derep_Rs <- derepFastq(fas_Rs_filtered[[i]], verbose = TRUE)
  dds_Fs <- dada(derep_Fs, err = error_F, multithread = TRUE, verbose = TRUE)
  dds_Rs <- dada(derep_Rs, err = error_R, multithread = TRUE, verbose = TRUE)
  merged_list[[i]] <- mergePairs(dds_Fs,
                                  derep_Fs,
                                  dds_Rs,
                                  derep_Rs,
                                  verbose = TRUE)
}

#### CONSTRUCT SEQUENCE TABLE ####

to.rm <- as.numeric(lapply(merged_list, function(x) sum(x$abundance)))
to.kp <- !(to.rm == 0)
merged_list_clean <- merged_list[to.kp]

seqtab <- makeSequenceTable(merged_list_clean)
dim(seqtab)
table(nchar(getSequences(seqtab)))

#### REMOVE CHIMERA ####

seqtab_nochim <- removeBimeraDenovo(seqtab,
                                    method = "consensus",
                                    multithread = TRUE,
                                    verbose = TRUE)
dim(seqtab_nochim)
table(nchar(getSequences(seqtab_nochim)))

#### TRACK READS THROUGH THE PIPELINE ####

reads.raw <- data.frame("sample" = rownames(seqtab),
                        "Merged" = rowSums(seqtab))
reads.nochim <- data.frame("sample" = rownames(seqtab_nochim),
                          "Nonchim" = rowSums(seqtab_nochim))

track <- data.frame("sample" = sample_names,
                    "Raw" = out_1[, 1],
                    "Cutadapt" = out_1[, 2],
                    "Filtered" = out_2[, 2])
track <- left_join(track, reads.raw) %>% left_join(reads.nochim, by = "sample")
track[is.na(track)] <- 0
track <- track[!(track$Nonchim > track$Raw), ]

#### ASSIGN TAXONOMY ####

# Here we download and use the Diat.barcode (last version) for DADA2.
# You can use your own local database if needed.
tax_fas <- diatbarcode::download_diatbarcode(flavor = "rbcl312_dada2_tax")
tax <- assignTaxonomy(seqtab_nochim, tax_fas$path, minBoot = 75,
                      taxLevels = c("Empire",
                                    "Kingdom",
                                    "Subkingdom",
                                    "Phylum",
                                    "Class",
                                    "Order",
                                    "Family",
                                    "Genus",
                                    "Species"),
                      outputBootstraps = TRUE,
                      verbose = TRUE,
                      multithread = TRUE)

spe_fas <- diatbarcode::download_diatbarcode(flavor = "rbcl312_dada2_spe")
exact_sp <- assignSpecies(seqtab_nochim, spe_fas$path)


#### CLEAN AND SAVE EVERYTHING ####

# I like to save in CSV :)
# I also transpose community matrices and turn rownames into columns for better integration with dplyr.
# Note that if you need to reload these files to process them later with DADA2 functions
# (eg. merge several runs), it is better to save raw DADA2 objects in .rds (see official tutorial).

prep_cdm <- function(x) {
  x <- t(x)
  x <- as.data.frame(x)
  asv_num <- vector(dim(x)[1], mode = "character")
  for (i in 1:dim(x)[1]) {
    asv_num[i] <- paste("ASV", i, sep = "_")
  }
  x <- cbind(rownames(x), asv_num, x)
  colnames(x)[1] <- "DNA_SEQ"
  colnames(x)[2] <- "ASV"
  return(x)
}

write.csv(prep_cdm(seqtab),
          file.path(path_results, "sequence_table.csv"),
          row.names = FALSE)
write.csv(prep_cdm(seqtab_nochim),
          file.path(path_results, "sequence_table_nochim.csv"),
          row.names = FALSE)

track <- as.data.frame(track)
track <- cbind(rownames(track), track)
colnames(track)[1] <- "Sample"
write.csv(track, file.path(path_results, "track_reads.csv"), row.names = FALSE)

# giving our seq headers more manageable names (ASV_1, ASV_2...)
asv_seqs <- colnames(seqtab_nochim)
asv_headers <- vector(dim(seqtab_nochim)[2], mode = "character")

for (i in 1:dim(seqtab_nochim)[2]) {
  asv_headers[i] <- paste(">ASV", i, sep = "_")
}

# making and writing out a fasta of our final ASV seqs:
asv_fasta <- c(rbind(asv_headers, asv_seqs))
write(asv_fasta, file.path(path_results, "ASVs.fa"))

tax <- as.data.frame(tax)
ASV <- gsub(pattern = ">", replacement = "", x = asv_headers)
tax <- cbind(rownames(tax), ASV, tax)
colnames(tax)[1] <- "DNA_SEQ"
colnames(tax) <- sub("^tax\\.", "", colnames(tax))
colnames(tax) <- sub("^boot\\.", "BOOT_", colnames(tax))
write.csv(tax,
          file.path(path_results, "seq_nochim_tax.csv"),
          row.names = FALSE)

exact_sp <- as.data.frame(exact_sp)
exact_sp <- cbind(rownames(exact_sp), exact_sp)
colnames(exact_sp)[1] <- "DNA_SEQ"
colnames(tax)[2] <- "ASV"
write.csv(exact_sp,
          file.path(path_results, "seq_nochim_exact_sp.csv"),
          row.names = FALSE)

##Removing likely contaminants

library(decontam)

# count table:
asv_tab <- prep_cdm(seqtab_nochim)
rownames(asv_tab) <- seq(from = 1, to = dim(asv_tab)[1], by = 1)

#Create bolean vector of blanks
vector_for_decontam <- c()

#Replace neg by the character identifier for your blank samples
for (i in 1:dim(asv_tab[, 3:dim(asv_tab)[2]])[2]) {
  vector_for_decontam[i] <- str_detect(colnames(asv_tab[, 3:dim(asv_tab)[2]])[i],
                                       "blank_identifier",
                                       negate = FALSE)
}

contam_df <- isContaminant(t(asv_tab[, 3:dim(asv_tab)[2]]),
                            neg = vector_for_decontam)

table(contam_df$contaminant)

# getting vector holding the identified contaminant IDs
contam_asvs <- row.names(contam_df[contam_df$contaminant == TRUE, ])
asv_tax[contam_asvs, ]

# making new fasta file
contam_indices <- which(asv_fasta %in% paste0(">ASV_", contam_asvs))
dont_want <- sort(c(contam_indices, contam_indices + 1))
asv_fasta_no_contam <- asv_fasta[- dont_want]

# making new count table
asv_tab_no_contam <- asv_tab[!asv_tab$ASV %in% paste0("ASV_", contam_asvs), ]

# making new taxonomy table
asv_tax_no_contam <- asv_tax[!asv_tax$ASV %in% paste0("ASV_", contam_asvs), ]

dim(asv_tab_no_contam)
dim(asv_tax_no_contam)
## and now writing them out to files
write(asv_fasta_no_contam,
      file.path(path_results,
      "ASVs-no-contam.fa"))
write.csv(asv_tab_no_contam,
          file.path(path_results,
          "ASVs_counts-no-contam.csv"),
          row.names = FALSE)
write.csv(asv_tax_no_contam,
          file.path(path_results,
          "ASVs_taxonomy-no-contam.csv"),
          row.names = FALSE)
