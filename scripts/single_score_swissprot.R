#!/usr/bin/env Rscript
# Ensure R uses the conda environment library
.libPaths(file.path(Sys.getenv("CONDA_PREFIX"), "Lib", "R", "library"))

# load required packages
library(epictope)
# ensure dplyr is available for data manipulation verbs used below
suppressPackageStartupMessages(library(dplyr))
rm(list = ls())

# helper functions
setup_files()
check_config()
args <- commandArgs(trailingOnly = TRUE); query <- args[1]
# download uniprot information 
uniprot_fields <- c("accession", "id", "gene_names", "xref_alphafolddb", "sequence", "organism_name", "organism_id")
uniprot_data <- query_uniProt(query = query, fields = uniprot_fields)

####
# download associated alphafold2 pdb
if(is.na(uniprot_data$AlphaFoldDB)){
  message("WARNING. NO CROSSREFERENCE ALPHAFOLD ENTRY FOUND; ATTEMPTING DIRECT LOOKUP. RESULTS MAY BE LOW CONFIDENCE.")
  uniprot_data$AlphaFoldDB <- query}
alphafold_file <- fetch_alphafold(gsub(";", "", uniprot_data$AlphaFoldDB))
if(is.na(alphafold_file)){stop("No alphafold file found for ", query)}
# calculate dssp on alphafold pdb file
dssp_res <- dssp_command(alphafold_file)
# parse and read in dssp
dssp_df <- parse_dssp(dssp_res)

# retrieve iupred/anchor2 disordered binding regions
iupred_df <- iupredAnchor(query)

# blast query aa sequence
seq <- Biostrings::AAStringSet(uniprot_data$Sequence, start=NA, end=NA, width=NA, use.names=TRUE)


# blast
blast_res <- protein_blast(seq, file.path(cds_folder, blast_db))

# get top protein per species
filtered_hits <- blast_res %>%
  group_by(staxid) %>%
  arrange(desc(bitscore), .by_group = TRUE) %>%
  slice(1) %>%
  ungroup()

entry_file <- tempfile("entries_", fileext = ".txt")
writeLines(unique(filtered_hits$sseqid), entry_file)

# fetch the AA sequence for blast best match
aa_seqs_text <- system2(
  "blastdbcmd",
  args = c(
    "-db", file.path(cds_folder, blast_db),
    "-entry_batch", entry_file,
    "-outfmt", "%f"
  ),
  stdout = TRUE
)

# clean names of sequences
aa_seqs_clean <- gsub("\\..*", "", aa_seqs_text)
aa_file <- tempfile("hits_", fileext = ".fasta")
writeLines(aa_seqs_clean, aa_file)
blast_stringset <- Biostrings::readAAStringSet(aa_file, format = "fasta")

# add query to list if missing
if (!(query %in% names(blast_stringset))) {
  blast_seqs <- as.character(blast_stringset)
  blast_seqs[query] <- seq
  blast_stringset <- Biostrings::AAStringSet(blast_seqs)
}

# multiple sequence alignment
msa_res <-  muscle(blast_stringset)
# write msa to file
# convert to XStringSet
msa_res <- Biostrings::AAStringSet(msa_res)
Biostrings::writeXStringSet(msa_res, file = paste0(outputFolder, "/", query, "_msa.fasta"))
# shannon entropy calculation
shannon_df <- shannon_reshape(msa_res, query)

# join tagging features in dataframe
features_df <- Reduce(function(x, y) merge(x, y, all=TRUE), list(shannon_df, dssp_df, iupred_df), accumulate=FALSE)
colnames(features_df)
# normalize features and calculate tagging score
norm_feats_df <- calculate_scores(features_df)

# write to file.
res_df <- merge(norm_feats_df, features_df)
write.csv(apply(res_df,2,as.character), file = paste0(outputFolder, "/", query, "_score.csv"), row.names = FALSE)
