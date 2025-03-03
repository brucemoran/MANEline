#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(optparse)
  library(stringr)
})

option_list <- list(
  make_option("--gtf_file", type="character", help="Path to MANE GTF"),
  make_option("--output", type="character", default="mane_mapping.tsv", help="Output TSV")
)
opt <- parse_args(OptionParser(option_list=option_list))

gtf_path  <- opt$gtf_file
out_tsv   <- opt$output

# We'll read just the transcript lines from the GTF
# GTF format: 9th column has attributes like: transcript_id "ENST..."; db_xref "RefSeq:NM_xxx"; ...
# This is a simple example that uses grep, so be mindful for large GTF files.

raw_gtf <- readLines(gtf_path)
transcript_lines <- grep("\\ttranscript\\t", raw_gtf, value=TRUE)

# We'll parse out transcript_id and RefSeq from the attributes
# Simplistic approach, you might want to use a real GTF parser
ensembl_ids <- c()
refseq_ids  <- c()

for(line in transcript_lines) {
  fields <- strsplit(line, "\t")[[1]]
  attrs  <- fields[9]
  # Extract transcript_id
  enst_match <- str_match(attrs, 'transcript_id \"(ENST[0-9\\.]+)\"')[,2]
  # Extract the RefSeq ID from db_xref=RefSeq:...
  refseq_match <- str_match(attrs, 'RefSeq:([^,;"]+)')[,2]  # e.g. NM_000546.6
  if(!is.na(enst_match)) {
    ensembl_ids <- c(ensembl_ids, enst_match)
    if(is.na(refseq_match)) {
      refseq_ids <- c(refseq_ids, NA)
    } else {
      refseq_ids <- c(refseq_ids, refseq_match)
    }
  }
}

mane_table <- data.frame(ensembl_id=ensembl_ids, refseq_id=refseq_ids, stringsAsFactors=FALSE)
mane_table <- unique(mane_table)  # remove duplicates

# Write output
write.table(mane_table, file=out_tsv, sep="\t", quote=FALSE, row.names=FALSE)
