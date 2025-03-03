#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(httr)
  library(jsonlite)
  library(optparse)
})

option_list <- list(
  make_option("--input_ids", type="character", help="File with old GRCh37 ENST IDs, one per line"),
  make_option("--output", type="character", default="converted_ids.tsv", help="Output TSV"),
  make_option("--batchSize", type="integer", default=50, help="Number of IDs to query per batch")
)
opt <- parse_args(OptionParser(option_list=option_list))

id_file   <- opt$input_ids
out_tsv   <- opt$output
batchSize <- opt$batchSize

old_ids <- readLines(id_file)
old_ids <- old_ids[old_ids != ""]  # remove empty lines

get_current_id <- function(ensembl_id) {
  # Query the Ensembl REST API for the current stable ID
  url <- paste0("https://rest.ensembl.org/archive/id/", ensembl_id, "?")
  resp <- GET(url, add_headers(`Content-Type`="application/json"))
  if (resp$status_code == 200) {
    json <- fromJSON(content(resp, as="text"))
    return(json$id)  # 'id' is the current stable ID
  } else {
    return(NA)
  }
}

results <- data.frame(old_id=character(), new_id=character(), stringsAsFactors=FALSE)

# We'll do a simple batching approach
n_ids <- length(old_ids)
idx_start <- seq(1, n_ids, by=batchSize)

for (start_i in idx_start) {
  end_i <- min(start_i + batchSize - 1, n_ids)
  batch <- old_ids[start_i:end_i]

  for (oid in batch) {
    cid <- get_current_id(oid)
    results <- rbind(results, data.frame(old_id=oid, new_id=ifelse(is.null(cid), NA, cid), stringsAsFactors=FALSE))
  }
  message(sprintf("Processed %d / %d IDs...", end_i, n_ids))
  Sys.sleep(0.5) # small pause to avoid overwhelming the server
}

# Write out final table
write.table(results, file=out_tsv, sep="\t", row.names=FALSE, quote=FALSE)
