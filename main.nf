#!/usr/bin/env nextflow

nextflow.enable.dsl=2

/*********************************************************************************
 * PARAMETERS
 *********************************************************************************/
params.ENST_list  = 'data/Oncomine_Focus_designs_063020_Reference.enst.txt'  // Path to the file containing old GRCh37 ENST IDs (one per line)
params.mane_gtf   = 'data/MANE.GRCh38.v1.2.ensembl_genomic.gtf.gz' // Path to the MANE GTF file
params.outdir     = 'results'                    // Output directory
params.batchSize  = 50                           // Example: how many IDs to query at once from the REST API (optional tweak)

/*********************************************************************************
 * WORKFLOW DEFINITION
 *********************************************************************************/
workflow {

    // Channels: input files
    ch_mane_gtf = Channel.fromPath(params.mane_gtf)
    ch_grch37_ids = Channel.fromPath(params.ENST_list)

    // Step 1: Parse the MANE GTF to build Ensembl->RefSeq mapping
    parseManeGtf(ch_mane_gtf)

    // Step 2: Convert old GRCh37 ENST to current GRCh38 IDs (via Ensembl REST)
    convertIds(ch_grch37_ids)

    // Step 3: Join the updated GRCh38 ENST with the MANE mapping to get RefSeq
    joinResults(parseManeGtf.out, convertIds.out)
}


/*********************************************************************************
 * MODULES / PROCESSES
 *********************************************************************************/

/**
 * parseManeGtf
 * Parses the MANE GTF file to extract a mapping of:
 *   ENST<version> -> RefSeq<version>
 * Output is a TSV with columns: ensembl_id, refseq_id
 */
process parseManeGtf {

    publishDir params.outdir, mode: 'copy'

    // The input channel is the single GTF file
    input:
    path gtf_file

    // The output is a single TSV
    output:
    path 'mane_mapping.tsv'


    script:
    """
    Rscript ${workflow.projectDir}/assets/parse_mane.R \\
        --gtf_file ${gtf_file} \\
        --output mane_mapping.tsv
    """
}


/**
 * convertIds
 * Reads a list of old (GRCh37) Ensembl transcript IDs and uses the Ensembl REST API
 * to fetch the *current* stable ID on GRCh38 (if it still exists).
 * Output is a TSV with columns: old_id, new_id
 */
process convertIds {

    publishDir params.outdir, mode: 'copy'

    // The input channel is the single text file with old ENST IDs
    input:
    path enst_list_file

    // The output is a single TSV
    output:
    path 'converted_ids.tsv'

    script:
    """
    Rscript ${workflow.projectDir}/assets/convert_ids.R \\
        --input_ids ${enst_list_file} \\
        --output converted_ids.tsv \\
        --batchSize ${params.batchSize}
    """
}


/**
 * joinResults
 * Joins the updated (GRCh38) ENST IDs from convertIds with
 * the mane_mapping from parseManeGtf to find the matched RefSeq ID.
 * Produces a final TSV with columns:
 *   old_enst, updated_enst, refseq_id, status
 */
process joinResults {

    publishDir params.outdir, mode: 'copy'
    input:
    path mane_map_tsv
    path converted_ids_tsv

    output:
    path 'final_mapping.tsv'


    script:
    """
    Rscript -e ' \
      library(dplyr); \
      mane_map <- read.delim(\"${mane_map_tsv}\", header=TRUE, sep=\"\\t\"); \
      converted <- read.delim(\"${converted_ids_tsv}\", header=TRUE, sep=\"\\t\"); \
      # Left join on updated_enst = ensembl_id
      result <- converted %>% left_join(mane_map, by=c(\"new_id\"=\"ensembl_id\")); \
      # If refseq_id is NA, then no direct MANE match
      result\$status <- ifelse(is.na(result\$refseq_id), \"No direct MANE match\", \"MANE match\"); \
      write.table(result, file=\"final_mapping.tsv\", sep=\"\\t\", quote=FALSE, row.names=FALSE); \
    '
    """
}
