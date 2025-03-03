#!/usr/bin/env nextflow

nextflow.enable.dsl=2


def helpMessage() {
  log.info"""
  ------------------------------------------------------------------------------
                            MANEline PIPELINE
  ------------------------------------------------------------------------------
  Usage:

  nextflow run brucemoran/MANEline

  Description:

  Convert MANE GRCh38 annotation to BED format, and lift backwards to GRCh37.
  Optionally supply a BED file for which exon-level annotations are returned at
  GRCh37 and GRCh38 coordinates.

  Mandatory arguments:

    -profile        [str]     Configuration profile (required: singularity)

    --runID         [str]     Identifier for the run

    --mane_version  [str]     Numeric version of MANE to use (default:
                              current; please only supply numeric e.g. 0.91,
                              not release_0.91)

    --email         [str]     Email address to send reports

    --bedFile       [str]     Path to BED file of regions to annotate with
                              MANE. Default: null.

    --bedAssembly   [str]     One of GRCh37, GRCh38 indicating the assembly used in --bedFile.
                              Default: null.

    --bedOtherAss   [bool]    Return the bedFile with the other assembly lifted over Default: true.
    """.stripIndent()
}

if (params.help) exit 0, helpMessage()
//require runID
if(params.runID == null){
  exit 1, "Please specify --runID"
}

/*********************************************************************************
 * WORKFLOW DEFINITION
 *********************************************************************************/
workflow {

    // Channels in
    ch_bedFile =Channel.fromPath(params.bedFile) 

    // Step 1: DL MANE GTF, Summary
    DLMANE(params.mane_base, params.mane_version)

    // Step 2: 
    GTFBED(DLMANE.out)

    // Step 3: Join the updated GRCh38 ENST with the MANE mapping to get RefSeq
    LIFTOVER(GTFBED.out)

    OVERLAP37(LIFTOVER.out, ch_bedFile)

    OTHER37(OVERLAP37.out)

    REPORT37(OTHER37.out)

}


/*********************************************************************************
 * MODULES / PROCESSES
 *********************************************************************************/

//Get GRCh38 (only assembly available) MANE GTF and summary files based on input versions required
process DLMANE {

    publishDir "${params.outDir}/download", mode: "copy"

    input:
    val mane_base
    val mane_vers

    output:
    tuple path("*.ensembl_genomic.gtf.gz"), path("*.summary.txt.gz"), val(mane_vers)

    script:
    mane_base = params.mane_base
    mane_vers = params.mane_version == "current" ? "current" : "release_" + params.mane_version
    """
    VERS=\$(curl -l ftp.ncbi.nlm.nih.gov/refseq/MANE/MANE_human/${mane_vers}/ | grep summary | cut -d "." -f3,4)
    wget ${mane_base}/${mane_vers}/MANE.GRCh38.\$VERS.ensembl_genomic.gtf.gz
    wget ${mane_base}/${mane_vers}/MANE.GRCh38.\$VERS.summary.txt.gz
    echo \$VERS > vers.txt
    """
}

//parse the GTF to BED
process GTFBED {

  label 'process_low'
  publishDir "${params.outDir}/bed", mode: "copy"

  input:
  tuple path(gtf_gz), path(sum_gtf), val(vers)

  output:
  tuple path("GRCh38.MANE.${vers}.exon.bed"), val(vers)

  script:
  """
  perl ${workflow.projectDir}/assets/pexon.pl ${gtf_gz} GRCh38.MANE.${vers}.exon.bed
  """
}

//Liftover MANE to GRCh37
process LIFTOVER {
  publishDir "${params.outDir}/bed", mode: "copy"

  input:
  tuple path(exon_bed), val(vers)

  output:
  tuple path("GRCh37.MANE.${vers}.exon.bed"), val(vers)

  script:
  """
  ##lift
  echo "nameserver 8.8.8.8" > /tmp/resolv.conf
  wget http://hgdownload.cse.ucsc.edu/goldenPath/hg38/liftOver/hg38ToHg19.over.chain.gz
  liftOver ${exon_bed} hg38ToHg19.over.chain.gz GRCh37.MANE.${vers}.exon.bed unmapped
  """
}

//if a BED is supplied overlapp the created Liftover
process OVERLAP37 {
    publishDir "${params.outDir}/bed", mode: "copy"

    input:
    tuple path(exon_lift), val(vers)
    path(bed_user)

    output:
    tuple path("GRCh37.MANE.${vers}.exon.overlap.${bed_user}"), path(bed_user), path(exon_lift)

    script:
    """
    ##if bed_user overlaps, return all exons for the gene overlapped
    perl ${workflow.projectDir}/assets/pover.pl ${exon_lift} ${bed_user}  GRCh37.MANE.${vers}.exon.overlap.${bed_user}
    """
}

process OTHER37 {
    publishDir "${params.outDir}/bed", mode: "copy"

    input:
    tuple path(bed_ass), path(bed_in), path(exon_lift)

    output:
    tuple path(bed_ass), path(bed_in), path(exon_lift), path("GRCh38*.${bed_in}")

    script:
    def bedname = "${bed_ass}".replace('GRCh37','GRCh38')
    """
    ##count total fields as liftOver needs to know this
    LC=\$(head -n5 ${bed_ass} | tail -n1 | awk -F"\\t" '{ print NF-2 }')
    echo "nameserver 8.8.8.8" > /tmp/resolv.conf
    wget http://hgdownload.cse.ucsc.edu/goldenPath/hg19/liftOver/hg19ToHg38.over.chain.gz
    liftOver -bedPlus=\$LC ${bed_ass} hg19ToHg38.over.chain.gz ${bedname} unmapped
    """
}

process REPORT37 {
    publishDir "${params.outDir}/bed", mode: "copy"

    input:
    tuple path(bed_ass), path(bed_in), path(exon_37), path(exon_38)

    output:
    path("${params.runID}.exns.csv")

    script:
    """
    ##parse into table format desired
    perl ${workflow.projectDir}/assets/ptabl.pl ${exon_37} ${exon_38} ${params.runID}.exns.csv
    """
}