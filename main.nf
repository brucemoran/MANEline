#!/usr/bin/env nextflow

def helpMessage() {
  log.info"""
  ------------------------------------------------------------------------------
                            MANEline PIPELINE
  ------------------------------------------------------------------------------
  Usage:

  nextflow run brucemoran/MANEline

  Mandatory arguments:

    -profile        [str]       Configuration profile (required: singularity)

    --assembly      [str]       Assembly of bedFile, either GRCh37 or GRCh38 (default)

    --mane_version  [str]       Numeric version of MANE to use (default:
                                current; please only supply numeric e.g. 0.91,
                                not release_0.91)

    --bedFile       [str]       Path to BED file of regions to annotate with
                                MANE (ENS, REF transcript ID, gene name). If none given, all 'gene' entires will be used to create a 'full' BED of MANE transcripts and their coordinates and annotations.

    --email         [str]       Email address to send reports

    """.stripIndent()
}

if (params.help) exit 0, helpMessage()
if (params.assembly != "GRCh38" && params.assembly != "GRCh37"){
  exit 1, "Please define --assembly as either GRCh37 or GRCh38"
}
process Download {

  label 'process_low'
  publishDir "${params.outDir}/downloads", mode: "copy"

  output:
  tuple file("*.ensembl_genomic.gtf.gz"), file("*.summary.txt.gz") into ( bed_gtf, liftover )
  file('vers.txt') into vers_get
  file('grch_vers.txt') into grch_get

  script:
  def mane_base = params.mane_base
  def mane_vers = params.mane_version == "current" ? "current" : "release_" + params.mane_version
  def grch_vers = params.assembly
  """
  VERS=\$(curl -l ftp.ncbi.nlm.nih.gov/refseq/MANE/MANE_human/${mane_vers}/ | grep summary | cut -d "." -f3,4)
  wget ${mane_base}/${mane_vers}/MANE.GRCh38.\$VERS.ensembl_genomic.gtf.gz
  wget ${mane_base}/${mane_vers}/MANE.GRCh38.\$VERS.summary.txt.gz
  echo \$VERS > vers.txt
  echo ${grch_vers} > grch_vers.txt
  """
}

vers_get.map { it.text.strip() }.set{ vers_mane }
grch_get.map { it.text.strip() }.set{ grch_vers }

process GtfBed {

  label 'process_low'
  publishDir "${params.outDir}/gtf", mode: "copy"

  input:
  tuple file(gtf_gz), file(sum_gtf) from bed_gtf
  val(vers) from vers_mane
  val(grch_vers) from grch_vers

  output:
  file("${grch_vers}.MANE.${vers}.gtf.bed") into lift_bed
  val(vers) into vers_mane_1
  val(grch_vers) into grchvers_1

  script:
  """
  ##create MANE bedfile
  gunzip -c ${gtf_gz} | sed 's/\"//g' | sed 's/;//g' | \\
    perl -ane 'chomp; if(\$F[2] eq "transcript"){
    print "\$F[0]\\t\$F[3]\\t\$F[4]\\t\$F[15];\$F[9];\$F[11];\$F[23];\$F[25]\\n";}' | \\
    sed 's/RefSeq://g' > ${grch_vers}.MANE.${vers}.gtf.bed
  """
}

process Liftover {

  label 'process_low'
  publishDir "${params.outDir}/liftover", mode: "copy"

  input:
  val(vers) from vers_mane_1
  file(bed) from lift_bed
  val(grch_vers) from grchvers_1

  output:
  file("${grch_vers}.lift.MANE.${vers}.gtf.bed") into lifted_bed
  val(vers) into vers_mane_2
  val(grch_vers) into grchvers_2

  script:
  def hgTohg = "hg38ToHg19"
  def hg = "hg19"
  if(params.assembly == "GRCh37")
    """
    ##lift
    echo "nameserver 8.8.8.8" > /tmp/resolv.conf
    wget http://hgdownload.cse.ucsc.edu/goldenPath/hg38/liftOver/${hgTohg}.over.chain.gz

    liftOver ${bed} ${hgTohg}.over.chain.gz ${grch_vers}.lift.MANE.${vers}.gtf.bed unmapped
    """
  else
    """
    cp ${bed} ${grch_vers}.lift.MANE.${vers}.gtf.bed
    """
}

if( params.bedFile != null ){
  Channel.fromPath( "${params.bedFile}" ).set( bed_file )

  process LiftOverlap {
    label 'process_low'
    publishDir "${params.outDir}/liftover", mode: "copy"

    input:
    val(vers) from vers_mane_2
    file(bed_lift) from lifted_bed
    file(bed_over) from bed_file
    val(grch_vers) from grchvers_2

    output:
    file("${grch_vers}.lift.overlap.MANE.${vers}.gtf.bed") into complete

    script:
    """
    ##overlap
    perl ${workflow.projectDir}/assets/pover.pl ${bed_lift} ${bed_over} ${grch_vers}.lift.overlap.MANE.${vers}.gtf.bed
    """
  }
} else {

  process Finish {
    label 'process_low'
    publishDir "${params.outDir}/liftover", mode: "copy"

    input:
    val(vers) from vers_mane_2
    file(bed_lift) from lifted_bed
    val(grch_vers) from grchvers_2

    output:
    file("${grch_vers}.lift.overlap.MANE.${vers}.gtf.bed") into complete

    script:
    """
    ##overlap
    perl ${workflow.projectDir}/assets/pover.pl ${bed_lift} ${bed_lift} ${grch_vers}.lift.overlap.MANE.${vers}.gtf.bed
    """
  }
}
