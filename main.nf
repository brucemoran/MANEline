#!/usr/bin/env nextflow

def helpMessage() {
  log.info"""
  ------------------------------------------------------------------------------
                            MANEline PIPELINE
  ------------------------------------------------------------------------------
  Usage:

  nextflow run brucemoran/MANEline

  Mandatory arguments:

    -profile        [str]       Configuration profile

    --assembly      [str]       Either GRCh37 or GRCh38 (default)

    --mane_version  [str]       Numeric version of MANE to use (default:
                                current; please only supply numeric e.g. 0.91,
                                not release_0.91)

    --bedFile       [str]       Path to BED file of regions to annotate with
                                MANE (ENS, REF transcript ID, gene name)

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
  tuple file("*.ensembl_genomic.gtf.gz"), file("*.summary.txt.gz") into liftover
  file('vers.txt') into vers_get

  script:
  def mane_base = params.mane_base
  def mane_vers = params.mane_version == "current" ? "current" : "release_" + params.mane_version
  def grch_vers = params.assembly
  """
  VERS=\$(curl -l ftp.ncbi.nlm.nih.gov/refseq/MANE/MANE_human/${mane_vers}/ | grep summary | cut -d "." -f3,4)
  wget ${mane_base}/${mane_vers}/MANE.GRCh38.\$VERS.ensembl_genomic.gtf.gz
  wget ${mane_base}/${mane_vers}/MANE.GRCh38.\$VERS.summary.txt.gz
  echo \$VERS > vers.txt
  """
}

vers_get.map { it.text.strip() }.set{ vers_mane }
if( params.bedFile != null ){
  Channel.fromPath( "${params.bedFile}" ).set( bed_file )
} else {
  process FastaBed {

    label 'process_low'
    publishDir "${params.outDir}/downloads", mode: "copy"

    output:
    file('*.bed') into bed_file
    file("README_GRCh_bed.txt") into readme

    script:
    def grch_vers = params.assembly
    def dlink = params.assembly == "GRCh37" ? "http://ftp.ensembl.org/pub/grch37/current/fasta/homo_sapiens/dna/Homo_sapiens.GRCh37.dna.toplevel.fa.gz" : "http://ftp.ensembl.org/pub/current_fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.toplevel.fa.gz"
    """
    wget -O fasta.gz ${dlink}
    gunzip -c fasta.gz | \\
      grep ">" | \\
      grep -v ">[GH]" | \\
      perl -ane 'chomp; @s=split(/:/); print "\$s[4]\\t\$s[5]\\t\$s[3]\\n";' \\
      > ${grch_vers}.bed
    echo -e "Command to download fasta from which BED was made:\\nwget ${dlink}" > README_GRCh_bed.txt
    """
  }
}

process Bedparse {

  label 'process_low'
  publishDir "${params.outDir}/bed", mode: "copy"

  input:
  file(bed) from bed_file

  output:
  file("parsed.bed") into bed_parse

  script:
  """
  ##remove any non-chr, coord lines in top of file
  CHR=\$(tail -n1 ${bed} | perl -ane 'print \$F[0];')
  if [[ \$CHR =~ "chr" ]]; then
    perl -ane 'if(\$F[0]=~m/^chr/){print \$_;}' ${bed} > parsed.bed
  else
    perl -ane 'if(\$F[0]=~m/^[0-9MXY]/){print \$_;}' ${bed} > parsed.bed
  fi
  """
}

process Liftover {

  conda 'ucsc-liftover'
  label 'process_low'
  publishDir "${params.outDir}/liftover", mode: "copy"

  input:
  tuple file(gtf_gz), file(sum_gz) from liftover
  val(vers_mane) from vers_mane
  file(bed) from bed_parse

  output:
  file("*") into complete

  script:
  def grch_vers = params.assembly
  def hgTohg = params.assembly == "GRCh37" ? "hg19ToHg38" : null
  def hg = params.assembly == "GRCh37" ? "hg19" : null
  """
  wget http://hgdownload.cse.ucsc.edu/goldenPath/${hg}/liftOver/${hgTohg}.over.chain.gz
  liftOver ${bed} ${hgTohg}.over.chain.gz ${grch_vers}.lift.bed unmapped

  echo -r"Liftover using:\\nwget http://hgdownload.cse.ucsc.edu/goldenPath/${hg}/liftOver/${hgTohg}.over.chain.gz
  liftOver ${bed} ${hgTohg}.over.chain.gz ${grch_vers}.lift.bed unmapped" >> README.txt
  """
}
