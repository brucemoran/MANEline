#!/usr/bin/env nextflow

def helpMessage() {
  log.info"""
  ------------------------------------------------------------------------------
                            MANEline PIPELINE
  ------------------------------------------------------------------------------
  Usage:

  nextflow run brucemoran/MANEline

  Description:

  Convert MANE GRCh38 annotation backwards to GRCH37. Allows supply of a BED file for which annotations are parsed and included in output BED.

  Mandatory arguments:

    -profile        [str]       Configuration profile (required: singularity)

    --mane_version  [str]       Numeric version of MANE to use (default:
                                current; please only supply numeric e.g. 0.91,
                                not release_0.91)

    --bedFile       [str]       Path to BED file of regions to annotate with
                                MANE (ENS, REF transcript ID, gene name). Default: null.

    --bedAssembly   [str]       One of GRCh37, GRCh38 indicating the assembly used in --bedFile.

    --combine       [bool]      Do you want to combine transcript and exon BEDs to a 'txps_exon' BED with both transcript and exon?

    --email         [str]       Email address to send reports

    """.stripIndent()
}

if (params.help) exit 0, helpMessage()

//Get GRCh38 MANE GTF and summary files based on input versions required
process Download {

  label 'process_low'
  publishDir "${params.outDir}/download", mode: "copy"

  output:
  tuple file("*.ensembl_genomic.gtf.gz"), file("*.summary.txt.gz") into ( bed_gtf, liftover )
  file('vers.txt') into vers_get

  script:
  def mane_base = params.mane_base
  def mane_vers = params.mane_version == "current" ? "current" : "release_" + params.mane_version
  """
  VERS=\$(curl -l ftp.ncbi.nlm.nih.gov/refseq/MANE/MANE_human/${mane_vers}/ | grep summary | cut -d "." -f3,4)
  wget ${mane_base}/${mane_vers}/MANE.GRCh38.\$VERS.ensembl_genomic.gtf.gz
  wget ${mane_base}/${mane_vers}/MANE.GRCh38.\$VERS.summary.txt.gz
  echo \$VERS > vers.txt
  """
}

//map version and assembly
vers_get.map { it.text.strip() }.set{ vers_mane }

//parse the GTF to BED
process GtfBed {

  label 'process_low'
  publishDir "${params.outDir}/bed", mode: "copy"

  input:
  tuple file(gtf_gz), file(sum_gtf) from bed_gtf
  val(vers) from vers_mane

  output:
  file("GRCh38.MANE.${vers}.transcript.bed") into ( txp_bed, just_txp_bed )
  file("GRCh38.MANE.${vers}.exon.bed") into ( exon_bed, just_exon_bed )
  val(vers) into vers_mane_1

  script:
  """
  ##create MANE bedfiles of transcripts (inc. gene, protein IDs), and exons (link by txp)
  gunzip -c ${gtf_gz} | sed 's/\\"//g' | sed 's/;//g' | \\
    perl -ane 'chomp; if(\$F[2] eq "transcript"){
    print "\$F[0]\\t\$F[3]\\t\$F[4]\\t\$F[15];\$F[9];\$F[11];\$F[23];\$F[25]\\n";}' | \\
    sed 's/RefSeq://g' > GRCh38.MANE.${vers}.transcript.bed

  gunzip -c ${gtf_gz} | sed 's/\\"//g' | sed 's/;//g' | \\
    perl -ane 'chomp; if(\$F[2] eq "exon"){
    print "\$F[0]\\t\$F[3]\\t\$F[4]\\t\$F[11];\$F[23];exon_\$F[21]\\n";}' | \\
    sed 's/RefSeq://g' >> GRCh38.MANE.${vers}.exon.bed
  """

}

//Liftover MANE to GRCh37
process Liftover {

  label 'process_low'
  publishDir "${params.outDir}/bed", mode: "copy"

  input:
  file(bed_txp) from txp_bed
  file(bed_exon) from exon_bed
  val(vers) from vers_mane_1

  output:
  file("GRCh37.MANE.${vers}.transcript.bed") into lifted_txp_bed
  file("GRCh37.MANE.${vers}.exon.bed") into lifted_exon_bed
  val(vers) into vers_mane_2

  script:
  """
  ##lift
  echo "nameserver 8.8.8.8" > /tmp/resolv.conf
  wget http://hgdownload.cse.ucsc.edu/goldenPath/hg38/liftOver/hg38ToHg19.over.chain.gz

  liftOver ${bed_txp} hg38ToHg19.over.chain.gz GRCh37.MANE.${vers}.transcript.bed unmapped
  liftOver ${bed_exon} hg38ToHg19.over.chain.gz GRCh37.MANE.${vers}.exon.bed unmapped
  """
}

//if a BED is supplied overlapp the created Liftover
if( params.bedFile != null ){
  if(params.bedAssembly == "GRCh38" ){
    process JustOverlap {
      label 'process_low'
      publishDir "${params.outDir}/bed", mode: "copy"

      input:
      file(txp_just) from just_txp_bed
      file(exon_just) from just_exon_bed
      file(bed_over) from Channel.fromPath( "${params.bedFile}" )
      val(vers) from vers_mane_2

      output:
      tuple file("*.overlap.MANE.${vers}.transcript.bed"), file("*.overlap.MANE.${vers}.exon.bed") into complete

      script:
      def bedname = "${bed_over}".split('\\.')[0]
      """
      ##overlap
      perl ${workflow.projectDir}/assets/pover.pl ${txp_just} ${bed_over} 1
      uniq 1 > ${bedname}.GRCh38.overlap.MANE.${vers}.transcript.bed
      perl ${workflow.projectDir}/assets/pover.pl ${exon_just} ${bed_over} 1
      uniq 1 > ${bedname}.GRCh38.overlap.MANE.${vers}.exon.bed
      """
    }
  } else {
    process LiftOverlap {
      label 'process_low'
      publishDir "${params.outDir}/bed", mode: "copy"

      input:
      file(txp_lift) from lifted_txp_bed
      file(exon_lift) from lifted_exon_bed
      file(bed_over) from Channel.fromPath( "${params.bedFile}" )
      val(vers) from vers_mane_2

      output:
      tuple file("*.overlap.MANE.${vers}.transcript.bed"), file("*.overlap.MANE.${vers}.exon.bed") into complete

      script:
      def bedname = "${bed_over}".split('\\.')[0]
      """
      ##overlap
      perl ${workflow.projectDir}/assets/pover.pl ${txp_lift} ${bed_over} 1
      uniq 1 > ${bedname}.GRCh37.overlap.MANE.${vers}.transcript.bed
      perl ${workflow.projectDir}/assets/pover.pl ${exon_lift} ${bed_over} 1
      uniq 1 > ${bedname}.GRCh37.overlap.MANE.${vers}.exon.bed
      """
    }
  }

  if( params.combine != null ){
    process CombineET {
      label 'process_low'
      publishDir "${params.outDir}/bed", mode: "copy"

      input:
      tuple file(txp_com), file(exon_com) from complete

      output:
      file('*') into completed

      script:
      def comname = "${exon_com}".replace('exon', 'txps_exon')[0]
      """
      ##overlap
      perl ${workflow.projectDir}/assets/pcomb.pl ${txp_com} ${exon_com} 1
      uniq 1 > ${comname}
      """
    }
  }
}

if( params.email != null ){
  workflow.onComplete {
    sleep(100)
    def subject = """\
      [brucemoran/MANEline] SUCCESS [$workflow.runName]
      """
      .stripIndent()
    if (!workflow.success) {
        subject = """\
          [brucemoran/MANEline] FAILURE [$workflow.runName]
          """
          .stripIndent()
    }

    def msg = """\
      Pipeline execution summary
      ---------------------------
      RunName     : ${workflow.runName}
      Completed at: ${workflow.complete}
      Duration    : ${workflow.duration}
      workDir     : ${workflow.workDir}
      exit status : ${workflow.exitStatus}
      """
      .stripIndent()

    sendMail(to: "${params.email}",
             subject: subject,
             body: msg)
  }
}
