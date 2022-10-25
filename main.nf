#!/usr/bin/env nextflow

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

  Optional arguments:

    --bedFile       [str]     Path to BED file of regions to annotate with
                              MANE. Default: null.

    --bedAssembly   [str]     One of GRCh37, GRCh38 indicating the assembly used in --bedFile.
                              Default: null.

    --bedOtherAss   [bool]    Return the bedFile with the other assembly lifted over? Default: true.
    """.stripIndent()
}

if (params.help) exit 0, helpMessage()
//require runID
if(params.runID == null){
  exit 1, "Please specify --runID"
}

//Get GRCh38 (only assembly available) MANE GTF and summary files based on input versions required
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

//map version
vers_get.map { it.text.strip() }.set{ vers_mane }

//parse the GTF to BED
process GtfBed {

  label 'process_low'
  publishDir "${params.outDir}/bed", mode: "copy"

  input:
  tuple file(gtf_gz), file(sum_gtf) from bed_gtf
  val(vers) from vers_mane

  output:
  file("GRCh38.MANE.${vers}.exon.bed") into ( feat_bed, just_feat_bed, sendmail_mane )
  val(vers) into vers_mane_1

  script:
  """
  perl ${workflow.projectDir}/assets/pexon.pl ${gtf_gz} GRCh38.MANE.${vers}.exon.bed
  """
}

//Liftover MANE to GRCh37
process Liftover {

  label 'process_low'
  publishDir "${params.outDir}/bed", mode: "copy"

  input:
  file(exon_bed) from feat_bed
  val(vers) from vers_mane_1

  output:
  file("GRCh37.MANE.${vers}.exon.bed") into lift_feat_bed
  val(vers) into vers_mane_2

  script:
  """
  ##lift
  echo "nameserver 8.8.8.8" > /tmp/resolv.conf
  wget http://hgdownload.cse.ucsc.edu/goldenPath/hg38/liftOver/hg38ToHg19.over.chain.gz
  liftOver ${exon_bed} hg38ToHg19.over.chain.gz GRCh37.MANE.${vers}.exon.bed unmapped
  """
}

//if a BED is supplied overlapp the created Liftover
if( params.bedFile != null ){
  process Overlap37 {
    label 'process_low'
    publishDir "${params.outDir}/bed", mode: "copy"

    input:
    file(exon_lift) from lift_feat_bed
    file(bed_user) from Channel.fromPath( "${params.bedFile}" )
    val(vers) from vers_mane_2

    output:
    tuple file("GRCh37.MANE.${vers}.exon.overlap.${bed_user}"), file(bed_user), file(exon_lift) into sendmail_over
    val(vers) into vers_mane_3

    script:
    """
    ##if bed_user overlaps, return all exons for the gene overlapped
    perl ${workflow.projectDir}/assets/pover.pl ${exon_lift} ${bed_user}  GRCh37.MANE.${vers}.exon.overlap.${bed_user}
    """
  }
  if(params.bedOtherAss){
    process OtherAss37 {
      label 'process_low'
      publishDir "${params.outDir}/bed", mode: "copy"

      input:
      tuple file(bed_ass), file(bed_in), file(exon_lift) from sendmail_over
      val(vers) from vers_mane_3

      output:
      tuple file(bed_ass), file(bed_in), file(exon_lift), file("GRCh38*.${bed_in}"), file('SVUHMolReport_exns_postgresql.csv') into sendmail_asss

      script:
      def bedname = "${bed_ass}".replace('GRCh37','GRCh38')
      """
      ##count total fields as liftOver needs to know this
      LC=\$(head -n5 ${bed_ass} | tail -n1 | awk -F"\\t" '{ print NF-2 }')
      echo "nameserver 8.8.8.8" > /tmp/resolv.conf
      wget http://hgdownload.cse.ucsc.edu/goldenPath/hg19/liftOver/hg19ToHg38.over.chain.gz
      liftOver -bedPlus=\$LC ${bed_ass} hg19ToHg38.over.chain.gz ${bedname} unmapped

      ##parse into table format desired
      perl ${workflow.projectDir}/assets/ptabl.pl ${bed_ass} ${bedname} SVUHMolReport_exns_postgresql.csv
      """
    }
    sendmail_asss.mix(sendmail_mane).set { sendmail_beds }
  } else {
    sendmail_over.mix(sendmail_mane).set { sendmail_beds }
  }
}


// ZIP for sending on sendmail
process zipup {

  label 'low_mem'
  publishDir "${params.outDir}/bed", mode: 'copy'

  input:
  file(send_beds) from sendmail_beds.collect()

  output:
  file("${params.runID}.MANEline.zip") into send_zip

  script:
  """
  zip -r ${params.runID}.MANEline.zip *
  LSL=\$(ls -l ${params.runID}.MANEline.zip | cut -d" " -f5)
  ##if zip too big
  if [[ \$LSL > 6000000 ]]; then
    rm ${params.runID}.MANEline.zip
    zip -r ${params.runID}.MANEline.zip * -x GRCh38*
  fi
  LSL=\$(ls -l ${params.runID}.MANEline.zip | cut -d" " -f5)
  ##if zip too big
  if [[ \$LSL > 6000000 ]]; then
    SWRID=\$(echo ${params.runID} | cut -c -6)
    rm ${params.runID}.MANEline.zip
    zip -r ${params.runID}.MANEline.zip \$SWRID*
  fi
  """
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

    def attachments = send_zip.toList().getVal()

    sendMail(to: "${params.email}",
             subject: subject,
             body: msg,
             attach: attachments)
  }
}
