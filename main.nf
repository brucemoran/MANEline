#!/usr/bin/env nextflow

def helpMessage() {
  log.info"""
  ------------------------------------------------------------------------------
                            MANEline PIPELINE
  ------------------------------------------------------------------------------
  Usage:

  nextflow run brucemoran/MANEline

  Description:

  Convert MANE GRCh38 annotation to BED format, and lift backwards to GRCH37.
  Optionally supply a BED file for which annotations are parsed and included.

  Mandatory arguments:

    -profile        [str]     Configuration profile (required: singularity)

    --runID         [str]     Identifier for the run

    --mane_version  [str]     Numeric version of MANE to use (default:
                              current; please only supply numeric e.g. 0.91,
                              not release_0.91)

    --feature       [str]     One of 'exon' (ENST, ENSE, exon_no.), 'transcript' (gene, ENSG, ENST, ENSP, NMID) or 'txps_exon' (default: transcript, then ENSE, exon_no.)

    --email         [str]     Email address to send reports

  Optional arguments:

    --bedFile       [str]     Path to BED file of regions to annotate with
                              MANE (ENS, REF transcript ID, gene name). Default: null.

    --bedAssembly   [str]     One of GRCh37, GRCh38 indicating the assembly used in --bedFile.
                              Default: null.

    --bedOtherAss [bool]    Return the bedFile with the other assembly lifted over? Default: true.
    """.stripIndent()
}

if (params.help) exit 0, helpMessage()
//require runID
if(params.runID == null){
  exit 1, "Please specify --runID"
}
//require correct feature selection
if(params.feature != "txps_exon" && params.feature != "transcript" && params.feature != "exon"){
  exit 1, "Please use one of --feature exon | transcript | txps_exon"
}

//Get GRCh38 (only assembly available) MANE GTF and summary files based on input versions required
process Download {

  label 'process_low'
  publishDir "${params.outDir}/${params.runID}/download", mode: "copy"

  output:
  tuple file("*.ensembl_genomic.gtf.gz"), file("*.summary.txt.gz") into ( bed_gtf, liftover )
  file('vers.txt') into vers_get
  file('feat.txt') into feat_get

  script:
  def mane_base = params.mane_base
  def mane_vers = params.mane_version == "current" ? "current" : "release_" + params.mane_version
  """
  VERS=\$(curl -l ftp.ncbi.nlm.nih.gov/refseq/MANE/MANE_human/${mane_vers}/ | grep summary | cut -d "." -f3,4)
  wget ${mane_base}/${mane_vers}/MANE.GRCh38.\$VERS.ensembl_genomic.gtf.gz
  wget ${mane_base}/${mane_vers}/MANE.GRCh38.\$VERS.summary.txt.gz
  echo \$VERS > vers.txt
  echo ${params.feature} > feat.txt
  """
}

//map version and feature
vers_get.map { it.text.strip() }.set{ vers_mane }
feat_get.map { it.text.strip() }.set{ feat_mane }

//parse the GTF to BED
process GtfBed {

  label 'process_low'
  publishDir "${params.outDir}/${params.runID}/bed", mode: "copy"

  input:
  tuple file(gtf_gz), file(sum_gtf) from bed_gtf
  val(vers) from vers_mane
  val(feat) from feat_mane

  output:
  file("GRCh38.MANE.${vers}.${feat}.bed") into ( feat_bed, just_feat_bed )
  val(vers) into vers_mane_1
  val(feat) into feat_mane_1

  script:
  """
  gunzip -c ${gtf_gz} | sed 's/\\"//g' | sed 's/;//g' | \\
    if [[ ${feat} == "transcript" ]]; then
      perl -ane 'chomp; if(\$F[2] eq "transcript"){
      print "\$F[0]\\t\$F[3]\\t\$F[4]\\t\$F[6];\$F[15];\$F[9];\$F[11];\$F[23];\$F[25]\\n";}'
    elif [[ ${feat} == "exon" ]]; then
      perl -ane 'chomp; if(\$F[2] eq "exon"){
      print "\$F[0]\\t\$F[3]\\t\$F[4]\\t\$F[6];\$F[11];\$F[23];exon_\$F[21]\\n";}'
    else
      perl -ane 'chomp; if(\$F[2] eq "exon"){
      print "\$F[0]\\t\$F[3]\\t\$F[4]\\t\$F[6];\$F[15];\$F[9];\$F[11];\$F[27];\$F[29];\$F[23];exon_\$F[21]\\n";}'
    fi | \\
  sed 's/RefSeq://g' >> GRCh38.MANE.${vers}.${feat}.bed
  """
}

//Liftover MANE to GRCh37
process Liftover {

  label 'process_low'
  publishDir "${params.outDir}/${params.runID}/bed", mode: "copy"

  input:
  file(bed_feat) from feat_bed
  val(vers) from vers_mane_1
  val(feat) from feat_mane_1

  output:
  file("GRCh37.MANE.${vers}.${feat}.bed") into lift_feat_bed
  val(vers) into vers_mane_2
  val(feat) into feat_mane_2

  script:
  """
  ##lift
  echo "nameserver 8.8.8.8" > /tmp/resolv.conf
  wget http://hgdownload.cse.ucsc.edu/goldenPath/hg38/liftOver/hg38ToHg19.over.chain.gz
  liftOver ${bed_feat} hg38ToHg19.over.chain.gz GRCh37.MANE.${vers}.${feat}.bed unmapped
  """
}

//if a BED is supplied overlapp the created Liftover
if( params.bedFile != null ){
  if(params.bedAssembly == "GRCh38" ){
    process Overlap38 {
      label 'process_low'
      publishDir "${params.outDir}/${params.runID}/bed", mode: "copy"

      input:
      file(feat_just) from just_feat_bed
      file(bed_over) from Channel.fromPath( "${params.bedFile}" )
      val(vers) from vers_mane_2
      val(feat) from feat_mane_2

      output:
      tuple file("*.overlap.MANE.${vers}.${feat}.bed"), file(bed_over), file(feat_just) into sendmail_over

      script:
      def bedname = "${bed_over}".split('\\.')[0]
      """
      ##overlap
      perl ${workflow.projectDir}/assets/pover.pl ${feat_just} ${bed_over} 1
      uniq 1 > ${bedname}.GRCh38.overlap.MANE.${vers}.${feat}.bed
      rm 1
      """
    }
    if(params.bedOtherAss){
      process OtherAss38 {
        label 'process_low'
        publishDir "${params.outDir}/${params.runID}/bed", mode: "copy"

        input:
        tuple file(bed_ass), file(bed_in), file(bed_just) from sendmail_over
        val(vers) from vers_mane_3
        val(feat) from feat_mane_3

        output:
        tuple file(bed_ass), file(bed_in), file("*.overlap.MANE.${vers}.${feat}.bed") into sendmail_ass

        script:
        def bedname = "${bed_ass}".split('\\.')[0]
        """
        echo "nameserver 8.8.8.8" > /tmp/resolv.conf
        wget http://hgdownload.cse.ucsc.edu/goldenPath/hg38/liftOver/hg38ToHg19.over.chain.gz
        liftOver ${bed_ass} hg38ToHg19.over.chain.gz ${bedname}.GRCh37.overlap.MANE.${vers}.${feat}.bed unmapped
        """
      }
      sendmail_ass.set { sendmail_bed }
    } else {
      sendmail_over.set { sendmail_bed }
    }
  } else {
    process Overlap37 {
      label 'process_low'
      publishDir "${params.outDir}/${params.runID}/bed", mode: "copy"

      input:
      file(feat_lift) from lift_feat_bed
      file(bed_over) from Channel.fromPath( "${params.bedFile}" )
      val(vers) from vers_mane_2
      val(feat) from feat_mane_2

      output:
      tuple file("*.overlap.MANE.${vers}.${feat}.bed"), file(bed_over), file(feat_lift) into sendmail_over
      val(vers) into vers_mane_3
      val(feat) into feat_mane_3

      script:
      def bedname = "${bed_over}".split('\\.')[0]
      """
      ##overlap
      perl ${workflow.projectDir}/assets/pover.pl ${feat_lift} ${bed_over} 1
      uniq 1 > ${bedname}.GRCh37.overlap.MANE.${vers}.${feat}.bed
      rm 1
      """
    }
    if(params.bedOtherAss){
      process OtherAss37 {
        label 'process_low'
        publishDir "${params.outDir}/${params.runID}/bed", mode: "copy"

        input:
        tuple file(bed_ass), file(bed_in), file(bed_lift) from sendmail_over
        val(vers) from vers_mane_3
        val(feat) from feat_mane_3

        output:
        tuple file(bed_ass), file(bed_in), file(bed_lift), file("*.overlap.MANE.${vers}.${feat}.bed") into sendmail_asss

        script:
        def bedname = "${bed_over}".split('\\.')[0]
        """
        echo "nameserver 8.8.8.8" > /tmp/resolv.conf
        wget http://hgdownload.cse.ucsc.edu/goldenPath/hg19/liftOver/hg19ToHg38.over.chain.gz
        liftOver ${bed_ass} hg19ToHg38.over.chain.gz ${bedname}.GRCh38.overlap.MANE.${vers}.${feat}.bed unmapped
        """
      }
     sendmail_asss.set { sendmail_beds }
    } else {
     sendmail_over.set { sendmail_beds }
    }
  }
}

// ZIP for sending on sendmail
process zipup {

  label 'low_mem'
  publishDir "${params.outDir}/${params.runID}/bed", mode: 'copy'

  input:
  file(send_beds) from sendmail_beds.collect()

  output:
  file("${params.runID}.MANEline.zip") into send_zip

  script:
  """
  zip -r $${params.runID}.MANEline.zip *
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
