#! bash

##convert Oncomine fusion to bed
##input is the Oncomine file Oncomine_Focus_designs_063020_Reference.bed
##it is not a bed
##it is a big mess
##we want to be able to report format:
##GENE1::GENE2 NM_TRANSC1:r.123_456::NM_TRANSC2:r.321_654
##for the r. we need the transcript info, we have breakpoint start-end
##in the file as FP/TP_START,END but this is genomic coordinates...

##read the data
BED=$1;

tail -n+2 $BED |\
  perl -ane 'chomp; if($F[5]=~m/^TYPE/){
    @s=split(/\;/, $F[5]);
    for $x(@s){
      if($x=~m/^FP_GENE_ID/){$x=~s/FP_GENE_ID=//g; $fp_g = $x}
      if($x=~m/^FP_CHROM/){$x=~s/FP_CHROM=//g; $fp_c = $x}
      if($x=~m/^FP_START/){$x=~s/FP_START=//g; $fp_s = $x}
      if($x=~m/^FP_END/){$x=~s/FP_END=//g; $fp_e = $x}
      if($x=~m/^TP_GENE_ID/){$x=~s/TP_GENE_ID=//g; $tp_g = $x}
      if($x=~m/^TP_CHROM/){$x=~s/TP_CHROM=//g; $tp_c = $x}
      if($x=~m/^TP_START/){$x=~s/TP_START=//g; $tp_s = $x}
      if($x=~m/^TP_END/){$x=~s/TP_END=//g; $tp_e = $x}
    }
    ##split on comma as some start end are multiples
    if($fp_s=~m/\,/){
      @sfs=split(/\,/, $fp_s);
      @sfe=split(/\,/, $fp_e);
      chomp @sfs; chomp @sfe;
      for($i=0;$i<@sfs;$i++){
        print "$fp_c\t$sfs[$i]\t$sfe[$i]\t$fp_g;$F[3]\n";
      }
    } else {
      print "$fp_c\t$fp_s\t$fp_e\t$fp_g;$F[3]\n";
    }
    if($tp_s=~m/\,/){
      @sts=split(/\,/, $tp_s);
      @ste=split(/\,/, $tp_e);
      chomp @sts; chomp @ste;
      for($i=0;$i<@sts;$i++){
        print "$tp_c\t$sts[$i]\t$ste[$i]\t$tp_g;$F[3]\n";
      }
    } else {
      print "$tp_c\t$tp_s\t$tp_e\t$tp_g;$F[3]\n";
    }}' > $2
