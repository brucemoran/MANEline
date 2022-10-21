#! perl

use strict;
use warnings;

open(GTF, "gunzip -c $ARGV[0] |");
open(OUT, ">$ARGV[1]");

while(<GTF>){
  chomp;
  my @sp=split("\t", $_, 9);
  if($sp[2] eq "exon"){
    my $out="$sp[0]\t$sp[3]\t$sp[4]\t$sp[6];";
    my @an=split(";", $sp[8]);
    my %ann;
    foreach my $k(@an){
      chomp $k;
      $k =~s/^\s+//;
      if($k =~m/gene_id/){$ann{"gene_id"}=parz($k);}
      if($k =~m/transcript_id/){$ann{"transcript_id"}=parz($k);}
      if($k =~m/gene_name/){$ann{"gene_name"}=parz($k);}
      if($k =~m/exon_number/){$ann{"exon_number"}="exon_" . parz($k);}
      if($k =~m/exon_id/){$ann{"exon_id"}=parz($k);}
      if($k =~m/protein_id/){$ann{"protein_id"}=parz($k);}
      if($k =~m/db_xref/){
        my $rfs=parz($k); $rfs=~s/RefSeq://; $ann{"db_xref"}=$rfs;
      }
    }
    print OUT $out . $ann{"gene_name"} . ";" . $ann{"gene_id"} . ";" . $ann{"transcript_id"} . ";" . $ann{"db_xref"} . ";" . $ann{"protein_id"} . ";" . $ann{"exon_id"} . ";" . $ann{"exon_number"} . "\n";
  }
}
close GTF;
close OUT;

sub parz {
  my @s = split(/ /, $_[0]);
  $s[1]=~s/\"//g;
  return $s[1];
}
