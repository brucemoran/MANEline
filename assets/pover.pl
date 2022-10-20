#! perl

use strict;
use warnings;

open(BED, $ARGV[0]);
my %regions;
while(<BED>){
  chomp;
  my @sp=split(/\t/, $_, 4);
  push(@{$regions{$sp[0]}},"$sp[1]\t$sp[2]\t$sp[3]");
}
close BED;

open(BED, $ARGV[1]);
my %genes;
while(<BED>){
  chomp;
  my ($chr,$str,$end,$anno)=split(/\t/, $_, 4);
  for my $key(@{$regions{$chr}}){
    my @an=split(/\t/, $key);
    my @ans=split(/;/, $an[2]);
    ##start gt anno start, start must also be lt anno end
    if($str>=$an[0] && $str<=$an[1]){
      push(@{$genes{$ans[1]}}, 1);
      next;
    }
    ##start lt anno start, end must be gt anno start
    if($str<=$an[0] && $end>$an[0]){
      push(@{$genes{$ans[1]}}, 1);
      next;
    }
  }
}
close BED;

open(BED, $ARGV[0]);
open(OUT, ">$ARGV[2]");
while(<BED>){
  chomp;
  my ($chr,$str,$end,$anno)=split(/\t/, $_, 4);
  my @ans=split(";", $anno);
  if(exists $genes{$ans[1]}){
    print OUT $_ . "\n";
  }
}
close OUT;
close BED;

sub uniq {
  my %seen;
  return grep { !$seen{$_}++ } @_;
}
