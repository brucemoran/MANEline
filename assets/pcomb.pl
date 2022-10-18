#! perl

use strict;
use warnings;

open(BED, $ARGV[0]);
my %regions;
while(<BED>){
    chomp;
    my @sp=split(/\t/);
    push(@{$regions{$sp[0]}},"$sp[1]\t$sp[2]\t$sp[3]");
}
close BED;

open(BED, $ARGV[1]);
open(OUT, ">$ARGV[2]");
while(<BED>){
  chomp;
  my ($chr,$str,$end,$anno)=split(/\t/);
  print $_ . "\n";
  for my $key(@{$regions{$chr}}){
    my @an=split(/\t/, $key);
    print OUT "$chr\t$str\t$end\t$an[2];$anno\n";
  }
}

close BED;
close OUT;
