#! perl

use strict;
use warnings;

open(BED, $ARGV[0]);
my %regions;
while(<BED>){
    chomp;
    my @sp=split(/\t/);
    push(@{$regions{"$sp[0],$sp[1],$sp[2]"}},$sp[3]);
}
close BED;

open(BED, $ARGV[1]);
open(OUT, ">$ARGV[2]");
while(<BED>){
  chomp;
  my ($chr,$str,$end,$anno)=split(/\t/);
  for my $key(@{$regions{"$chr,$str,$end"}}){
    print OUT "$chr\t$str\t$end\t$key;$anno\n";
  }
}

close BED;
close OUT;
