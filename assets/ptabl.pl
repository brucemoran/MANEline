#! perl

use strict;
use warnings;

# chr VARCHAR ( 100 ) NOT NULL,
# start_37 VARCHAR ( 100 ) NOT NULL,
# end_37 VARCHAR ( 100 ) NOT NULL,
# start_38 VARCHAR ( 100 ) NOT NULL,
# end_38 VARCHAR ( 100 ) NOT NULL,
# strand VARCHAR ( 100 ) NOT NULL,
# gene VARCHAR ( 100 ) NOT NULL,
# ensg VARCHAR ( 100 ) NOT NULL,
# enst VARCHAR ( 100 ) NOT NULL,
# rfid VARCHAR ( 100 ) NOT NULL,
# ensp VARCHAR ( 100 ) NOT NULL,
# ense VARCHAR ( 100 ) NOT NULL,
# exon VARCHAR ( 100 ) NOT NULL

open(B37, $ARGV[0]);
my %b37h;
while(<B37>){
  chomp;
  my @sp=split("\t");
  push(@{$b37h{$sp[3]}},"$sp[0]\t$sp[1]\t$sp[2]");
}
close B37;

open(B38, $ARGV[1]);
open(OUT, ">$ARGV[2]");
while(<B38>){
  chomp;
  my ($chr,$str,$end,$anno)=split(/\t/, $_, 4);
  if(exists $b37h{$anno}){
    for my $val(@{$b37h{$anno}}){
      $anno=~s/;/,/g;
      $val=~s/\t/,/g;
      print OUT "$val,$str,$end,$anno\n";
    }
  }
}
close B38;
close OUT;
