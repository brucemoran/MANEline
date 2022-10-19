#! perl

use strict;
use warnings;

open(BED, $ARGV[0]);
my %regions;
while(<BED>){
    chomp;
    my @sp=split(/\t/, 4);
    push(@{$regions{$sp[0]}},"$sp[1]\t$sp[2]\t$sp[3]");
}
close BED;

open(BED, $ARGV[1]);
open(OUT, ">$ARGV[2]");
while(<BED>){
    chomp;
    my ($chr,$str,$end, $anno)=split(/\t/, $_, 4);
    for my $key(@{$regions{$chr}}){
        my @an=split(/\t/, $key);
        ##start gt anno start, start must also be lt anno end
        if($str>=$an[0] && $str<=$an[1]){
            print OUT "$chr\t$str\t$end\t$anno\t$an[2]\n";
            next;
        }
        ##start lt anno start, end must be gt anno start
        if($str<=$an[0] && $end>$an[0]){
            print OUT "$chr\t$str\t$end\t$anno\t$an[2]\n";
            next;
        }
    }
}
close BED;
close OUT;
