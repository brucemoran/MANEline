#! perl

##annotate with HGVS from FA into BED
use strict;
use warnings;

open(BED, $ARGV[0]);
my %txps;
while(<BED>){
  chomp;
  if($_!~m/\t/){next;}
  my @sp = split(/\s+/, $_, 4);
  my @an = split(/[;\s]/, $sp[3]);
  for my $k(@an){
    if($k=~m/^ENST/){
      push(@{$txps{$k}},"$sp[0]\t$sp[1]\t$sp[2]\t" . join("\t", @an));
    }
  }
}
close BED;
print "User bed read...\n";

##BED from GTF from MANE to allow finding variants as this has the
open(MNBD, $ARGV[1]);
my %mpos;
while(<MNBD>){
  chomp;
  if($_!~m/\t/){next;}
  my @sp = split(/[;\s]/, $_, 4);
  my @an = split(/;/, $sp[3]);
  if(exists $txps{$an[3]}){
    ##appends the genomic location of the transcript
    for my $bd(@{$txps{$an[3]}}){
      push(@{$txps{$an[3]}}, $bd . "$sp[0]\t$sp[1]\t$sp[2]");
      print $bd . "\t$sp[0]:$sp[1]-$sp[2]\n";
    }
  }
}
close MNBD;
print "MANE bed read...\n";

open(FAGZ, "gunzip -c $ARGV[2]|");
#open(OUT, ">$ARGV[3]");
local $/ = ">";
while(<FAGZ>){
  my($id_line, $sq) = split("\n", $_, 2);
  my @id = split(" ", $id_line);
  my $txp = substr($id_line, 0, 17);
  if(exists $txps{$txp}){
    print "$txp exists\n";
    $sq=~s/\n//g;
    for my $bd(@{$txps{$txp}}){
      my ($chr, $str, $end, $ann) = split("\t", $bd, 4);
      my @an = split("\t", $ann);
      my @txloc = split("[:-]", $an[-1]);
      print "@txloc exists\n";
      ##split sequence into individual characters, then return
      ##before the var, the var (end), and after the var
      ##count into sq the number of fstrt-$end (variant is $end)
      my @sqs = split("", $sq);
      my $varsq = $end-$txloc[1];
      print $end . " - " . $txloc[1] . " = " . $varsq . "\n";
      print $sqs[$varsq-1] . $sqs[$varsq] . $sqs[$varsq+1] . " -> " . $an[3] . "\n";
    }
  }
}

#   for my $key(@{$txps{$chr}}){
#     my @an=split(/\t/, $key);
#     ##start gt anno start, start must also be lt anno end
#     if($str>=$an[0] && $str<=$an[1]){
#         print OUT "$chr\t$str\t$end\t$anno\t$an[2]\n";
#         next;
#     }
#     ##start lt anno start, end must be gt anno start
#     if($str<=$an[0] && $end>$an[0]){
#         print OUT "$chr\t$str\t$end\t$anno\t$an[2]\n";
#         next;
#     }
#   }
# }
close FAGZ;
close OUT;
