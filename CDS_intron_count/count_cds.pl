#!/usr/bin/perl
##by Li Lei, 201907011,Walnut Creek;
#This is to count the proportion of different LRT from LTR_retriever;
#usage: 

use strict;
use warnings;
use Data::Dumper;
my $file= $ARGV[0];

my %tehash;


open(SNPID,  "$file") or die "Could not open $file";
#my $header = <SNPID>;
#print "$header";
my $sum=0;
foreach my $row (<SNPID>){
        chomp $row;
        my @rtemp= split(/\t/,$row);
        my $delta = abs($rtemp[4]-$rtemp[3]);
        $sum = $sum + $delta; 
}
print "Total length for the sequence: $sum bp\n";
close (SNPID);

