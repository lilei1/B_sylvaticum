#!/usr/bin/perl
##by Li Lei, 20190620,Walnut Creek;
#this is to masked all of the number. Treat all of the non-zero as 1. This file is for making vann diagram
#usage: 
use strict;
use warnings;
use Data::Dumper;
my $file = $ARGV[0];

my %gidhash;


open(SNPID,  "$file") or die "Could not open $file";
#my $header = <SNPID>;
#print "$header";
foreach my $row (<SNPID>){
        chomp $row;
        my @rtemp= split(/\t/,$row);
        my $key = $rtemp[0]."_".$rtemp[1]."_".$rtemp[2]."_".$rtemp[3]."_".$rtemp[4]."_".$rtemp[5];
        push @{$gidhash{$key}}, $rtemp[6];
}
close (SNPID);

print "Bdistachyon\tBstacei\tBsylvaticum\tOsativaKitaake\tSbicolor\tPhallii\tcount\n";
my $count=0;
foreach my $key (keys %gidhash){
 			my @arrary = @{$gidhash{$key}};
 			foreach my $ele (@arrary){
 				$count = $count + $ele;
 			}
 			my @tt = split(/\_/,$key);
 			print "$tt[0]\t$tt[1]\t$tt[2]\t$tt[3]\t$tt[4]\t$tt[5]\t$count\n";
 			$count=0;
}