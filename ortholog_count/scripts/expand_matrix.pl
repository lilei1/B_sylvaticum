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
my $header = <SNPID>;
print "Sbicolor\tOsativaKitaake\tPhalliiHAL\tBsylvaticum\tBstacei\tBdistachyon\tcount\n";
foreach my $row (<SNPID>){
        chomp $row;
        my @rtemp= split(/\t/,$row);
        my $nb = $rtemp[6] - 1;
         for (my $i=0; $i <= $nb; $i++ ){
         	     print "$rtemp[0]\t$rtemp[1]\t$rtemp[2]\t$rtemp[3]\t$rtemp[4]\t$rtemp[5]\n";

         }
}
close (SNPID);
