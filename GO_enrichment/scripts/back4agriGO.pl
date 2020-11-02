#!/usr/bin/perl
##by Li Lei, 20200915, in El Cerrito, CA;
#this is to change the format of the gene-GO for the format uploading for AgriGO; 
#usage: 
use strict;
use warnings;

my $file = $ARGV[0];

#my %hash; #define a hash to store the GO and the foldchanges for each GO;
open(F1,  "$file") or die "Could not open $file";
my $header = <F1>;
foreach my $row (<F1>){
        chomp $row;
        my @rtemp = split(/\t/,$row);
        my @tt = split (/\,/,$rtemp[1]); 
        foreach my $ele (@tt){
            print "$rtemp[0]\t$ele\n";
        }
}
close (F1);
