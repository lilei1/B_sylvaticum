#!/usr/bin/perl
##by Li Lei, 201907010,Walnut Creek;
#this is to count the proportion of different LRT from LTR_retriever;
#usage: 
use strict;
use warnings;
use Data::Dumper;
my $file= $ARGV[0];

my %tehash;


open(SNPID,  "$file") or die "Could not open $file";
my $header = <SNPID>;
#print "$header";
foreach my $row (<SNPID>){
        chomp $row;
        my @rtemp = split(/\t/,$row);
        my @tmp = split(/;/,$rtemp[8]);
        #print "$tmp[0]\n";
        my @tt =split(/\./,$tmp[0]);
        #print "$tt[1]\n";
        if ($tt[1] == 1){
           print "$row\n";
        }
}
close (SNPID);
