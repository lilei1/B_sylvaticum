#!/usr/bin/perl
##by Li Lei, 20190829,Walnut Creek;
#this is to count the genes for each GO id;
#usage: 
use strict;
use warnings;
use Data::Dumper;
my $gene = $ARGV[0];

my %hash;


open(GID,  "$gene") or die "Could not open $gene";
my $header = <GID>;
foreach my $row (<GID>){
        chomp $row;
        my @rtemp = split(/\t/,$row);
        my @tt = split(/,/,$rtemp[9]);
        
        if($#tt == 0){
            print "$tt[0]\t$rtemp[2]\n";
        }
        else{
            foreach my $ele (@tt){
                print "$ele\t$rtemp[2]\n";      
            }
        }
}

close (GID);

