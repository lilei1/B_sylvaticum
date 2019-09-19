#!/usr/bin/perl
##by Li Lei, 20190829,Walnut Creek;
#this is to extract the transcript with .1 to do the GO analysis
#usage: 
use strict;
use warnings;
use Data::Dumper;
my $file = $ARGV[0];

my %hash;


open(GID,  "$file") or die "Could not open $file";
my $header = <GID>;
print "$header";
foreach my $row (<GID>){
        chomp $row;
        my @rtemp = split(/\t/,$row);
        my @tt = split(/\./,$rtemp[2]);
        if ($tt[1]==1){
            print "$row\n";
        }
}
close (GID);

