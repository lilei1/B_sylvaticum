#!/usr/bin/perl
##by Li Lei, 20200831, in El Cerrito, CA;
#this is to creat a matrix to display the GO term enrichment for each module defined by WGCNA;
#
#usage: 
use strict;
use warnings;

my ($color,$GO_id) = @ARGV;

my %hash; #define a hash to store the GO and the foldchanges for each GO;
open(OUT,  "$color") or die "Could not open $color";
my $header = <OUT>;
foreach my $row (<OUT>){
        chomp $row;
        my @rtemp1 = split(/\t/,$row);
           #print "$rtemp1[0]\n";
           $hash{$rtemp1[0]} = $rtemp1[1]; #key is the GO id and value is the foldchange 
}
close (OUT);

open(OUT,  "$GO_id") or die "Could not open $GO_id";
foreach my $row (<OUT>){
        chomp $row;
        my @rtemp = split(/\t/,$row);
        if (exists $hash{$rtemp[0]}){
        	print "$row\t$hash{$rtemp[0]}\n";
        }
        else{
        	print "$row\tNA\n";
        }

}
close (OUT);