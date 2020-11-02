#!/usr/bin/perl
##by Li Lei, 20200915, in El Cerrito, CA;
#this is to match sylvaticum gene ID with distachyon gene ID
#usage: 
use strict;
use warnings;

my ($anno,$matrix) = @ARGV;

my %hash; #define a hash to store the GO and the foldchanges for each GO;
open(F1,  "$anno") or die "Could not open $anno";
my $header = <F1>;
foreach my $row (<F1>){
        chomp $row;
        my @rtemp1 = split(/\t/,$row);
           #print "$rtemp1[3]\t$rtemp1[1]\n";
           $hash{$rtemp1[3]} = $rtemp1[1]; #key is the syl gene id and value is distachyon gene id.
}
close (F1);

open(F2,  "$matrix") or die "Could not open $matrix";
my $header1 = <F2>;
#print "$header1";
foreach my $row (<F2>){
        chomp $row;
        my @rtemp = split(/\t/,$row);
        if (exists $hash{$rtemp[0]}){
        	print "$hash{$rtemp[0]}\n";
        }
        else{
        	#print "NA\n";
        }

}
close (F2);