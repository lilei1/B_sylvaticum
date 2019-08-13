#!/usr/bin/perl
##by Li Lei, 20190620,Walnut Creek;
#this is to masked all of the number. Treat all of the non-zero as 1. This file is for making vann diagram
#usage: 
use strict;
use warnings;
#use Data::Dumper;
my $file = $ARGV[0];

my %gidhash;


open(SNPID,  "$file") or die "Could not open $file";
#my $header = <SNPID>;
#print "$header";
foreach my $row (<SNPID>){
        chomp $row;
        my @rtemp= split(/\t/,$row);
        my @tmp = @rtemp[0..5];
        #print "$rtemp[5]\n";
        foreach my $ele (@tmp){
                   if ($ele !=0){
                       print "1\t";
                   }
                   else{
                       print "0\t";
                   }
         #print "$ele\t"
        }
        print "$rtemp[6]\n";
}
close (SNPID);
