#!/usr/bin/perl
##by Li Lei, 20190819,Walnut Creek;
#this is to count the members in each species for each cluster:
use strict;
use warnings;
use Data::Dumper;
my $file = $ARGV[0];

#initiate a hash:
#my %shash = (
#    "Sbicolor"  => "annual",
#    "Osativa" => "annual",
#    "Bdistachyon"  => "annual",
#    "Bsylvaticum" => "perennial"
#    "Bstacei" => "perennial"
#   "Phallii" => "perennial"
#);
#initiate a array:

my @species = qw(Sbicolor Osativa Bdistachyon Bsylvaticum Bstacei Phallii);

open(SNPID,  "$file") or die "Could not open $file";
#my $header = <SNPID>;
#print "$header";
foreach my $row (<SNPID>){
        chomp $row;
        my @rtemp = split(/\t/,$row);
        my $id = $rtemp[0];
        print "$id\t";
        my @tmp = split(/\,/,$rtemp[1]);
        my @arr;
        foreach my $ele (@tmp){
                my @tt = split(/:/,$ele);
                push @arr, $tt[0];
        }
        for my $mnz (@species) {
            my $zs; 
            $zs = grep { $_ eq $mnz } @arr; #count the occrence of an element in an array
            print "$mnz:$zs\t";
            $zs = 0; #clear the counting
        }
        print "\n";
        @arr = ''; #clear the array
}
close (SNPID);