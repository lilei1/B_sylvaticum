#!/usr/bin/perl
##by Li Lei, 20190830,Walnut Creek;
#this is to count the GO id share by the six species:
#usage: 
use strict;
use warnings;
use Data::Dumper;
my @(file1,file2)= @ARGV;

my %hash;


open(GID,  "$file1") or die "Could not open $file1";
my $header = <GID>;
foreach my $row (<GID>){
        chomp $row;
        my @rtemp = split(/\t/,$row);
        $hash{$rtemp[0]} = $rtemp[1];
}

close (GID);

open(ID,  "$file2") or die "Could not open $file2";
my $header = <ID>;
foreach my $row (<ID>){
        chomp $row;
        my @tt = split(/\t/,$row);
        if (exists ($hash{$tt[0]})){
            print "$row\t$hash{$tt[0]}\n";
        }
}

close (ID);
