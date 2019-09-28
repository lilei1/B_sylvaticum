#!/usr/bin/perl
##by Li Lei, 20190918,Walnut Creek;
#this is to match the target gene back to the GO term and split the GO list as target and background sets;
#usage: 
use strict;
use warnings;
use Data::Dumper;
my ($file1,$file2,$result1,$result2)= @ARGV;

my %hash;


open(GID,  "$file1") or die "Could not open $file1";
my $header = <GID>;
foreach my $row (<GID>){
        chomp $row;
        my @rtemp = split(/\t/,$row);
        $hash{$rtemp[0]} = 999;
}
close (GID);

open(ID,  "$file2") or die "Could not open $file2";
open(OUT1,  ">$result1") or die "Could not open $result1";
open(OUT2,  ">$result2") or die "Could not open $result2";
my $header1 = <ID>;
foreach my $row (<ID>){
        chomp $row;
        my @tt = split(/\t/,$row);
        if (exists ($hash{$tt[0]})){
            print OUT1 "$row\n";
        }
        else{
            print OUT2 "$row\n";
        }
}

close (ID);
close (OUT1);
close (OUT2);
