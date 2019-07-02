#!/usr/bin/perl
##by Li Lei, 20190702,Walnut Creek;
#this is to count the repetitive length for all of the repetive segments
#usage: 
use strict;
use warnings;
#use Data::Dumper;
my ($file,$total)= @ARGV;

my %tehash;


open(SNPID,  "$file") or die "Could not open $file";
my $header = <SNPID>;
#print "$header";
foreach my $row (<SNPID>){
        chomp $row;
        my @rtemp= split(/\t/,$row);
        my $key = $rtemp[2];
        push @{$tehash{$key}}, $rtemp[1];
}
close (SNPID);

print "Categories\tLen\tProp\n";
my $count=0;
foreach my $key (keys %tehash){
 			my @arrary = @{$tehash{$key}};
 			foreach my $ele (@arrary){
 				$count = $count + $ele;
 			}
 			my $prop = ($count/$total)*100;
 			print "$key\t$count\t$prop\n";
 			$count=0;
}