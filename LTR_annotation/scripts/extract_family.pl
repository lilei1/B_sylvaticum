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
        my @rtemp= split(/\t/,$row);
        if($rtemp[0] !~ "^Summary:"){
            my $key = $rtemp[10];
            if ($key ne ''){
            push @{$tehash{$key}}, $rtemp[8];
            }
         }
}
close (SNPID);
#print Dumper(\%tehash);
#print "Categories\tProp\n";
my $prop = 0;
foreach my $key (keys %tehash){
 			my @arrary = @{$tehash{$key}};
            my $ele = 0;
 			foreach $ele (@arrary){
 				$prop = $prop + $ele;
 			}
 			#my $prop = ($count/$total)*100;
 			print "$key\t$prop\n";
 			$prop=0;
}
