#!/usr/bin/perl
##by Li Lei, 20190830,Walnut Creek;
#this is to count the genes for each GO id;
#usage: 
use strict;
use warnings;
use Data::Dumper;
my $gene = $ARGV[0];

my %hash;


open(GID,  "$gene") or die "Could not open $gene";
#my $header = <GID>;
foreach my $row (<GID>){
        chomp $row;
        my @rtemp = split(/\t/,$row);
        push @{$hash{$rtemp[0]}}, $rtemp[1];
}

close (GID);
print "GO_id\tgene_nb\n";

foreach my $key (keys %hash){
 		my @arr = @{$hash{$key}};
        my $len = $#arr + 1;
        print "$key\t$len\n";

}