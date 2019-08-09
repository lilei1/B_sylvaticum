#!/usr/bin/perl
##by Li Lei, 2019/08/09,Walnut Creek;
#this is to extract the genes with DGE based on the pair comparisons 
#usage: 
use strict;
use warnings;
#use Data::Dumper;
my $file = $ARGV[0];


open(GID,  "$file") or die "Could not open $file";
my $header = <GID>;
foreach my $row (<GID>){
        chomp $row;
        my @rtemp = split(/\t/,$row);
        my @arr = ($rtemp[3],$rtemp[6],$rtemp[9],$rtemp[12],$rtemp[15],$rtemp[18],$rtemp[21],$rtemp[24],$rtemp[27],$rtemp[30],$rtemp[33],$rtemp[36],$rtemp[39],$rtemp[42],$rtemp[45]);
        my $count = 0;
        foreach my $ele (@arr){
                    if ($ele eq "TRUE"){
                        $count++;
                    }
        }
        if ($count >= 1){
        print "$rtemp[0]\n";
        
        }
        
}
close (GID);
