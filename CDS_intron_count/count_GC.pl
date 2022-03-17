#!/usr/bin/perl

#read file in from input line

$infile = $ARGV[0];
open(TXT, "$infile") or die "Could not open $infile";
#read in the DNA string using the fasta subfunction
$DNA = &read_fasta();
$len = length($DNA);
print "\n DNA Length is: $len \n";

my $numG=0;
my $numC=0;
my $numT=0;
my $numA=0;
@bases=split(//,$DNA);

foreach $bp(@bases){
    if($bp =~ m/G/i){$numG++};
    if($bp =~ m/C/i){$numC++};
    if($bp =~ m/T/i){$numT++};
    if($bp =~ m/A/i){$numA++};
}
print "\n Number of G bases: $numG";
print "\n Number of C bases: $numC";
print "\n Number of T bases: $numT";
print "\n Number of A bases: $numA";

$GC_content = (($numG+$numC)/$len)*100;
print "\n\n GC Content is: $GC_content % \n";

close(TXT);

sub read_fasta{
    my $sequence = "";
    while(<TXT>)
    {   
        $line = $_;
        #print " $line \n";
       #remove newline characters 
       chomp($line);
       # discard fasta header line
        if($line =~ />/){ next }
       # append the line to the DNA sequence
        else{$sequence .= $line }
    } 
    return($sequence);
}
