#!/usr/bin/perl
my ($file1, $file2)= @ARGV;
###20210115 by Li Lei, El Cerrito, CA
###aim: This is for replace the header of each fasta sequence with the taxonomy
#This script is download from https://github.com/MorrellLAB/Barley_Inversions/blob/master/analyses/BAD_Mutations/script/fasta_splitter.pl

my %hash;
open (IN, "< $file1")or die "Can't open $file1";
my $head = <IN>;
while (<IN>) {
		$line = $_;
		chomp $line;
			#close OUTFILE;
            #print "$line\n";
			#my $seq_name = substr($line,1);
			my @array = split(/\t/, $line);
                $hash{$array[1]} = 999;  
                #print "$array[1]\n";
}
close IN;
#print Dumper(\%hash);


open (INFILE, "< $file2")or die "Can't open $file2";
while (<INFILE>) {
		$line = $_;
		chomp $line;
		if ($line =~ /\>/) { #if has fasta >
			#close OUTFILE;
			my $seq_name = substr($line,1);
			my @array = split(/\|/, $seq_name);
            $array[0] =~ s/\s+$//;
            my $gid = $array[0];
                #print "$gid\n";
            if (exists $hash{$gid}){
                #print "$line\n";
                		next;
            }
            else{
                #print "$line\n";
                print "$gid\n"
            }
		}
        
}
close INFILE;

#| [Brachypodium sylvaticum] | [Eukaryota; Viridiplantae; Streptophyta; Streptophytina; Embryophyta; Tracheophyta; Euphyllophyta; Spermatophyta; Magnoliopsida; Mesangiospermae; Liliopsida; Petrosaviidae; commelinids; Poales; Poaceae; BOP clade; Pooideae; Brachypodieae; Brachypodium]