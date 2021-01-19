#!/usr/bin/perl
my ($file1, $file2)= @ARGV;
###20210115 by Li Lei, El Cerrito, CA
###aim: This is for replace the header of each fasta sequence with the taxonomy
#This script is download from https://github.com/MorrellLAB/Barley_Inversions/blob/master/analyses/BAD_Mutations/script/fasta_splitter.pl
use Data::Dumper;


my %hash;
open (IN, "< $file1")or die "Can't open $file1";
while (<IN>) {
		$line = $_;
		chomp $line;
			#close OUTFILE;
			#my $seq_name = substr($line,1);
			my @array = split(/\|/, $line);
            #print "$array[1]\t$array[2]\n";
            $array[1] =~ s/^\s+//;
            $array[1] =~ s/\s+$//;
            $array[2] =~ s/^\s+//;
            $array[2] =~ s/\s+$//;
            my $key = $array[1];
            my $value = $array[2];
            if ($key ne "" and $value ne ""){
               $hash{$key} = $value;   
            }
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
			my @array = split(/\s+/, $seq_name);
            my $gid = $array[0];
            my @matches = $seq_name =~ /\[ ( [^\]]* )\]/xg;#capture the text in the "[]"
               $matches[0]=~ s/^\s+//;
               $matches[0]=~ s/\s+$//;
            my $tax = $matches[0];
            if (exists $hash{$tax}){
            my $new_head = $gid." | "."[".$tax."]"." | "."[".$hash{$tax}."]";
            $line = ">".$new_head
            #print "$line\n";
            #print "$seq_name\n";
            }
        }
        print "$line\n";
 }
close INFILE;


                
            #print "$key\t$value\n";
            #my $new_head = $gid." | [Brachypodium sylvaticum] | [Eukaryota; Viridiplantae; Streptophyta; Streptophytina; Embryophyta; Tracheophyta; Euphyllophyta; Spermatophyta; Magnoliopsida; Mesangiospermae; Liliopsida; Petrosaviidae; commelinids; Poales; Poaceae; BOP clade; Pooideae; Brachypodieae; Brachypodium]";
            #print "$new_head\n";
            #$line =~ s/$line/\>$new_head/g; 
			#my @tem = split (/\s+/,$array[0]);
			#my @rtem = split (/\./,$tem[0]);
			#my $new_file = $rtem[0]."_".$rtem[1];
			   #$new_file =~ s/^\s+|\s+$//g; #remove all the whitespace;
			#print "$new_file\n";
			#$new_file .= ".fasta";
			#print "$new_file\n";
			#open (OUTFILE, ">$new_file")or die "Can't open: $new_file $!";
		#print "$line\n";


#| [Brachypodium sylvaticum] | [Eukaryota; Viridiplantae; Streptophyta; Streptophytina; Embryophyta; Tracheophyta; Euphyllophyta; Spermatophyta; Magnoliopsida; Mesangiospermae; Liliopsida; Petrosaviidae; commelinids; Poales; Poaceae; BOP clade; Pooideae; Brachypodieae; Brachypodium]