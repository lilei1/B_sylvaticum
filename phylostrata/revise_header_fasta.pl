#!/usr/bin/perl
my $file = shift;
###20210115 by Li Lei, El Cerrito, CA
###aim: This is for replace the header of each fasta sequence with the taxonomy
#This script is download from https://github.com/MorrellLAB/Barley_Inversions/blob/master/analyses/BAD_Mutations/script/fasta_splitter.pl
open (INFILE, "< $file")or die "Can't open $file";
while (<INFILE>) {
		$line = $_;
		chomp $line;
		if ($line =~ /\>/) { #if has fasta >
			#close OUTFILE;
			my $seq_name = substr($line,1);
			my @array = split(/\.p/, $seq_name);
            my $gid = $array[0];
            my $new_head = $gid." | [Brachypodium sylvaticum] | [Eukaryota; Viridiplantae; Streptophyta; Streptophytina; Embryophyta; Tracheophyta; Euphyllophyta; Spermatophyta; Magnoliopsida; Mesangiospermae; Liliopsida; Petrosaviidae; commelinids; Poales; Poaceae; BOP clade; Pooideae; Brachypodieae; Brachypodium]";
            #print "$new_head\n";
            $line =~ s/$line/\>$new_head/g; 
			#my @tem = split (/\s+/,$array[0]);
			#my @rtem = split (/\./,$tem[0]);
			#my $new_file = $rtem[0]."_".$rtem[1];
			   #$new_file =~ s/^\s+|\s+$//g; #remove all the whitespace;
			#print "$new_file\n";
			#$new_file .= ".fasta";
			#print "$new_file\n";
			#open (OUTFILE, ">$new_file")or die "Can't open: $new_file $!";
		}
		print "$line\n";
}
close INFILE;

#| [Brachypodium sylvaticum] | [Eukaryota; Viridiplantae; Streptophyta; Streptophytina; Embryophyta; Tracheophyta; Euphyllophyta; Spermatophyta; Magnoliopsida; Mesangiospermae; Liliopsida; Petrosaviidae; commelinids; Poales; Poaceae; BOP clade; Pooideae; Brachypodieae; Brachypodium]