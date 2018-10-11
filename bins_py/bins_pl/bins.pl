#!/usr/bin/perl
use warnings;
use 5.12.4;
use File::Basename;
#I modified the script little bit and pasted below. If you want bigger chunks,
#just change "my $record_per_file = 4" value to any chunk size you like.
#Split FASTA files into chunks determined by user.
#by RNAeye

my $file = "input.fa"; 		#enter name of your FASTA file here
my $record_per_file = 4;	#Enter how many record you want per file / chunk size
my $file_number = 1;		#this is going to be a part of your new file names.
my $counter = 0;			#counts number of records

open (FASTA,  "<", "$file" )   or die "Cannot open file $file $!";

while (<FASTA>) {
	if (/^>/) {
		if ($counter++ % $record_per_file == 0) {
			my $basename = basename($file);
			my $new_file_name = $basename. $file_number++ . ".fa";
			close(NEW_FASTA);
			open(NEW_FASTA, ">", $new_file_name) or die "Cannot open file $new_file_name $!";
		}
	}
	print NEW_FASTA $_;
}

#ORIGINAL

#!/usr/bin/perl -w

my $fasta_file = "something.fasta";
my $seqs_per_file = 100;  # whatever your batch size

my $file_number = 1;  # our files will be named like "something.fasta.1"
my $seq_ctr = 0;

open(FASTA, $fasta_file) || die("can't open $fasta_file");

while(<FASTA>) {

    if(/^>/) {

       # open a new file if we've printed enough to one file
       if($seq_ctr++ % $seqs_per_file == 0) {
         close(OUT);
         open(OUT, "> " . $fasta_file . "." . $file_number++);
       }

    }

    print OUT $_;

 }
