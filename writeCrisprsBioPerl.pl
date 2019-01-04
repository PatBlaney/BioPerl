#!/usr/bin/perl
use warnings;
use strict;
use diagnostics;

# username: PatrickBlaney

# Original file framework belongs to Chuck Roesel
# Use BioPerl to read a fasta file, parse the sequences, identify all possible CRISPR sequences
# and write them to a new fasta output file

use Bio::Seq;
use Bio::SeqIO;
use Getopt::Long;
use Pod::Usage;

# Initialize variable that hold the user-specifed option for the program
my $fastaIn = '';
my $usage   = "\n$0 [options] \n
Options:
	-fastaIn			FASTA file with sequence
	-help 				Show this message
\n";

# Utilize the GetOptions subroutine from the Getopt::Long module to store the program
# option to accept
GetOptions(
	'fastaIn=s' => \$fastaIn,
	'help'      => sub { pod2usage($usage); },
) or pod2usage($usage);

# Use the BioPerl module SeqIO to read the fasta file
my $seqioObj = Bio::SeqIO->new(
	-file   => "$fastaIn",
	-format => 'fasta'
);

# Hash to store kmers
my %kMerHash = ();

# Hash to store occurrences of last 12 positions
my %last12Counts = ();

# Loop through each entry within the fasta file and use the subroutine to identify all k-mers that
# could potentially be CRISPR sequences to then store in a hash
while ( my $seqLine = $seqioObj->next_seq ) {
	kmerFinder( $seqLine->seq, \%kMerHash, \%last12Counts );
}

# Subroutine used to slide along each sequence read from the fasta file to identify
# k-mers of a specific length
sub kmerFinder {

	# Use the sequence read from fasta file as the parameter
	my ( $sequenceString, $kMerHash, $last12Counts ) = @_;

	# Set the size of the sliding window
	my $windowSize = 21;

	# Set the step size
	my $stepSize  = 1;
	my $seqLength = length($sequenceString);

	# Test to see if subroutine works with loop
	#print $seqLength, "\n";

	# For loop to increment the starting position of the sliding window
	# Starts at position zero; doesn't move past end of file; advance the window
	# by step size
	for (
		my $windowStart = 0 ;
		$windowStart <= ( $seqLength - $windowSize ) ;
		$windowStart += $stepSize
	  )
	{

	  	# Get a 21-mer substring from sequenceRef (two $ to deference reference to
	  	# sequence string) starting at the window start for length $windowStart
		my $crisprSeq = substr( $sequenceString, $windowStart, $windowSize );

		# If the 21-mer ends in GG, create a hash with key=last 12 of k-mer and value
		# is 21-mer. Regex where $1 is the crispr, and $2 contains the last 12 of crispr.
		if ( $crisprSeq =~ /([ATGC]{9}([ATGC]{10}GG))$/ ) {

			# Put the crispr in the hash with last 12 as key, full 21 as value.
			$$kMerHash{$2} = $1;
			$$last12Counts{$2}++;

			# Test to see if correct string is being stored
			#print $last12Counts{$2}, "\n";

		}
	}
}

# Test to see hashes are properly populated
#my $hashsize = %kMerHash;
#print $hashsize, "\n";
#my $hashsize2 = %last12Counts;
#print $hashsize2, "\n";
#print values(%last12Counts);

# Open a new fasta file to write the CRISPR sequence to
my $seqioObjOutput = Bio::SeqIO->new(
	-file   => '>crisprs1.fasta',
	-format => 'fasta'
);

# Initialize the CRISPR count to zero
my $crisprCount = 0;

# Loop through the hash of last 12 counts
for my $last12Seq ( sort ( keys %last12Counts ) ) {

	# Check if count == 1 for this sequence
	if ( $last12Counts{$last12Seq} == 1 ) {

		# The last 12 seq of this CRISPR is unique in the genome.
		# Increment the CRISPR count.
		$crisprCount++;

		# Create a unique object for each CRISPR sequence
		my $seqObj = Bio::Seq->new(
			-seq        => "$kMerHash{$last12Seq}",
			-alphabet   => "dna",
			-display_id => ">crispr_$crisprCount",
			-desc       => "CRISPR"
		);

		# Write the CRISPR sequence to the output file
		$seqioObjOutput->write_seq($seqObj);
	}
}
