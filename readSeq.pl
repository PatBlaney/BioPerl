#!/usr/bin/perl
use warnings;
use strict;

use Bio::Seq;
use Bio::SeqIO;

my $seqio_obj = Bio::SeqIO->new(
		-file => 'dmel-all-chromosome-r6.17.fasta',
		-format => 'fasta'
);

while (my $seq_obj = $seqio_obj->next_seq){
	print $seq_obj->desc, "\n";
}