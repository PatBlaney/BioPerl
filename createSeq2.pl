#!/usr/bin/perl
use warnings;
use strict;

use Bio::Seq;
use Bio::SeqIO;

my $seq_obj = Bio::Seq->new(
	-seq        => "aaaatgggggggggggccccgtt",
	-display_id => "#12345",
	-desc       => "example 1",
	-alphabet   => "dna"
);

my $seqio_obj = Bio::SeqIO->new(
	-file   => '>sequence.fasta',
	-format => 'fasta'
);

$seqio_obj->write_seq($seq_obj);
