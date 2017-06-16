#!/usr/bin/perl
# extract reads mapped around coordinate if
# they or their mates map in region
use strict;
use warnings;

my $usage = 'cat aln.sam | extract.pl <chr> <coord> <range>';

my $chr = shift @ARGV or die $usage;
my $coord = shift @ARGV or die $usage;
my $pad = shift @ARGV or die $usage;

my @reads;
while (<>) {
	if (/^\@/) {
		print;
		next;
	}

	my ($read, $chra, $coorda, $chrb, $coordb) = (split)[0,2,3,6,7];
	$chrb = $chrb eq '='? $chra : $chrb;	

	my $self = ($chra eq $chr) &&
		($coorda<=$coord+$pad && $coorda>=$coord-$pad);
	my $mate = ($chrb eq $chr) &&
		($coordb<=$coord+$pad && $coordb>=$coord-$pad);

	if ($self || $mate) {
		print;
		push @reads, $read;
		next;
	}
}
