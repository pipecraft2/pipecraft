#!/usr/bin/perl
# Keep the longest ORF per sequence from ORFfinder -outfmt 1 output and
# retain it only when its length is within the user-set min and max.
# Based on parse_orfs3.plx from MetaWorks v1.12.0 by Teresita M. Porter.
# USAGE: perl parse_ORFs.plx ORFs.nt.fasta out.fasta min_length max_length

use strict;
use warnings;

die "USAGE: perl parse_ORFs.plx ORFs.nt.fasta out.fasta min_length max_length\n" unless @ARGV == 4;
my ($infile, $outfile, $min, $max) = @ARGV;

my %ntLength; # key1 = seq ID, key2 = orf, value = length
my %ntSeq;    # key1 = seq ID, key2 = orf, value = nt sequence

open(my $in, "<", $infile) || die "Can't open $infile: $!\n";
my ($orf, $id);
while (my $line = <$in>) {
	chomp $line;
	if ($line =~ /^>/) {
		($orf, $id, my $length) = parse_nt_header($line);
		$ntLength{$id}{$orf} = $length;
		$ntSeq{$id}{$orf} = "";
	}
	else {
		$ntSeq{$id}{$orf} .= $line;
	}
}
close $in;

open(my $out, ">", $outfile) || die "Can't open $outfile: $!\n";
foreach $id (sort keys %ntLength) {
	my ($longest) = sort { $ntLength{$id}{$b} <=> $ntLength{$id}{$a} } keys %{ $ntLength{$id} };
	my $length = $ntLength{$id}{$longest};
	if ($length >= $min && $length <= $max) {
		print $out ">".$id."\n".$ntSeq{$id}{$longest}."\n";
	}
}
close $out;

# Header example: >lcl|ASV_1:3-308 ORF1_ASV_1:2:307 (minus strand ORFs have start > stop)
sub parse_nt_header {
	my ($header) = @_;
	my $orfid = (split(/ /, $header))[1];
	my ($orfotu, $start, $stop) = $orfid =~ /^(.+):(\d+):(\d+)$/;
	my ($orf, $id) = split(/_/, $orfotu, 2);
	return ($orf, $id, abs($stop - $start) + 1);
}
