#!/usr/bin/perl
# Teresita M. Porter, Feb. 28, 2020 [MetaWorks v1.12.0]
	# slightly modified for PipeCraft2 by Sten Anslan, August 30, 2023.
# Script to parse hmm output in ORFfinder_hmm.sh

use strict;
use warnings;
use Data::Dumper;
use Statistics::Descriptive;

# declare var
my $i=0;
my $line;
my $id;
my $score;
my $size;
my $stat;
my $percentile25;
my $percentile75;
my $iqr;
my $lc;
my $uc;
my $j;
my $seq;
my $record;

# declare array
my @hmm;
my @orfs;
my @rdp;
my @line;
my @scores;
my @ids;

# declare hash
my %hmm; #key=id, value=score
my %orfs; #key=id, value=seq

open (IN1, "<", $ARGV[0]) || die "Error can't open infile1: $!\n";
@hmm=<IN1>;
close IN1;

open (IN2, "<", $ARGV[1]) || die "Error can't open infile2: $!\n";
@orfs=<IN2>;
close IN2;

# parse HMMER output and hash scores
while ($hmm[$i]) {
	$line = $hmm[$i];
	chomp $line;

	if ($line =~ /^#/) {
		$i++;
		next;
	}
	else {
		@line = split ' ', $line; # split on whitespace (of any size)
		$id = $line[2];
		$score = $line[5];
		push(@scores, $score);
		$hmm{$id} = $score;

		$i++;
		next;
	}
}
$i=0;

# figure out cutoff for outlier hmmer scores
$stat = Statistics::Descriptive::Full->new();
$stat->add_data(@scores);
$percentile25 = $stat -> percentile(25);
$percentile75 = $stat -> percentile(75);
$iqr = $percentile75-$percentile25;
$lc = $percentile25-($iqr*1.5);
$uc = $percentile75+($iqr*1.5);

# get list of good ids (skip over ids with outlier scores)
while ( ($id, $score) = each(%hmm) ) {
	unless ($score < $lc || $score > $uc) {
		push(@ids, $id);
		print $id."\n";
	}
}

# parse filtered nt orfs and hash seqs
while ($orfs[$i]) {
	$line = $orfs[$i];
	chomp $line;
	if ($line =~ /^>/) {
		$line =~ s/^>//;
		$id = $line;
		$j = $i+1;
		$seq = $orfs[$j];
		chomp $seq;
		$orfs{$id} = $seq;
		$i+=2;
		next;
	}
	else {
		$i++;
	}
}
$i=0;
