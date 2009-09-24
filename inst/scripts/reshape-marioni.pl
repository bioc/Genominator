#!/usr/bin/perl
use strict;

my $infile = $ARGV[0];
# my $lane = $ARGV[1]; # We will add lanes at a later time
my ($line, $name, $chr, $start, $unique, $strand);
my $sep="|";

open INFILE, "< $infile" or die "Cannot open file!";

while($line = <INFILE>) {
    if($line =~ m/\tU.\t/) {
	chomp($line);
	($name, $chr, $start, $unique, $strand) = split /\t/, $line;
	$strand =~ s/[Ff+]/1/;
	$strand =~ s/[Rr-]/-1/;
	$chr =~ s/chr//;
	$chr =~ s/.fa//;
	if($chr =~ m/_/) {
	    next
	    }
	$chr =~ s/^X/23/;
	$chr =~ s/^Y/24/;
	$chr =~ s/^M/25/;
	$unique =~ s/U//;
	print "$chr$sep$strand$sep$start$sep$unique\n";
    }
}

close INFILE;
exit 0;
