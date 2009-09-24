#!/usr/bin/perl
use strict;

my $infile = $ARGV[0];
my $outfile = $ARGV[1];
my ($line, $name, $chr, $start, $unique, $strand);
my $sep="|";
my @args;

open INFILE, "< $infile" or die "Cannot open file!";
open OUTFILE, "> $outfile.reshaped" or die "Cannot open file!";

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
	print OUTFILE "$chr$sep$strand$sep$start$sep$unique\n";
    }
}

close INFILE;
close OUTFILE;

@args = ("sqlite3", "$outfile.db", "CREATE TABLE rawReads (chr INTEGER, strand INTEGER, location INTEGER, quality INTEGER);");
system(@args) == 0 || die "system @args failed: $?\n";

@args = ("sqlite3", "$outfile.db", ".import $outfile.reshaped rawReads");
system(@args) == 0 || die "system @args failed: $?\n";

@args = ("sqlite3", "$outfile.db", "CREATE INDEX rawReadsIDX ON rawReads (chr, strand, location)");
system(@args) == 0 || die "system @args failed: $?\n";

@args = ("rm", "$outfile.reshaped");
system(@args) == 0 || die "system @args failed: $?\n";

exit 0;
