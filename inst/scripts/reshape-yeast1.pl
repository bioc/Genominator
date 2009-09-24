#!/usr/bin/perl
use strict;

my $infile = $ARGV[0];
my ($line, $name, $chr, $start, $unique, $strand, $seq, $crap1, $crap2, $crap3);
my $sep="|";

open INFILE, "< $infile" or die "Cannot open file!";

while($line = <INFILE>) {
    if($line =~ m/\tU.\t/) {
	chomp($line);
	($name, $seq, $unique, $crap1, $crap2, $crap3, $chr, $start, $strand) = split /\t/, $line;
	if($unique =~ m/[NnQqRr]/) {
	    next
	    }
	$strand =~ s/[Ff+]/1/;
	$strand =~ s/[Rr-]/-1/;
	if($stand == -1) {
	    $start = $start - 11;
	}
	$chr =~ s/chr//;
	$chr =~ s/.fsa//;
	$chr =~ s/^0//;
	$chr =~ s/mt/17/;
	$unique =~ s/U//;
	print "$chr$sep$strand$sep$start$sep$unique\n";
    }
}

close INFILE;
exit 0
