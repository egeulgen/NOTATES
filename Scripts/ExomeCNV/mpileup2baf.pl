#!/usr/bin/perl
use strict;
use warnings;
use Getopt::Long qw(:config pass_through no_ignore_case);

my ($min_reads) = (10);
GetOptions (
        "min-reads:s" => \$min_reads,
);

while (<>) {
        my $line = $_;
        my @columns = split("\t",$line);
        my $chr = $columns[0];
        my $position = $columns[1];
        my $num_reads = $columns[3];
        my $calls = $columns[4];
        if($num_reads < $min_reads){ # not enough coverage to have good confidence in the call
                next;
        }
        my $num_ref = 0;
        while ($calls =~ /[,.]/g) { $num_ref++ }
        my $num_var = $num_reads - $num_ref;
        print("$chr\t$position\t$num_reads\t$num_var\n");
}