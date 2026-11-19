#!/usr/bin/env perl
use strict;
use warnings;

# We need to know the library type
my $library_type = shift @ARGV or die "Usage: $0 <forward|reverse|unstranded>\n";
if ($library_type !~ /^(forward|reverse|unstranded)$/) {
    die "library_type must be: forward, reverse, or unstranded\n";
}

# Print the header with the required columns
print join("\t", qw(chr pos ref coverage variant_cov A C G T a c g t ref_count_plus ref_count_minus editing_level)), "\n";

# Process each line of the input
while (<STDIN>) {
    chomp;
    my ($chr, $pos, $ref, $depth, $bases, $quals, @rest) = split("\t", $_);

    # Skip if ref is somehow neither A nor T
    next unless uc($ref) eq 'A' || uc($ref) eq 'T';

    # Initialize base counts for alt bases only
    my %counts = (
        'A' => 0, 'C' => 0, 'G' => 0, 'T' => 0,
        'a' => 0, 'c' => 0, 'g' => 0, 't' => 0
    );

    my $ref_count_plus = 0;
    my $ref_count_minus = 0;

    # Clean up pileup read bases
    $bases =~ s/\^.//g;                          # remove start of read segment + mapping quality
    $bases =~ s/\$//g;                           # remove end of read marker
    while ($bases =~ /([+-])(\d+)/) {            # remove indels
        my $len = $2;
        $bases =~ s/[+-]\d+[ACGTNacgtn]{$len}//;
    }
    $bases =~ s/\*//g;                           # remove any *

    # Process each base in the pileup line
    foreach my $b (split('', $bases)) {
        if ($b eq '.') {
            $ref_count_plus++;  # Plus strand, reference base
        } elsif ($b eq ',') {
            $ref_count_minus++;  # Minus strand, reference base
        } elsif ($b eq 'A' || $b eq 'C' || $b eq 'G' || $b eq 'T' || $b eq 'a' || $b eq 'c' || $b eq 'g' || $b eq 't') {
            # Increment alt base counts (ignoring reference base counts for A, C, G, T)
            $counts{$b}++;
        }
    }

    # Calculate coverage based on strand and reference base
    my ($coverage, $variant_cov, $editing_level) = (0, 0, "NA");

    if (uc($ref) eq 'A') {
        my ($G, $refcount);

        if ($library_type eq "forward") {
            $G = $counts{'G'};
            $refcount = $ref_count_plus;
        }
        elsif ($library_type eq "reverse") {
            $G = $counts{'g'};
            $refcount = $ref_count_minus;
        }
        else { # unstranded
            $G = $counts{'G'} + $counts{'g'};
            $refcount = $ref_count_plus + $ref_count_minus;
        }

        $coverage = $G + $refcount;
        $variant_cov = $G;
        $editing_level = $coverage > 0 ? sprintf("%.4f", $G / $coverage) : "NA";

    }
    elsif (uc($ref) eq 'T') {

        my ($C, $refcount);

        if ($library_type eq "forward") {
            $C = $counts{'c'};
            $refcount = $ref_count_minus;
        }
        elsif ($library_type eq "reverse") {
            $C = $counts{'C'};
            $refcount = $ref_count_plus;
        }
        else { # unstranded
            $C = $counts{'C'} + $counts{'c'};
            $refcount = $ref_count_plus + $ref_count_minus;
        }

        $coverage = $C + $refcount;
        $variant_cov = $C;
        $editing_level = $coverage > 0 ? sprintf("%.4f", $C / $coverage) : "NA";

    }

    # Print the results in the desired format
    print join("\t", $chr, $pos, $ref, $coverage, $variant_cov,
        @counts{qw(A C G T a c g t)},  # Base counts for alt bases (A, C, G, T, a, c, g, t)
        $ref_count_plus, $ref_count_minus,  # Reference counts for plus and minus strands
        $editing_level  # Editing level (if applicable)
    ), "\n";
}

