#!/usr/bin/env perl
# =============================================================================
# mismatch_spectrum_count.pl
#
# Strand-aware base-substitution counts from a single-sample `samtools mpileup`
# stream (stdin), for stranded 10x 3' libraries (read strand = RNA strand,
# i.e. parse_pileup_stranded.pl "forward"). Every aligned base that passed the
# mpileup -q/-Q filters is counted once, on the RNA strand:
#   + strand read ('.' or upper case): RNA ref = ref,             RNA base = base
#   - strand read (',' or lower case): RNA ref = complement(ref), RNA base = complement(base)
# Deletions (*, #), N bases, and reference N positions are skipped.
#
# Usage : samtools mpileup ... | perl mismatch_spectrum_count.pl <label>
# Output: TSV to stdout: label, rna_ref, rna_base, n  (rna_ref == rna_base = matches)
# Used by: human/pm/scripts/run_gt_truth_rebuild_v6.sh (stages spectrum*)
# =============================================================================
use strict;
use warnings;

my $label = shift @ARGV // "sample";
my %comp = (A => "T", C => "G", G => "C", T => "A");
my %n;

while (my $line = <STDIN>) {
  chomp $line;
  my @f = split /\t/, $line;
  next if @f < 5;
  my $ref = uc $f[2];
  next unless exists $comp{$ref};
  my $s = $f[4];
  # drop read-start markers (^ + mapping-quality char) and read-end markers ($)
  $s =~ s/\^.//g;
  $s =~ s/\$//g;
  # drop indel sequences (+3ACG / -2TT): the length digits say how many to skip
  $s =~ s/[+-](\d+)(??{ ".{$^N}" })//g if $s =~ /[+-]/;
  my $rc = $comp{$ref};
  # count by character class ('*', '#', 'N', 'n', '>', '<' are not counted)
  $n{$ref}{$ref} += ($s =~ tr/.//);
  $n{$rc}{$rc}   += ($s =~ tr/,//);
  $n{$ref}{A} += ($s =~ tr/A//); $n{$ref}{C} += ($s =~ tr/C//);
  $n{$ref}{G} += ($s =~ tr/G//); $n{$ref}{T} += ($s =~ tr/T//);
  $n{$rc}{T}  += ($s =~ tr/a//); $n{$rc}{G}  += ($s =~ tr/c//);
  $n{$rc}{C}  += ($s =~ tr/g//); $n{$rc}{A}  += ($s =~ tr/t//);
}

for my $r (sort keys %n) {
  for my $b (sort keys %{$n{$r}}) {
    print join("\t", $label, $r, $b, $n{$r}{$b}), "\n";
  }
}
