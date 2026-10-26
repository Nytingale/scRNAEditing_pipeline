#!/usr/bin/env perl
# =============================================================================
# mismatch_site_count.pl
#
# Per-site, strand-aware counts for ALL 12 mismatch types from a single-sample
# `samtools mpileup` stream (stdin), for stranded 10x 3' libraries (read strand
# = RNA strand; parse_pileup_stranded.pl "forward"). Read-string cleaning is the
# same as parse_pileup_stranded.pl; the per-type definitions generalise its A>G:
#   '+' strand: RNA ref = ref,             ref reads = '.', alt Y reads = upper-case Y
#   '-' strand: RNA ref = complement(ref), ref reads = ',', alt Y reads = lower-case complement(Y)
#   coverage = ref reads + Y reads on that strand; variant_cov = Y reads;
#   editing_level = variant_cov / coverage
# so for type AG this reproduces parse_pileup_stranded.pl forward exactly
# (ref A: '+'; ref T: '-').
#
# Usage : samtools mpileup ... | perl mismatch_site_count.pl scan [MIN_COV MIN_ALT MIN_EDIT]
#           -> candidate rows chr pos strand rna_ref rna_alt coverage variant_cov editing_level
#              for every (strand, alt) with coverage >= MIN_COV (5), variant_cov >= MIN_ALT (1)
#              and editing_level >= MIN_EDIT (0.01)   [the v6 scan floor, per type]
#         samtools mpileup ... | perl mismatch_site_count.pl site
#           -> chr pos ref ref_plus ref_minus A C G T a c g t (raw counts, every position)
# Used by: human/pm/scripts/run_gt_spectrum12_v6.sh
# =============================================================================
use strict;
use warnings;

my $mode = shift @ARGV or die "Usage: $0 <scan|site> [MIN_COV MIN_ALT MIN_EDIT]\n";
die "mode must be scan or site\n" unless $mode eq "scan" || $mode eq "site";
my ($min_cov, $min_alt, $min_edit) = (shift(@ARGV) // 5, shift(@ARGV) // 1, shift(@ARGV) // 0.01);
my %comp = (A => "T", C => "G", G => "C", T => "A");

if ($mode eq "scan") {
  print join("\t", qw(chr pos strand rna_ref rna_alt coverage variant_cov editing_level)), "\n";
} else {
  print join("\t", qw(chr pos ref ref_plus ref_minus A C G T a c g t)), "\n";
}

while (my $line = <STDIN>) {
  chomp $line;
  my ($chr, $pos, $ref, $depth, $bases) = split /\t/, $line;
  next unless defined $bases;
  $ref = uc $ref;
  next unless exists $comp{$ref};

  $bases =~ s/\^.//g;
  $bases =~ s/\$//g;
  while ($bases =~ /([+-])(\d+)/) {
    my $len = $2;
    $bases =~ s/[+-]\d+[ACGTNacgtn]{$len}//;
  }
  $bases =~ s/\*//g;

  my $rp = ($bases =~ tr/.//);
  my $rm = ($bases =~ tr/,//);
  my %c = (A => ($bases =~ tr/A//), C => ($bases =~ tr/C//), G => ($bases =~ tr/G//), T => ($bases =~ tr/T//),
           a => ($bases =~ tr/a//), c => ($bases =~ tr/c//), g => ($bases =~ tr/g//), t => ($bases =~ tr/t//));

  if ($mode eq "site") {
    print join("\t", $chr, $pos, $ref, $rp, $rm, @c{qw(A C G T a c g t)}), "\n";
    next;
  }

  # '+' strand: RNA ref = ref, alt Y = upper-case Y
  for my $y (qw(A C G T)) {
    next if $y eq $ref;
    my $alt = $c{$y};
    next if $alt < $min_alt;
    my $cov = $rp + $alt;
    next if $cov < $min_cov || $alt / $cov < $min_edit;
    printf "%s\t%s\t+\t%s\t%s\t%d\t%d\t%.4f\n", $chr, $pos, $ref, $y, $cov, $alt, $alt / $cov;
  }
  # '-' strand: RNA ref = comp(ref), RNA alt Y = lower-case comp(Y)
  my $rref = $comp{$ref};
  for my $y (qw(A C G T)) {
    next if $y eq $rref;
    my $alt = $c{lc $comp{$y}};
    next if $alt < $min_alt;
    my $cov = $rm + $alt;
    next if $cov < $min_cov || $alt / $cov < $min_edit;
    printf "%s\t%s\t-\t%s\t%s\t%d\t%d\t%.4f\n", $chr, $pos, $rref, $y, $cov, $alt, $alt / $cov;
  }
}
