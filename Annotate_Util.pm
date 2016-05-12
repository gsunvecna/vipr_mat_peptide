package Annotate_Util;

use strict;
use warnings;
use English;
use Carp;
use Data::Dumper;

use version; our $VERSION = qv('1.1.1'); # October 13 2010
use File::Temp qw/ tempfile tempdir /;
use Bio::SeqIO;
use Bio::Seq;
use Bio::AlignIO;
use Bio::Tools::Run::StandAloneBlast;
use IO::String;

use Annotate_Verify;

my $debug_all = 1;

####//README//####
#
# Annotate_Util contains the core functions to perform annotation based on both MUSCLE and bl2seq alignment
#
#    Authors Chris Larsen, clarsen@vecna.com; Guangyu Sun, gsun@vecna.com
#    April 2010
#
##################

## //EXECUTE// ##

# Get refseq object if $refseq_fn is given
#  my $dbh_ref;
#  my $refseq = &get_refseq($refseq_fn, $dir_path);

my $taxon = { loaded =>0, fn =>"viprbrc_taxon_records.txt" };

=head2 generate_fasta

Takes an array of SeqFeature.
 Write the translations to a virtual file
 Returns the fasta file in a string

=cut

sub generate_fasta {
    my ($feats_all) = @_;

    my $debug = 0 && $debug_all;

    my $seq_out;
    my $faa1 = '';
    my $vfile = IO::String->new($faa1);
    $seq_out = Bio::SeqIO->new(
                  '-fh'     => $vfile,
                  '-format' => 'fasta'
                              );

    foreach my $feat (@$feats_all) {

        next if ($feat->primary_tag eq 'CDS'); # Exclude CDS
        my @values = $feat->get_tag_values('translation');
        my $s = $values[0];
        my $desc;
        if ($feat->has_tag('note')) {
            @values = $feat->get_tag_values('note');
            foreach my $value (@values) {
                $desc = $1 if ($value =~ /^Desc:(.+)$/i);
            }
        } elsif($feat->has_tag('product')) {
            $desc = 'ACC='.$feat->seq->accession_number;
            @values = $feat->get_tag_values('product');
            $desc .= '|product='.$values[0];
            $desc .= '|Loc='.$feat->location->to_FTstring;
        }
        my $f = Bio::PrimarySeq->new(
                         -seq      => $s,
                         -id       => '',	# id can't contain space
                         -desc     => $desc,	# desc can contain space
                         -alphabet => 'protein'
                                     );
        if ($f->alphabet eq 'dna') {
            $seq_out->write_seq($f->translate());
        } elsif($f->alphabet eq 'protein') {
            $seq_out->write_seq($f);
        }

    }

    return $faa1;
} # sub generate_fasta


=head2 msa_get_feature_loc_bl2seq

Takes reffeat, gaps in the reffeat, and gaps (in DNA coordinate) in the target
 calculate the location of the feature in target genome
 Return a Bio::Location object for the new annotation

Note: We can only deal with Bio::Location::Simple, not more complicated location types

=cut

sub msa_get_feature_loc_bl2seq {
    my ($reffeat, $refcds_loc, $cds_loc, $gaps_q, $gaps_h) = @_;

    my $debug = 0 && $debug_all;
    my $subname = 'msa_get_feature_loc_bl2seq';

    my $refcds_start = $refcds_loc->start;
    my $cds_start = $cds_loc->start;
    my @gaps_q = @$gaps_q;
    my @gaps_h = @$gaps_h;
#    $debug && print STDERR "$subname: \@gaps_q = '@gaps_q'\n";
#    $debug && print STDERR "$subname: \@gaps_h = '@gaps_h'\n";

    my $loc2 = undef;
    if (!$reffeat->location->isa('Bio::Location::Simple')) {
        print STDERR "$subname: not a Bio::Location::Simple, skip.\n";
        return $loc2;
    }

    my $start;	# these location start and end refer to original DNA seq
    my $end;
    $start = $reffeat->location->start;
    $end = $reffeat->location->end;
    $debug && print STDERR "$subname: \$start=$start \$end=$end\n";

    # need to take care of any gaps in the alignment, basically add 1 if there is a gap in query,
    # or subtract 1 if there is a gap in hit(target)
    my $adjust = {
                   'start_q' => 0,
                   'end_q'   => 0,
                   'start_h' => 0,
                   'end_h'   => 0,
                 }; # hash, containing numbers of gaps
    if ( 0 ) {
#      if ($#gaps_q || $#gaps_h) {
        $debug && $gaps_q->[1] && print STDERR "$subname: \@gaps_q='@gaps_q'\n";
        $debug && $gaps_h->[1] && print STDERR "$subname: \@gaps_h='@gaps_h'\n";

        $debug && print STDERR "$subname: Checking \$start=$start with \@gaps_q=@gaps_q\n";
        $adjust->{start_q} = &adjust_from_gaps($start, $gaps_q);

        $start = $start - $refcds_start + $cds_start;

        $debug && print STDERR "$subname: Checking \$start=$start with \@gaps_h=@gaps_h\n";
        $adjust->{start_h} = &adjust_from_gaps($start, $gaps_h);

        $start = $start + $adjust->{start_q} - $adjust->{start_h};

        $debug && print STDERR "$subname: Checking \$end=$end with \@gaps_q=@gaps_q\n";
        $adjust->{end_q} = &adjust_from_gaps($end, $gaps_q);

        $end   = $end   - $refcds_start + $cds_start;

        $debug && print STDERR "$subname: Checking \$end=$end with \@gaps_h=@gaps_h\n";
        $adjust->{end_h} = &adjust_from_gaps($end, $gaps_h);

        $end = $end + $adjust->{end_q} - $adjust->{end_h};

        $debug && print STDERR "$subname: \$adjust=\n".Dumper($adjust)."end of \$adjust\n\n";

 #     }
    } else {
#      if ($#gaps_q || $#gaps_h) {
        $debug && $gaps_q->[1] && print STDERR "$subname: \@gaps_q='@gaps_q'\n";
        $debug && $gaps_h->[1] && print STDERR "$subname: \@gaps_h='@gaps_h'\n";

        for (my $i = 1; $i<=$#{@$gaps_q}; $i++) {
#            $debug && print STDERR "$subname: \$start=$start  + $adjust->{'start_q'}\n";
#            $debug && print STDERR "$subname: \$start=$start \t\$i=$i \$gap=$gaps_q->[$i]\n";
            if ($gaps_q->[$i] < ($start + $adjust->{'start_q'})) { # This $position doesn't have any gap
                $adjust->{'start_q'} = $i;
            } else {
                last;
            }
        }
        $start = $start - $refcds_start + $cds_start;
        for (my $i = 1; $i<=$#{@$gaps_h}; $i++) {
#            $debug && print STDERR "$subname: \$start=$start  + $adjust->{'start_q'} - $adjust->{'start_h'}\n";
#            $debug && print STDERR "$subname: \$start=$start \t\$i=$i \$gap=$gaps_h->[$i]\n";
            if ($gaps_h->[$i] < ($start + $adjust->{'start_q'})) { # This $position doesn't have any gap
                $adjust->{'start_h'} = $i;
            } else {
                last;
            }
        }

        for (my $i = 1; $i<=$#{@$gaps_q}; $i++) {
#            $debug && print STDERR "$subname: \$end=$end + $adjust->{'end_q'}\n";
#            $debug && print STDERR "$subname: \$end=$end \t\$i=$i \$gap=$gaps_q->[$i]\n";
            if ($gaps_q->[$i] < ($end + $adjust->{'end_q'})) { # This $position doesn't have any gap
                $adjust->{'end_q'} = $i;
            } else {
                last;
            }
        }
        $end   = $end   - $refcds_start + $cds_start;
        for (my $i = 1; $i<=$#{@$gaps_h}; $i++) {
#            $debug && print STDERR "$subname: \$end=$end  + $adjust->{'end_q'} - $adjust->{'end_h'}\n";
#            $debug && print STDERR "$subname: \$end=$end \t\$i=$i \$gap=$gaps_h->[$i]\n";
            if ($gaps_h->[$i] < ($end + $adjust->{'end_q'})) { # This $position doesn't have any gap
                $adjust->{'end_h'} = $i;
            } else {
                last;
            }
        }

#      }
    }
    $debug && print STDERR "$subname: \$start=$start \$end=$end\n";
    $debug && print STDERR "$subname: $start += $adjust->{start_q} - $adjust->{start_h}\n";
    $debug && print STDERR "$subname: $end   += $adjust->{end_q}   - $adjust->{end_h}\n";
    $start = $start + $adjust->{start_q} - $adjust->{start_h};
    $end   = $end   + $adjust->{end_q}   - $adjust->{end_h};
    $debug && print STDERR "$subname: \$start=$start \$end=$end\n";


    # finally, adjust the beginning of refcds and cds
    $debug && print STDERR "$subname: \$start = $start - $refcds_start + $cds_start\n";
    $debug && print STDERR "$subname: \$end   = $end   - $refcds_start + $cds_start\n";

    $debug && print STDERR "$subname: \$cds_loc is ".ref($cds_loc).' '.$cds_loc->to_FTstring."\n";

    # This just checks if $loc2 is within the CDS
    if ($cds_start<=$start && $cds_start<$end && $start<$cds_loc->end && $end<=$cds_loc->end) {
        $loc2 = Bio::Location::Simple->new();
        $loc2->start($start);
        $loc2->end  ($end);
        $debug && print STDERR "$subname: \$cds_loc=\n".Dumper($cds_loc)."end of \$cds_loc\n\n";
        $debug && print STDERR "$subname: \$loc2=\n".Dumper($loc2)."end of \$loc2\n\n";
    } else {
        $debug && print STDERR "$subname: ($start, $end) ($cds_start, ".$cds_loc->end.")\n";
        $debug && print STDERR "$subname: \$loc2=doesn't look good, skip.\n";
        $loc2 = undef;
    }

    # further check the validity of $loc2, this ensures the new annotated feature lie within the cds,
    # and doesn't have more gap than 3 dna residues
    if ($loc2) {
        $loc2 = &msa_check_feature_loc($loc2, $reffeat, $refcds_loc, $cds_loc, $adjust, $gaps_q, $gaps_h);
    }

    return ($loc2);

} # sub msa_get_feature_loc_bl2seq


=head2 get_adjust_start

Takes reffeat, gaps in the reffeat, and gaps (in DNA coordinate) in the target
 calculate the location of the feature in target genome
 Return a Bio::Location object for the new annotation

Note: We can only deal with Bio::Location::Simple, not more complicated location types

=cut
=head1
sub get_adjust_start {
    my ($start, $refcds_start, $cds_start, $gaps_q, $gaps_h) = @_;

    my $debug = 1 && $debug_all;
    my $subname = 'get_adjust_start';
    my ($under_water, $above_water);

    my $adjust = {'start_q' => 0, 'start_h' => 0};

    return undef if (!$start);

    $debug && print STDERR "$subname: \$start=$start, \$refcds_start=$refcds_start, \$cds_start=$cds_start, $gaps_q, $gaps_h\n";
    $debug && print STDERR "$subname: \$gaps_q:$#{@$gaps_q}, \$gaps_h:$#{@$gaps_h}\n";

    for (my $i = 1; $i<=$#{@$gaps_q}; $i++) {
        $debug && print STDERR "$subname: \$start=$start  + $adjust->{'start_q'}";
        $debug && print STDERR " \$start=$start \t\$i=$i \$gap=$gaps_q->[$i]\n";
        if ($gaps_q->[$i] <= ($start + $adjust->{'start_q'})) {
            $adjust->{'start_q'} = $i;
        } else {
            last;
        }
    }
    ($under_water, $above_water) = (0, 0);
    for (my $i=1; $i<=$#{@$gaps_h}; $i++) {
        if ($gaps_h->[$i]<=0) {
            $under_water++;
        } elsif ($gaps_h->[$i]+$under_water==$i) {
            $above_water++;
        }
    }
    $debug && print STDERR "$subname: \$under_water=$under_water \$above_water=$above_water\n";
    $start = $start - $refcds_start + $cds_start + $adjust->{'start_q'} -$under_water;
    for (my $i = 1; $i<=$#{@$gaps_h}; $i++) {
        $debug && print STDERR "$subname: \$start=$start  - $adjust->{'start_h'}";
        $debug && print STDERR " \$start=$start \t\$i=$i \$gap=$gaps_h->[$i]\n";
        # Don't change this. EF407494 with refseq=NC_004102
        if ($gaps_h->[$i] < ($start)) {
            $adjust->{'start_h'} = $i;
        } else {
            last;
        }
    }
    $debug && print STDERR "\n";

    return ($adjust->{'start_q'}, $adjust->{'start_h'});

} # sub get_adjust_start
=cut

=head2 get_adjust_end

Takes reffeat, gaps in the reffeat, and gaps (in DNA coordinate) in the target
 calculate the location of the feature in target genome
 Return a Bio::Location object for the new annotation

Note: We can only deal with Bio::Location::Simple, not more complicated location types

=cut
=head1
sub get_adjust_end {
    my ($end, $refcds_start, $cds_start, $gaps_q, $gaps_h, $adjust) = @_;

    my $debug = 1 && $debug_all;
    my $subname = 'get_adjust_end';
    my ($under_water, $above_water);

    $adjust->{'end_q'} = 0;
    $adjust->{'end_h'} = 0;
    return undef if (!$end);

    for (my $i = 1; $i<=$#{@$gaps_q}; $i++) {
#        $debug && print STDERR "$subname: \$end=$end + $adjust->{end_q} \$end=$end \t\$i=$i \$gap=$gaps_q->[$i]\n";
        if ($gaps_q->[$i] <= ($end + $adjust->{'end_q'})) {
            $adjust->{'end_q'} = $i;
        } else {
            last;
        }
    }
    $debug && print STDERR "$subname: \n";
    ($under_water, $above_water) = (0, 0);
    for (my $i=1; $i<=$#{@$gaps_h}; $i++) {
        if ($gaps_h->[$i]<=0) {
            $under_water++;
        } elsif ($gaps_h->[$i]+$under_water==$i) {
            $above_water++;
        }
#        $debug && print STDERR "$subname: \$gaps_h->[$i]=$gaps_h->[$i] \$under_water=$under_water \$above_water=$above_water\n";
    }
    $debug && print STDERR "$subname: \$under_water=$under_water \$above_water=$above_water\n";
    $end   = $end   - $refcds_start + $cds_start + $adjust->{'end_q'} -$under_water;
    for (my $i = 1; $i<=$#{@$gaps_h}; $i++) {
        $debug && print STDERR "$subname: \$end=$end  - $adjust->{'end_h'} ";
        $debug && print STDERR " \$end=$end \t\$i=$i \$gap=$gaps_h->[$i]\n";

        if ($gaps_h->[$i] <= ($end)) { # AF054258 with refseq=NC_004102
            $adjust->{'end_h'} = $i;
        } else {
            last;
        }
    }

    return ($adjust->{'end_q'}, $adjust->{'end_h'});

} # sub get_adjust_end
=cut


sub adjust_start {
    my ($qstart, $refcds_loc, $cds, $aln, $aln_q, $aln_h) = @_;

    return undef if (!$qstart);
    my $debug = 0 && $debug_all;
    my $subname = 'adjust_start';

    my $cds_loc = $cds->location;
    # Convert to protain from DNA
    $qstart = ($qstart - $refcds_loc->start)/3 +1;
    $debug && print STDERR "$subname: \$qstart=$qstart\n";
#    $debug && print STDERR "$subname: \$aln_q=\n".Dumper($aln_q)."end of $aln_q\n";
    my $qseq = $aln_q->seq;
    $debug && print STDERR "$subname: \$qseq='$qseq' \$qseq=".length($qseq)."\n";
    my $gap_char = $aln->gap_char;
    my $c;
    for (my $n=1; $n<=$qstart && $n<=length($aln_q->seq); $n++) {
        $c = substr($qseq, $n-1, 1);
        $qstart++ if ($c eq $gap_char);
        if (($n-1)%10 ==0) {
            $debug && print STDERR "\n$subname: \$n=$n \$gap_char=$gap_char ";
        }
        $debug && print STDERR "\$c=$c \$qstart=$qstart ";
    }
    $debug && print STDERR "\n";

    my $hstart = $qstart;
    my $hseq = $aln_h->seq;
    $hseq = substr($hseq, 0, $hstart-1);
    $debug && print STDERR "$subname: \$hseq='$hseq' \$hseq=".length($hseq)."\n";
    for (my $n=length($hseq);  $n>=1; $n--) {
        $c = substr($hseq, $n-1, 1);
        $hstart-- if ($c eq $gap_char);
        if (($n-0)%10 ==0 || $n==length($hseq)) {
            $debug && print STDERR "\n$subname: \$n=$n \$gap_char=$gap_char ";
        }
        $debug && print STDERR "\$c=$c \$hstart=$hstart ";
    }
    $debug && print STDERR "\n";
    $debug && print STDERR "$subname: \$hstart=$hstart\n";

    # Convert back to DNA
    $hstart = $cds_loc->start + ($hstart-1)*3;

    return $hstart;
} # sub adjust_start

sub adjust_end {
    my ($qend, $refcds_loc, $cds, $aln, $aln_q, $aln_h) = @_;

    return undef if (!$qend);
    my $debug = 0 && $debug_all;
    my $subname = 'adjust_end';

    my $last_mat_peptide = 0;
    $last_mat_peptide = 1 if ($qend == $refcds_loc->end || $qend == $refcds_loc->end-3);
    $debug && print STDERR "$subname: \$last_mat_peptide=$last_mat_peptide\n";

    my $cds_loc = $cds->location;
    # Convert to protain from DNA
    $qend = ($qend+1 - $refcds_loc->start)/3 -1 +1;
    $debug && print STDERR "$subname: \$qend=$qend\n";
#    $debug && print STDERR "$subname: \$aln_q=\n".Dumper($aln_q)."end of $aln_q\n";
    my $qseq = $aln_q->seq;
#    $qseq = substr($qseq, 0, $qend);
    $debug && print STDERR "$subname: \$qseq='$qseq' \$qseq=".length($qseq)."\n";
    my $gap_char = $aln->gap_char;
    my $c;
    for (my $n=1; $n<=$qend && $n<=length($aln_q->seq); $n++) {
        $c = substr($qseq, $n-1, 1);
        $qend++ if ($c eq $gap_char);
        if (($n-1)%10 ==0) {
            $debug && print STDERR "\n$subname: \$n=$n \$gap_char=$gap_char ";
        }
        $debug && print STDERR "\$c=$c \$qend=$qend ";
    }
    $debug && print STDERR "\n";
    $debug && print STDERR "$subname: \$qend=$qend\n";

    # If this is the last mat_peptide within refseq CDS, see if there is any trailing short sequence in target.
    # If so, and the length is 5 AA or less, include it in the annotation.
    if ($last_mat_peptide) {
        if ($qend<length($qseq) && length($qseq)-$qend<=5) {
            print STDERR "$subname: \$qend=$qend changed to include the trailing sequence after last mat_peptide\n";
            $qend = length($qseq);
            print STDERR "$subname: \$qend=$qend changed to include the trailing sequence after last mat_peptide\n";
        }
    }

    my $hend = $qend;
    my $hseq = $aln_h->seq;
    $hseq = $aln_h->seq;
    $hseq = substr($hseq, 0, $hend);
    $debug && print STDERR "$subname: \$hseq='$hseq' \$hseq=".length($hseq)."\n";
    for (my $n=length($hseq);  $n>=1; $n--) {
        $c = substr($hseq, $n-1, 1);
        $hend-- if ($c eq $gap_char);
        if (($n-0)%10 ==0 || $n==length($hseq)) {
            $debug && print STDERR "\n$subname: \$n=$n \$gap_char=$gap_char ";
        }
        $debug && print STDERR "\$c=$c \$hend=$hend ";
    }
    $debug && print STDERR "\n";

    # Convert back to DNA
    $hend = $cds_loc->start + ($hend-1)*3 +2;

    my %good_tail = ( # if CDS ends with any of these DNA sequences, an extra amino acid is decided
                      'CU'=>1, 'CT'=>1, # (Leu/L) Leucine
                      'GU'=>1, 'GT'=>1, # (Val/V) Valine
                      'UC'=>1, 'TC'=>1, # (Ser/S) Serine
                      'CC'=>1, # (Pro/P) Proline
                      'AC'=>1, # (Thr/T) Threonine
                      'GC'=>1, # (Ala/A) Alanine
                      'CG'=>1, # (Arg/R) Arginine
                      'GG'=>1, # (Gly/G) Glycine
                    );
    if ($hend == $cds_loc->end+1) {
            my $tail = $cds->seq->seq;
            $tail = substr($tail, length($tail)-2, 2);
            $debug && print STDERR "$subname: \$tail=$tail\n";
            $hend-- if (exists($good_tail{$tail}) && $good_tail{$tail});
    }

    $debug && print STDERR "$subname: \$qend=$hend\n";
    return $hend;

} # sub adjust_end

=head2 msa_get_feature_loc

Takes reffeat, gaps in the reffeat, and gaps (in DNA coordinate) in the target
 calculate the location of the feature in target genome
 Return a Bio::Location object for the new annotation

Note: We can only deal with Bio::Location::Simple, not more complicated location types

=cut

sub msa_get_feature_loc {
    my ($reffeat_loc, $refcds_loc, $cds, $aln, $aln_q, $aln_h) = @_;

    my $cds_loc = $cds->location;
    my $debug = 0 && $debug_all;
    my $subname = 'msa_get_feature_loc';

    my $loc2 = undef;
    my $errcode = {};

    $debug && print STDERR "$subname:\n";
    if ($reffeat_loc->isa('Bio::Location::Simple') ) {
        my $qstart = $reffeat_loc->start;
        my $qend = $reffeat_loc->end;
        $debug && print STDERR "$subname: \$qstart=$qstart \$qend=$qend\n";

        $qstart = Annotate_Util::adjust_start( $qstart, $refcds_loc, $cds, $aln, $aln_q, $aln_h);
        $qend = Annotate_Util::adjust_end( $qend, $refcds_loc, $cds, $aln, $aln_q, $aln_h);
        my $codon_start = 1;
        if ($cds->has_tag('codon_start')) {
            $codon_start = [$cds->get_tag_values('codon_start')];
            $codon_start = $codon_start->[0];
        }

        $qstart = $qstart + $codon_start -1;
        $qend = $qend + $codon_start -1;

        $debug && print STDERR "$subname: \$qstart=$qstart \$qend=$qend\n";

        if ($qstart<=$qend) {
            $loc2 = Bio::Location::Simple->new();
            $loc2->start($qstart);
            $loc2->end  ($qend);
        }

    } elsif ($reffeat_loc->isa('Bio::Location::Fuzzy')) {


        my $codon_start = 1;
        if ($cds->has_tag('codon_start')) {
            $codon_start = [$cds->get_tag_values('codon_start')];
            $codon_start = $codon_start->[0];
        }

        my $start = 9999999;
        my $end   = -1;
        my $qstart_min = $reffeat_loc->min_start;
        if (defined($qstart_min)) {
            $debug && print STDERR "$subname: \$qstart_min=$qstart_min\n";
            $qstart_min = Annotate_Util::adjust_start( $qstart_min, $refcds_loc, $cds, $aln, $aln_q, $aln_h);
            $qstart_min = $qstart_min + $codon_start -1;
            $start = $qstart_min;
        }

        my $qend_min = $reffeat_loc->min_end;
        if (defined($qend_min)) {
            $debug && print STDERR "$subname: \$qend_min=$qend_min\n";
            $qend_min = Annotate_Util::adjust_end( $qend_min, $refcds_loc, $cds, $aln, $aln_q, $aln_h);
            $qend_min = $qend_min + $codon_start -1;
            $end = $qend_min;
        }

        my $qstart_max = $reffeat_loc->max_start;
        if (defined($qstart_max)) {
            $qstart_max = Annotate_Util::adjust_start( $qstart_max, $refcds_loc, $cds, $aln, $aln_q, $aln_h);
            $qstart_max = $qstart_max + $codon_start -1;
            $start = $qstart_max if (not $start<=$qstart_max);
        }

        my $qend_max = $reffeat_loc->max_end;
        if (defined($qend_max)) {
            $qend_max = Annotate_Util::adjust_end( $qend_max, $refcds_loc, $cds, $aln, $aln_q, $aln_h);
            $qend_max = $qend_max + $codon_start -1;
            $end = $qend_max if ($end<=$qend_max);
        }

        if ($start<=$end) {
            $loc2 = Bio::Location::Fuzzy->new();
            $loc2->location_type($reffeat_loc->location_type);
            $loc2->start_pos_type($reffeat_loc->start_pos_type);
            $loc2->end_pos_type  ($reffeat_loc->end_pos_type  );
            $loc2->min_start($qstart_min);
            $loc2->min_end  ($qend_min);
            $loc2->max_start($qstart_max);
            $loc2->max_end  ($qend_max);
        }
    }
    # check the validity of $loc2, this ensures the new annotated feature lies within the cds,
    # and doesn't have more gap than 3 dna residues
    $errcode = &msa_check_loc($loc2, $reffeat_loc, $refcds_loc, $cds_loc, $aln, $aln_q, $aln_h);

    $debug && print STDERR "$subname: \$loc2=\n".Dumper($loc2)."end of \$loc2\n";
    $debug && print STDERR "$subname: \$errcode=\n".Dumper($errcode)."end of \$errcode\n";

    return ($loc2, $errcode);

} # sub msa_get_feature_loc


=head2 msa_check_loc

Takes reffeat, gaps in the reffeat, and gaps (in DNA coordinate) in the target
 calculate the location of the feature in target genome
 Return a Bio::Location object for the new annotation

=cut

sub msa_check_loc {
    my ($loc2, $reffeat_loc, $refcds_loc, $cds_loc, $aln, $aln_q, $aln_h) = @_;

    my $debug = 0 && $debug_all;
    my $subname = 'msa_check_loc';
    # error code: OutsideCDS; partial_mat_peptide; GapAtCleaveSite
    my $errcode = {OutsideCDS=>0, partial_mat_peptide=>0};

    return ($errcode) if (!$loc2);

    my $msg;
    my $max_gap_end = 3;
    my $max_gap_middle = 3;

    # See if $loc2 is outside of $cds_loc
    if ($loc2->start < $cds_loc->start) {  # if $loc2 is outside of $cds_loc
            print STDERR "$subname: a ".$loc2->start.' '.$cds_loc->start."\n";
    }
    if ($loc2->end <= $cds_loc->start) {  # if $loc2 is outside of $cds_loc
            print STDERR "$subname: b ".$loc2->end.' '.$cds_loc->start."\n";
    }
    if ($loc2->start >= $cds_loc->end) {  # if $loc2 is outside of $cds_loc
            print STDERR "$subname: c ".$loc2->start.' '.$cds_loc->end."\n";
    }
    if ($loc2->end > $cds_loc->end) {  # if $loc2 is outside of $cds_loc
            print STDERR "$subname: d ".$loc2->end.' '.$cds_loc->end."\n";
    }
    if ($loc2->start < $cds_loc->start
         || $loc2->end <= $cds_loc->start
         || $loc2->start >= $cds_loc->end
         || $loc2->end > $cds_loc->end
       ) {  # if $loc2 is outside of $cds_loc
            print STDERR "$subname: e \$loc2=".$loc2->to_FTstring." doesn't match \$cds=".$cds_loc->to_FTstring."\n";
    # error code: OutsideCDS; partial_mat_peptide; GapAtCleaveSite
            $errcode->{OutsideCDS} = 1;
            return ($errcode);
    }

    $errcode->{partial_mat_peptide} = &msa_check_loc_start($loc2, $reffeat_loc, $refcds_loc, $cds_loc, $aln, $aln_q, $aln_h);

    $errcode->{partial_mat_peptide} = $errcode->{partial_mat_peptide} || &msa_check_loc_end($loc2, $reffeat_loc, $refcds_loc, $cds_loc, $aln, $aln_q, $aln_h);

    $errcode->{long_internal_gap} = &msa_check_loc_internal($loc2, $reffeat_loc, $refcds_loc, $cds_loc, $aln, $aln_q, $aln_h);

    $debug && print STDERR "$subname: \$errcode=\n".Dumper($errcode)."\n";

    return ($errcode);

} # sub msa_check_loc

sub msa_check_loc_internal {
    my ($loc2, $reffeat_loc, $refcds_loc, $cds_loc, $aln, $aln_q, $aln_h) = @_;

    my $debug = 0 && $debug_all;
    my $subname = 'msa_check_loc_internal';
    my $errcode = { long_internal_gap => 0 };

        my $gap_char = $aln->gap_char;
        my $qstart = $reffeat_loc->start; # get start of reffeat
        $qstart = ($qstart - $refcds_loc->start)/3 +1; # convert to AA coordinate in CDS
        my $qseq = $aln_q->seq;
#        $debug && print STDERR "$subname: \$qseq='$qseq' \$qseq=".length($qseq)."\n";
        my $char;
        for (my $n=1; $n<=$qstart && $n<=length($qseq); $n++) {
            $char = substr($qseq, $n-1, 1);
            $qstart++ if ($char eq $gap_char);
#            $debug && print STDERR "$subname: \$n=$n \$char=$char \$gap_char=$gap_char \$qstart=$qstart\n";
        }

        my $qend = $reffeat_loc->end; # get end of reffeat
        $qend = ($qend - $refcds_loc->start +1)/3; # convert to AA coordinate in CDS
        $qseq = $aln_q->seq;
#        $debug && print STDERR "$subname: \$qseq='$qseq' \$qseq=".length($qseq)."\n";
        for (my $n=1; $n<=$qend && $n<=length($qseq); $n++) {
            $char = substr($qseq, $n-1, 1);
            $qend++ if ($char eq $gap_char);
#            $debug && print STDERR "$subname: \$n=$n \$char=$char \$gap_char=$gap_char \$qend=$qend\n";
        }

        # Count the number of gaps in the hit sequence
        my $hseq = $aln_h->seq;
        $debug && print STDERR "$subname: \$hseq='$hseq' \$hseq=".length($hseq)."\n";
        $hseq = substr($hseq, $qstart, $qend-$qstart+1); # Any gap at the beginning of CDS, plus the 1st AA of current feat.
        $debug && print STDERR "$subname: \$hseq='$hseq' \$hseq=".length($hseq)."\n";
        $hseq =~ s/^[$gap_char]*//; # remove any gap at start
        $hseq =~ s/[$gap_char]*$//; # remove any gap at end
        my $count = 0;
        my @chars = split(//, $hseq);
        foreach my $c (@chars) {
            $count++ if ($c eq $gap_char);
        }
        $debug && print STDERR "$subname: \$hseq='$hseq' \$hseq=".length($hseq)." \$count=$count\n";
        if ($count > 0.333*length($hseq)) {
            print STDERR "$subname: ERROR: \$qstart=$qstart \$qend=$qend \$gap_char=$gap_char\n";
            print STDERR "$subname: ERROR: \$hseq='$hseq' \$hseq=".length($hseq)." \$count=$count\n";
            $errcode->{long_internal_gap} = 1;
        }
    $debug && print STDERR "$subname: \$errcode=\n".Dumper($errcode)."\n";

    return ($errcode->{long_internal_gap});

} # sub msa_check_loc_internal


sub msa_check_loc_start {
    my ($loc2, $reffeat_loc, $refcds_loc, $cds_loc, $aln, $aln_q, $aln_h) = @_;

    my $debug = 0 && $debug_all;
    my $subname = 'msa_check_loc_start';
    # error code: OutsideCDS; partial_mat_peptide; GapAtCleaveSite
    my $errcode = { partial_mat_peptide => 0 };

        my $qstart = $reffeat_loc->start; # get start of reffeat
        $qstart = ($qstart - $refcds_loc->start)/3 +1; # convert to AA coordinate in CDS
        my $qseq = $aln_q->seq;
        $debug && print STDERR "$subname: \$qseq='$qseq' \$qseq=".length($qseq)."\n";
        my $gap_char = $aln->gap_char;
        my $char;
        for (my $n=1; $n<=$qstart && $n<=length($qseq); $n++) {
            $char = substr($qseq, $n-1, 1);
            $qstart++ if ($char eq $gap_char);
            $debug && print STDERR "$subname: \$n=$n \$char=$char \$gap_char=$gap_char \$qstart=$qstart\n";
        }

        # Any gap in hseq at the cleavage site makes the mat_peptide partial
        my $hseq = $aln_h->seq;
        $debug && print STDERR "$subname: \$hseq='$hseq' \$hseq=".length($hseq)."\n";
        $char = substr($hseq, 0, $qstart); # Any gap at the beginning of CDS, plus the 1st AA of current feat.
        $debug && print STDERR "$subname: \$char='$char' \$char=".length($char)."\n";
        if ($char =~ /[$gap_char]+$/) {
            $debug && print STDERR "$subname: \$qstart=$qstart \$char=$char \$gap_char=$gap_char\n";
            $errcode->{partial_mat_peptide} = 1;
        }
    $debug && print STDERR "$subname: \$errcode=\n".Dumper($errcode)."\n";

    return ($errcode->{partial_mat_peptide});

} # sub msa_check_loc_start

sub msa_check_loc_end {
    my ($loc2, $reffeat_loc, $refcds_loc, $cds_loc, $aln, $aln_q, $aln_h) = @_;

    my $debug = 0 && $debug_all;
    my $subname = 'msa_check_loc_end';
    # error code: OutsideCDS; partial_mat_peptide; GapAtCleaveSite
    my $errcode = { partial_mat_peptide => 0 };

        my $qend = $reffeat_loc->end; # get end of reffeat
        $qend = ($qend - $refcds_loc->start +1)/3; # convert to AA coordinate in CDS
        my $qseq = $aln_q->seq;
        $debug && print STDERR "$subname: \$qseq='$qseq' \$qseq=".length($qseq)."\n";
        my $gap_char = $aln->gap_char;
        my $char;
        for (my $n=1; $n<=$qend && $n<=length($qseq); $n++) {
            $char = substr($qseq, $n-1, 1);
            $qend++ if ($char eq $gap_char);
            $debug && print STDERR "$subname: \$n=$n \$char=$char \$gap_char=$gap_char \$qend=$qend\n";
        }
        # Any gap in hseq at the cleavage site makes the mat_peptide partial
        my $hseq = $aln_h->seq;
        $debug && print STDERR "$subname: \$hseq='$hseq' \$hseq=".length($hseq)."\n";
        $char = substr($hseq, $qend-1, length($hseq)-$qend+1); # Any gap at the end of CDS, plus the last AA of current feat.
        $debug && print STDERR "$subname: \$char='$char' \$char=".length($char)."\n";
        if ($char =~ /^[$gap_char]+/) {
            $debug && print STDERR "$subname: \$qend=$qend \$char=$char \$gap_char=$gap_char\n";
            $errcode->{partial_mat_peptide} = 1;
        }

    $debug && print STDERR "$subname: \$errcode=\n".Dumper($errcode)."\n";

    return ($errcode->{partial_mat_peptide});

} # sub msa_check_loc_end


=head2 msa_check_loc_fuzzy

Takes reffeat, gaps in the reffeat, and gaps (in DNA coordinate) in the target
 calculate the location of the feature in target genome
 Return a Bio::Location object for the new annotation

Note: We can only deal with Bio::Location::Simple, not more complicated location types

=cut
=head1
sub msa_check_loc_fuzzy {
    my ($loc2, $reffeat, $refcds_loc, $cds_loc, $aln, $aln_q, $aln_h) = @_;

    my $debug = 1 && $debug_all;
    my $subname = 'msa_check_loc_fuzzy';
    my $errcode = {};
    my $msg;
    my $max_gap_end = 3;
    my $max_gap_middle = 3;

    if ($loc2->start < $cds_loc->start
       ) {  # if $feat starts before CDS starts
            print STDERR "$subname: 1 ".$loc2->start.' '.$cds_loc->start."\n";
    }
    if ($loc2->end <= $cds_loc->start
       ) {  # if $feat end before CDS starts
            print STDERR "$subname: 2 ".$loc2->end.' '.$cds_loc->start."\n";
    }
    if ($loc2->start >= $cds_loc->end
       ) {  # if $feat starts after CDS ends
            print STDERR "$subname: 3 ".$loc2->start.' '.$cds_loc->end."\n";
    }
    if ($loc2->end > $cds_loc->end
       ) {  # if $feat ends after CDS ends
            print STDERR "$subname: 4 ".$loc2->end.' '.$cds_loc->end."\n";
    }
    if ($loc2->start < $cds_loc->start
         || $loc2->end <= $cds_loc->start
         || $loc2->start >= $cds_loc->end
         || $loc2->end > $cds_loc->end
       ) {  # if $feat is outside of CDS
            print STDERR "$subname: \$loc2=".$loc2->to_FTstring." doesn't match \$cds=".$cds_loc->to_FTstring."\n";
    # error code: 0: none; 1: gaps around cleavage site; 2: partial mat_peptide; 3: outside of CDS
            $errcode->{OutsideCDS} = 1;
            return ($errcode);
    }

    $debug && print STDERR "$subname: \n";


    return ($errcode);

} # sub msa_check_loc_fuzzy
=cut

=head2 adjust_from_gaps

Takes a position, and an array of gaps, compare if any gap is before the position.
 First array element is ignored.
 Return the number of gaps before the start

=cut

sub adjust_from_gaps {
        my ($position, $gaps_q) = @_;

        my $debug = 0 && $debug_all;
        my $adjust; # number of gaps
        $adjust = 0;
        foreach my $i (1 .. $#{@$gaps_q}) {
            $debug && print STDERR "adjust_from_gaps: \$position=$position \t\$i=$i \$gap=$gaps_q->[$i]\n";
            if ($gaps_q->[$i] < ($position + $adjust)) { # This $position doesn't have any gap
                $adjust = $i;
            } else {
                last;
            }
        }

        return ($adjust);
} # sub adjust_from_gaps


=head2 msa_get_gaps

Takes an alignment, and an id, return the gaps within the alignment with such id

=cut

sub msa_get_gaps {
    my ($aln, $id) = @_;

    my $debug = 0 && $debug_all;
    my $subname = 'msa_get_gaps';

    my $seq = &msa_get_aln($aln, $id);
    $debug && print "$subname: \$refaln = ".$seq->display_id()."\n";
    $debug && print "$subname: \$refaln = ".$seq->start()."\n";
    $debug && print "$subname: \$refaln = ".$seq->end()."\n";
    $debug && print "$subname: \$refaln = ".$seq->alphabet()."\n";
    $debug && print "$subname: \$refaln = ".$seq->length."\n";
    $debug && print "$subname: \$refaln = ".$seq->seq."\n";

    my @gaps = (''); # we don't need element 0, but need it to suppress error from print
    for (my $i = 1; $i<= $seq->length; $i++) {
        my $j;
        $j = ($i==1) ? 1+ index($seq->seq,$aln->gap_char) : 1+ index($seq->seq,$aln->gap_char,$gaps[$#gaps]);
        if ( $j>0 ) {
            push @gaps, $j;
            $i = $j;
        } elsif ( $j==0 ) {
            $i = $seq->length;
            last;
        }
#        $gaps[1] && print "annotate_msa: i=$i j=$j\n";
    }
    $debug && print "$subname: \@gaps($#gaps)='@gaps'\n";

    return (\@gaps);
} # sub msa_get_gaps


=head2 msa_get_aln

Takes an alignment, and an id, return the gaps within the alignment with such id

=cut

sub msa_get_aln {
    my ($aln, $id) = @_;

    my $debug = 0 && $debug_all;
    my $subname = 'msa_get_aln';
    my $refaln = [ $aln->each_seq_with_id($id) ]; # 
    $debug && print "$subname: \$refaln = \n".Dumper($refaln)."End of \$refaln\n\n";
    if ($#{@$refaln} < 0) {
        $debug && print "$subname: Couldn't find id=$id in alignment file.\n";
        return undef;
    }
    my $seq = $refaln->[0];
    $debug && print "$subname: \$refaln = ".$seq->display_id()."\n";
    $debug && print "$subname: \$refaln = ".$seq->start()."\n";
    $debug && print "$subname: \$refaln = ".$seq->end()."\n";
    $debug && print "$subname: \$refaln = ".$seq->alphabet()."\n";
    $debug && print "$subname: \$refaln = ".$seq->length."\n";
    $debug && print "$subname: \$refaln = ".$seq->seq."\n";

    my @gaps = (''); # we don't need element 0, but need it to suppress error from print
    $seq = $refaln->[0];
    $debug && print "$subname: \$seq='$seq'\n";

    return ($seq);
} # sub msa_get_aln


=head2 get_refpolyprots

Takes a sequence object,
 return the [CDS_id, CDS, mat_peptide 1, mat_peptide 2, ...] in the refseq. Empty if there is any problem

=cut

sub get_refpolyprots {
    my ( $refseqs, $inseq, $exe_dir) = @_;

    my $debug = 0 && $debug_all;
    my $subname = 'get_refpolyprots';

    my $refpolyprots = [];
    my $refseq;

    my @ann = $inseq->annotation->get_Annotations('date_changed');
    $debug && print STDERR "$subname: accession=".$inseq->accession_number."\tdate=$ann[0]->{'value'}\n";

#    $debug && print STDERR "$subname: \$debug=$debug\n";
#    $debug && print STDERR "$subname: \$refseq=".ref($refseq)."\n";
#    $debug && print STDERR "$subname: \$refseq->isa=".$refseq->isa('Bio::Seq::RichSeq')."\n" if ($refseq);
#    $debug && print STDERR "$subname: \$refseq_fit=".Annotate_Util::check_refseq($refseq, $inseq)."\n";
#    $debug && print STDERR "$subname: \$refseq_required=$refseq_required\n";
    if (!$refseq || !$refseq->isa('Bio::Seq::RichSeq') || !check_refseq($refseq, $inseq)) {

        $refseq = Annotate_Util::get_refseq( $inseq, $exe_dir);
        if (!$refseq) {
            my $acc = $inseq->accession;
            my $id  = $inseq->species->ncbi_taxid;
            print STDERR "$subname: ERROR genome=".$acc." w/ taxid=".$id." is not covered in V$VERSION.\n";
            print STDOUT "$subname: ERROR genome=".$acc." w/ taxid=".$id." is not covered in V$VERSION.\n";
            print STDERR "$subname: Please contact script author for any update.\n";
            return undef;
        } else {
#            print "$subname: \$refseq=$refseq\n\n";
            my $acc = $refseq->accession_number;
            my $id  = $refseq->species->ncbi_taxid;
            $debug && print STDOUT "$subname: Refseq changed to \$refseq=".$acc.' w/ taxid='.$id."\n";
            $debug && print STDERR "$subname: Refseq changed to \$refseq=".$acc.' w/ taxid='.$id."\n";
        }
    }
    if (exists($refseqs->{$refseq->accession_number})) {
        return $refseqs->{$refseq->accession_number};
    }

    # going through all features in refseq, look for polyproptein CDS
    my $reffeats = [ $refseq->get_SeqFeatures ];

    for (my $refct = 0; $refct<=$#{@$reffeats}; $refct++) {
        my $reffeat = $reffeats->[$refct];
        $debug && print STDERR "$subname: \$refct=$refct \$reffeat=".$reffeat->primary_tag."\n";
        next if ($reffeat->primary_tag ne 'CDS');

        my $is_poly = 0;
        $is_poly = &is_polyprotein($reffeat, ['product', 'note']);
        if (!$is_poly) {
            $debug && print STDERR "$subname: \$refct=$refct \$is_poly=$is_poly is not a polyprotein, skipping\n";
            next; # only run for CDS labeled as polyprotein in refseq
        } else {
            $debug && print STDERR "$subname: \$refct=$refct \$reffeat=".$reffeat->primary_tag."\n";
        }

        my $refcds = $reffeat;
        $debug && print STDERR "$subname: Found polyprotein at \$refct=$refct\n";

        # Found a polyprotein CDS, now collect all mat_peptide and sig_peptide following it
        my $refmatps = [];
        my %allowed_tags = ( 'mat_peptide' => 1, 'sig_peptide' => 1, );
        while ($refct<$#{@$reffeats} && $allowed_tags{$reffeats->[$refct+1]->primary_tag}) {
            $debug && print STDERR "$subname: \$refct=$refct \$reffeat=".$reffeats->[$refct+1]->primary_tag.' '.$reffeats->[$refct+1]->location->to_FTstring."\n";
            push @$refmatps, $reffeats->[++$refct];
        }

        my $cds_id = '';
        my @allowed_tags = ( 'db_xref', 'protein_id');
        foreach my $tag (@allowed_tags) {
            next if (!$refcds->has_tag($tag));
            my @values = $refcds->get_tag_values($tag);
            foreach my $v (@values) {
                if ($v =~ /^(GI:.+)$/i) {
                    $cds_id = $1;
                } elsif ($v =~ /^(NP_\d+)[.]\d+$/i) {
                    $cds_id = $1;
                }
                $debug && $cds_id && print STDERR "$subname: \$cds_id=$cds_id\n";
                last if ($cds_id);
            }
            last if ($cds_id);
        }
        $cds_id = 'unknown' if (!$cds_id);

        push @$refpolyprots, [
                               $cds_id,
                               $refcds,
                               @$refmatps,
                             ];
        $debug && print STDERR "$subname: \$#refpolyprots=$#{@$refpolyprots}\n";
        $debug && print STDERR "$subname: \$refpolyprots=\n".Dumper($refpolyprots)."end of \$refpolyprots\n\n";
    }

    return $refpolyprots;
} # sub get_refpolyprots


sub is_polyprotein {
    my ($feat, $allowed_tags) = @_;

    my $debug = 0 && $debug_all;

    my $is_poly = 0;
    foreach my $tag (@$allowed_tags) {
#        $debug && print STDERR "is_polyprotein: \$refct=$refct \$tag=$tag\n";
        next if (!$feat->has_tag($tag));
#        $debug && print STDERR "is_polyprotein: \$refct=$refct Found \$tag=$tag\n";
        my @tag = $feat->get_tag_values($tag);
        foreach my $t (@tag) {
#            $debug && print STDERR "is_polyprotein: \$refct=$refct \$t=$t\n";
            if ($t =~ /polyprotein/i) {
                $is_poly = 1;
                last;
            }
        }
    }

    return $is_poly;
} # sub is_polyprotein


=head2 get_polyprots

Takes an inseq object, and array of [CDS, mat_peptide 1, mat_peptide 2, ...]
 return the {CDS_id => [[CDS, mat_peptide 1, mat_peptide 2, ...] [] ...] in the sequence and refseq. Empty if there is any problem

=cut

sub get_polyprots {
    my ($inseq, $refpolyprots) = @_;

    my $debug = 0 && $debug_all;
    my $subname = 'get_polyprots';
    my $polyprots = {};
    if ($#{@$refpolyprots} < 0) {
            return $polyprots;
    }

    my $acc = $inseq->accession_number;
    my $num_cds = -1;
    for (my $i = 0; $i<=$#{@$refpolyprots}; $i++) {
        my $refset = $refpolyprots->[$i];
        my $reffeat = $refset->[1]; # [0] is CDS_id
        $debug && print STDERR "$subname: \$refseq=".$reffeat->seq->accession_number."   \$inseq=".$inseq->accession_number."\n";
        if ($reffeat->primary_tag ne 'CDS') {
            $debug && print STDERR "$subname: \$i=$i \$reffeat=".$reffeat->primary_tag." is not CDS as expected\n";
            $debug && print STDERR "$subname: \$i=$i Skipped.\n";
            next;
        }
        my $refcds = $reffeat;
        my $feats = [ $inseq->get_SeqFeatures ];
        for (my $ct = 0; $ct<=$#{@$feats}; $ct++) {
            my $feat = $feats->[$ct];
#            $debug && print STDERR "$subname: \$ct=$ct \$feat=".$feat->primary_tag."\n";
            if ($feat->primary_tag ne 'CDS') {
                # Skip those features such as source, gene, 5'URT
                $debug && print STDERR "$subname: \$ct=$ct \$feat=".$feat->primary_tag." is not CDS\n";
                next;
            } else {
                $debug && print STDERR "$subname: \$ct=$ct \$feat=".$feat->primary_tag."\n";
            }

            my $cds = $feat;
            $debug && print STDERR "$subname: Found polyprotein at \$ct=$ct\n";

            # double check to ensure the translation of CDS is right. For instance, DQ430819
            my $tr2 = &get_new_translation($cds, $cds);
            if (!$tr2) {
                print STDERR "$subname: get_new_translation returned translation='$tr2'. Skip.\n";
                next; # Skip if no translation can be obtained
            }
            $debug && print STDERR "$subname: translation='$tr2'\n";

            my $match_refcds = 0;
            $match_refcds = &match_cds_bl2seq($cds, $refcds);
            if (!$match_refcds) {
                print STDERR "$subname: \$ct=$ct \$feat=".$feat->primary_tag." doesn't match refcds\n";
                next;
            }

            # Found a polyprotein CDS, now collect all mat_peptide and sig_peptide following it
            my $matps = [];
            my %allowed_tags = ( 'mat_peptide' => 1, 'sig_peptide' => 1, 'misc_feature' => 1,);
            while ($ct<$#{@$feats} && $allowed_tags{$feats->[$ct+1]->primary_tag}) {
                ++$ct;
                next if ($feats->[$ct]->primary_tag eq 'misc_feature');
                $debug && print STDERR "$subname: \$ct=$ct \$feat=".$feats->[$ct]->primary_tag.' '.$feats->[$ct]->location->to_FTstring."\n";
                push @$matps, $feats->[$ct];
            }

            $debug && print STDERR "$subname: \$#polyprots=\n".Dumper($polyprots)."End of \%\$polyprots\n";

            # check if this is a new CDS, e.g. M55506 has a 2nd CDS that's part of 1st one
            my $seen=0;
            foreach my $j (keys %$polyprots) {
              foreach my $k (0 .. $#{%{$polyprots->{$j}}}) {
                my $old_cds = $polyprots->{$j}->[$k]->[0];
                $debug && print STDERR "$subname: \$j=$j \$k=$k \$old_cds=\n".Dumper($old_cds)."end of \$old_cds\n\n";
                $debug && print STDERR "$subname: \$j=$j \$k=$k \$cds=\n".Dumper($cds)."end of \$cds\n\n";
                my $s1 = &get_new_translation($old_cds, $old_cds);
                my $s2 = &get_new_translation($cds, $cds);
                $debug && print STDERR "$subname: \$s1=$s1\n";
                $debug && print STDERR "$subname: \$s2=$s2\n";
                if ($s1 =~ /$s2/i) { # $s1 is same to or longer than $s2
                    $debug && print STDERR "$subname: \$s1 is same to or longer than \$s2\n";
                    $seen=1;
                    push @{$polyprots->{$j}->[$k]}, @$matps;
                    last;
                } elsif ($s2 =~ /$s1/i) {
                    # if $s1 is shorter than $s2, need to update $old_cds, and append any mat_peptide to list
                    $debug && print STDERR "$subname: \$s2 is same to or longer than \$s1\n";
                    $seen=1;
                    $old_cds = $cds;
                    push @{$polyprots->{$j}->[$k]}, @$matps;
                    last;
                } else {
                    # if $s1 is not related to $s2
                    $debug && print STDERR "$subname: \$s1 is not related to \$s2\n";
                }
              }
            }

            if (!$seen) {
                push @{$polyprots->{$refset->[0]}}, [$cds, @$matps]; # $refset->[0] is CDS_id of refseq
                $num_cds++;
            }

            $debug && print STDERR "$subname: \$#polyprots=". keys(%$polyprots)."End of \$#polyprots\n";
            $debug && print STDERR "$subname: \$polyprots=\n".Dumper($polyprots)."end of \$polyprots\n\n";

        }

        print STDERR "$subname: \$refseq=".$reffeat->seq->accession_number."   \$inseq=".$acc;
        print STDERR "  \$num_cds=$num_cds\n";

    } #     for (my $i = 0; $i<=$#{@$refpolyprots}; $i++)

    my $n = [keys %$polyprots];
    $debug && print STDERR "$subname: polyprotein CDS for acc=$acc \$n=$#{@$n}\n";
    $debug && print STDERR "$subname: \$polyprots=\n".Dumper($polyprots)."end of \$polyprots\n\n";
    if ($#{@$n} <0) {
        print STDERR "$subname: ERROR: Can't find any polyprotein CDS for acc=$acc \$n=$#{@$n}.\n";
        print STDOUT "$subname: ERROR: Can't find any polyprotein CDS for acc=$acc \$n=$#{@$n}.\n";
#        return;
    }

    return ($polyprots, $num_cds);
} # sub get_polyprots


=head2 match_cds_bl2seq

Takes 2 CDS feature
 Compares the translation of 2 CDS by bl2seq, see if the conserved length is at least 80% of the length of the shorter seq
 Return 1 if the 2 seqs are similar
=cut

sub match_cds_bl2seq {
    my ($cds, $refcds) = @_;

    my $debug = 0 && $debug_all;
    my $subname = 'match_cds_bl2seq';

    my $match_cds = 0;

    my @values = $refcds->get_tag_values('translation');
    my $s1 = $values[0];

    @values = $cds->get_tag_values('translation');
    my $s2 = $values[0];
    my $bl2seq_result = &run_bl2seq_search($s1, $s2);
    $s1 = length($s1);
    $s2 = length($s2);
    $debug && print STDERR "$subname: \$bl2seq_result=\n".Dumper($bl2seq_result)."end of \$bl2seq_result\n\n";
    $debug && print STDERR "$subname: \$s1=$s1 \$s2=$s2\n\n";

    # Look at the hit
    my $hit_cds = $bl2seq_result->next_hit;
    $debug && print STDERR "$subname: \$hit_cds=\n".Dumper($hit_cds)."end of \$hit_cds\n\n";
    while (my $hsp_cds = $hit_cds->next_hsp) {
        $debug && print STDERR "$subname: \$hsp_cds=\n".Dumper($hsp_cds)."end of \$hsp_cds\n\n";

        $s1 = $s1<$s2 ? $s1 : $s2;
        $debug && print STDERR "$subname: \$s1=$s1 length of conserved=".$hsp_cds->{'CONSERVED'}."\n";
        my $conserved_residues_required = 0.66667;
        $conserved_residues_required = 0.5 if ($debug);
        if ($hsp_cds->{'CONSERVED'} >= $s1 * $conserved_residues_required) {
            $match_cds = 1;
            $debug && print STDERR "$subname: CDS length=$s1 conserved=".$hsp_cds->{'CONSERVED'}.".";
            $debug && print STDERR " Length of conserved meets required $conserved_residues_required\n";
            $debug && print STDERR "$subname: \$match_cds=$match_cds\n";
            last;
        } else {
            print STDERR "$subname: CDS length=$s1 conserved=".$hsp_cds->{'CONSERVED'}.".";
            print STDERR " Length of conserved doesn't meet required $conserved_residues_required\n";
            print STDERR "$subname: \$match_cds=$match_cds\n";
            print STDERR "$subname: QUERY_SEQ='".$hsp_cds->query_string."'\n";
            print STDERR "$subname:           '".$hsp_cds->homology_string."'\n";
            print STDERR "$subname:   HIT_SEQ='".$hsp_cds->hit_string."'\n";
        }
    }

    return $match_cds;
} # sub match_cds_bl2seq


=head2 get_refseq

Takes either a file name, or a Bio::Seq object based on genbank file, returns refseq in Bio::Seq.
 For a file name, simply load the genbank file
 For a Bio::Seq object, look at its taxon. Load the refseq either from file system (or MySQL)
 Return the refseq in Bio::Seq object.
=cut

sub get_refseq {
    my ($inseq, $exe_dir) = @_;

    my $debug = 0 && $debug_all;
    my $refseq = undef;
    return $refseq if (!$inseq);

    $debug && print STDERR "get_refseq: \$inseq=$inseq is a ".ref($inseq)."\n";
    my $refseq_fn;

    if ($inseq->isa('Bio::Seq::RichSeq')) {
            my $refacc = undef;
            $refacc = &get_refseq_acc($inseq);
            return undef if (!$refacc);
            $refseq_fn = $refacc.'_matpeptide.gb';
            $debug && print STDERR "get_refseq: \$exe_dir=${exe_dir}refseq/$refseq_fn\n";
            $debug && print STDERR "get_refseq: REFSEQ is $refseq_fn\n";
            if (!-e "${exe_dir}refseq/$refseq_fn") {
                $refseq_fn = $refacc.'.gb';
            }
        $debug && print STDERR "get_refseq: REFSEQ is $refseq_fn\n";

    } elsif ($inseq) {
        $refseq_fn = $inseq;
    }

    # Get refseq object if $refseq_fn is given, search in ./refseq, ./, and $dir_path/$ticket/
    # First look for refseq in file system
    $exe_dir = './' if (!$exe_dir);
    if ($refseq_fn) {
        $debug && print STDERR "get_refseq: REFSEQ is $refseq_fn\n";
        if (-e "$exe_dir/refseq/$refseq_fn") {
            $refseq = Bio::SeqIO->new( -file => "$exe_dir/refseq/$refseq_fn")->next_seq();
        } elsif (-e "$exe_dir/$refseq_fn") {
            $refseq = Bio::SeqIO->new( -file => "$exe_dir/$refseq_fn")->next_seq();
        } elsif (-e "./$refseq_fn") {
            $refseq = Bio::SeqIO->new( -file => "./$refseq_fn")->next_seq();
        } else {
            print STDERR "get_refseq: Can't find refseq gb file: $refseq_fn\n";
        }
    }

    return $refseq;
} # sub get_refseq


=head2 get_refseq_acc

Takes either a file name, or a Bio::Seq object based on genbank file, returns refseq in Bio::Seq.
 For a file name, simply load the genbank file
 For a Bio::Seq object, look at its taxon. Load the refseq either from file system (or MySQL)
 Return the refseq in Bio::Seq object.
=cut

sub get_refseq_acc {
    my ($inseq) = @_;

    my $debug = 0 && $debug_all;
    my $subname = 'get_refseq_acc';

    my $refseq_acc = undef;
    return $refseq_acc if (!$inseq);

        # list of refseqs
        my $refseq_list = {
           # Flavivirus
##           11079 => 'NC_000943', # Murray Valley encephalitis virus
           11072 => 'NC_001437', # Japanese encephalitis virus
##           11073 => 'NC_001437', # Japanese encephalitis virus. Refseq has 13 mat_peptides, w/ fuzzy ranges
##           11076 => 'NC_001437', # Japanese encephalitis virus. Refseq has 13 mat_peptides, w/ fuzzy ranges
           11060 => 'NC_001474', # Dengue virus 2
           11069 => 'NC_001475', # Dengue virus 3
           11053 => 'NC_001477', # Dengue virus 1
           12637 => 'NC_001477', # Dengue virus
           11082 => 'NC_001563', # West Nile virus (lineage II strain 956)
##           31658 => 'NC_001564', # Cell fusing agent virus
           11084 => 'NC_001672', # Tick-borne encephalitis virus. Refseq has 13 mat_peptides, w/ fuzzy ranges
##           11092 => 'NC_001672', # Tick-borne encephalitis virus. Refseq has 13 mat_peptides, w/ fuzzy ranges
##           11094 => 'NC_001672', # Tick-borne encephalitis virus. Refseq has 13 mat_peptides, w/ fuzzy ranges
##           47300 => 'NC_001672', # Tick-borne encephalitis virus. Refseq has 13 mat_peptides, w/ fuzzy ranges
##          638787 => 'NC_001672', # Tick-borne encephalitis virus. Refseq has 13 mat_peptides, w/ fuzzy ranges
##           11086 => 'NC_001809', # Louping ill virus
           11089 => 'NC_002031', # Yellow fever virus (YFV). This refseq has 13 mat_peptides, w/ fuzzy ranges
           11070 => 'NC_002640', # Dengue virus 4
##           64300 => 'NC_003635', # Modoc virus
##           64285 => 'NC_003675', # Rio Bravo virus
##           64280 => 'NC_003676', # Apoi virus
##           11083 => 'NC_003687', # Powassan virus
##           11085 => 'NC_003690', # Langat virus
##          161675 => 'NC_003996', # Tamana bat virus
           11103 => 'NC_004102', # Hepatitis C virus genotype 1
           31646 => 'NC_004102', # Hepatitis C virus genotype 1a
           31647 => 'NC_004102', # Hepatitis C virus genotype 1b
           31649 => 'NC_004102', # Hepatitis C virus genotype 2a
           31650 => 'NC_004102', # Hepatitis C virus genotype 2b
           31655 => 'NC_004102', # Hepatitis C virus genotype 6a
          356426 => 'NC_004102', # Hepatitis C virus genotype 3a
##           64312 => 'NC_004119', # Montana myotis leukoencephalitis virus
##          172148 => 'NC_004355', # Alkhurma hemorrhagic fever virus
##           64294 => 'NC_005039', # Yokose virus
##           12542 => 'NC_005062', # Omsk hemorrhagic fever virus
##          218849 => 'NC_005064', # Kamiti River virus
##           64286 => 'NC_006551', # Usutu virus
##           64287 => 'NC_006947', # Karshi virus
           11080 => 'NC_007580', # St. Louis encephalitis virus
##           11080 => 'NC_001563', # St. Louis encephalitis => West Nile virus (lineage II strain 956)
##          390844 => 'NC_008604', # Culex flavivirus
##           64283 => 'NC_008718', # Entebbe bat virus
##           44026 => 'NC_008719', # Sepik virus
##           64303 => 'NC_009026', # Bussuquara virus
##           59563 => 'NC_009028', # Ilheus virus
##           44024 => 'NC_009029', # Kokobera virus
#           40271 => 'NC_009823', # Hepatitis C virus genotype 2
#          356114 => 'NC_009824', # Hepatitis C virus genotype 3
#           33745 => 'NC_009825', # Hepatitis C virus genotype 4
#           33746 => 'NC_009826', # Hepatitis C virus genotype 5
#           42182 => 'NC_009827', # Hepatitis C virus genotype 6
#           11082 => 'NC_009942', # West Nile virus (lineage I strain NY99), missing 2k

           # SARS coronavirus
##          694009 => 'NC_004718', # SARS coronavirus
           # Caliciviridae
           11983 => 'NC_001959', # Norwalk virus
           95340 => 'NC_001959', # Norwalk virus
#          150080 => 'NC_001959', # Norwalk virus
           };
        $debug && print STDERR "$subname: \$refseq_list=\n".Dumper($refseq_list)."end of \$refseq_list\n\n";

        $debug && print STDERR "$subname: Taxon=$taxon->{loaded}\n";
        if (!$taxon->{loaded} && -f $taxon->{fn} ) {
            $debug && print STDERR "$subname: found taxon_fn=$taxon->{fn}\n";

            open my $taxon_file, '<', $taxon->{fn}
               or croak("$subname: Found $taxon->{fn}, but couldn't open: $OS_ERROR");
            while (<$taxon_file>) {
                s/(^[|]*\s+|\s+$)//x;
                my $words = [split(/\s*\|\s*/x)];
                $debug && print STDERR "$subname: \$words='@$words'\n";
                next if (!$words->[0] || $words->[0] !~ /^\d+$/x);
                next if (!$words->[1] || $words->[1] !~ /^\d+$/x); # ignore species=-1
                my ($strainid, $speciesid) = @{$words}[0..1];
                $debug && print STDERR "$subname: \$strainid=$strainid \$speciesid=$speciesid\n";
                if (exists($refseq_list->{$strainid})) {
                    $taxon->{speciesid}->{$speciesid}->{$strainid} = $refseq_list->{$strainid};
                } else {
                    $taxon->{speciesid}->{$speciesid}->{$strainid} = 1;
                }
            }
            close $taxon_file or croak "$subname: Couldn't close $taxon->{fn}: $OS_ERROR";
            my @keys = keys %{$taxon->{speciesid}};
            $taxon->{loaded} = 1 if ($#keys>1);

            $debug && print STDERR "$subname: finished reading list file: $taxon->{fn}.\n";
            $debug && print STDERR "$subname: \$taxon = \n".Dumper($taxon)."End of \$taxon\n\n";


        }
        $debug && print STDERR "$subname: \$taxon->{loaded}=$taxon->{loaded}\n";

        my $taxid = $inseq->species->ncbi_taxid;
        $debug && print "$subname: \$taxid=$taxid \$inseq=$inseq is a ".ref($inseq)."\n";


        # Determine the refseq by a 2-step process
        if (exists($refseq_list->{$taxid})) {
            $refseq_acc = $refseq_list->{$taxid};
        } elsif($taxon->{loaded}) {
            my @speciesids = keys %{$taxon->{'speciesid'}};
            foreach my $speciesid (@speciesids) {
                next if (!exists($taxon->{'speciesid'}->{$speciesid}->{$taxid}));
                # Could use 0 to signal the strains not covered
#                next if ($taxon->{'speciesid'}->{$speciesid}->{$taxid}!=1);
                my $refacc;
                $refacc = $taxon->{'speciesid'}->{$speciesid}->{$speciesid};
                $debug && print "$subname: \$taxid=$taxid \$speciesid=$speciesid \$refacc=$refacc\n";
                $refseq_acc = $refacc if ($refacc && $refacc =~ /^NC_0/i);
                $debug && print "$subname: \$refacc=$refacc \$refseq_acc=$refseq_acc\n";
                last;
            }

        }

    return $refseq_acc;
} # sub get_refseq_acc


=head2 get_gene_symbol

Takes a mat_peptide feature from refseq object, returns the gene symbol. The gene symbols are
 defined as hash (accession) of hash (location of matpeptide) here, since the values in genbank
 file are inconsistent or unavailable for every mat_peptide

Returns the gene symbol in a string

=cut

sub get_gene_symbol {
    my ($reffeat) = @_;

    my $debug = 0 && $debug_all;
#    print STDERR "get_gene_symbol: \$reffeat=\n".Dumper($reffeat)."end of \$reffeat\n\n";

    # These gene symbols are defined by CLarsen, after considering refseqs, e.g. NC_001477 & NC_009942
    my $gene_symbols;
    $gene_symbols = {
        # Dengue virus 1
        'NC_001477' => {
               '95..436' => 'ancC',
               '95..394' => 'C',
              '437..934' => 'preM',
              '710..934' => 'M',
             '935..2419' => 'E',
            '2420..3475' => 'NS1',
            '3476..4129' => 'NS2a',
            '4130..4519' => 'NS2b',
            '4520..6376' => 'NS3',
            '6377..6757' => 'NS4a',
            '6758..6826' => '2k',
            '6827..7573' => 'NS4b',
           '7574..10270' => 'NS5',
        },
        # Dengue virus 2
        'NC_001474' => {
               '97..438' => 'ancC',
               '97..396' => 'C',
              '439..936' => 'preM',
              '712..936' => 'M',
             '937..2421' => 'E',
            '2422..3477' => 'NS1',
            '3478..4131' => 'NS2a',
            '4132..4521' => 'NS2b',
            '4522..6375' => 'NS3',
            '6376..6756' => 'NS4a',
            '6757..6825' => '2k',
            '6826..7569' => 'NS4b',
           '7570..10269' => 'NS5',
        },
        # Dengue virus 3
        'NC_001475' => {
               '95..436' => 'ancC',
               '95..394' => 'C',
              '437..934' => 'preM',
              '710..934' => 'M',
             '935..2413' => 'E',
            '2414..3469' => 'NS1',
            '3470..4123' => 'NS2a',
            '4124..4513' => 'NS2b',
            '4514..6370' => 'NS3',
            '6371..6751' => 'NS4a',
            '6752..6820' => '2k',
            '6821..7564' => 'NS4b',
           '7565..10264' => 'NS5',
        },
        # Dengue virus 4
        'NC_002640' => {
              '102..440' => 'ancC',
              '102..398' => 'C',
              '441..938' => 'preM',
              '714..938' => 'M',
             '939..2423' => 'E',
            '2424..3479' => 'NS1',
            '3480..4133' => 'NS2a',
            '4134..4523' => 'NS2b',
            '4524..6377' => 'NS3',
            '6378..6758' => 'NS4a',
            '6759..6827' => '2k',
            '6828..7562' => 'NS4b',
           '7563..10262' => 'NS5',
        },

        # Japanese encephalitis virus
        'NC_001437' => {
              '96..>476' => 'ancC',
              '96..>410' => 'C',
            '<477..>977' => 'preM',
            '<753..>977' => 'M',
           '<978..>2477' => 'E',
          '<2478..>3533' => 'NS1',
          '<3534..>4214' => 'NS2a',
          '<4215..>4607' => 'NS2b',
          '<4608..>6464' => 'NS3',
          '<6465..>6842' => 'NS4a',
          '<6843..>6911' => '2k',
          '<6912..>7676' => 'NS4b',
           '7677..10391' => 'NS5',
        },

        # West Nile virus
        'NC_001563' => {
              '97..>465' => 'ancC',
              '97..>411' => 'C',
            '<466..>966' => 'preM',
            '<742..>966' => 'M',
           '<967..>2457' => 'E',
          '<2458..>3513' => 'NS1',
          '<3514..>4206' => 'NS2a',
          '<4207..>4599' => 'NS2b',
          '<4600..>6456' => 'NS3',
          '<6457..>6834' => 'NS4a',
          '<6835..>6903' => '2k',
          '<6904..>7671' => 'NS4b',
           '7672..10386' => 'NS5',
        },

        # Tick-borne encephalitis virus
        'NC_001672' => {
             '133..>468' => 'ancC',
             '133..>420' => 'C',
            '<469..>972' => 'preM',
            '<748..>972' => 'M',
           '<973..>2460' => 'E',
          '<2461..>3516' => 'NS1',
          '<3517..>4206' => 'NS2a',
          '<4207..>4599' => 'NS2b',
          '<4600..>6462' => 'NS3',
          '<6463..>6840' => 'NS4a',
          '<6841..>6909' => '2k',
          '<6910..>7665' => 'NS4b',
           '7666..10374' => 'NS5',
        },

        # Yellow fever virus (YFV)
        'NC_002031' => {
             '119..>481' => 'ancC',
             '119..>421' => 'C',
            '<482..>973' => 'preM',
            '<749..>973' => 'M',
           '<974..>2452' => 'E',
          '<2453..>3508' => 'NS1',
          '<3509..>4180' => 'NS2a',
          '<4181..>4570' => 'NS2b',
          '<4571..>6439' => 'NS3',
          '<6440..>6817' => 'NS4a',
          '<6818..>6886' => '2k',
          '<6887..>7636' => 'NS4b',
           '7637..10351' => 'NS5',
        },

        # St. Louis encephalitis virus
        'NC_007580' => {
              '99..>461' => 'ancC',
              '99..>407' => 'C',
            '<462..>962' => 'preM',
            '<738..>962' => 'M',
           '<963..>2465' => 'E',
          '<2466..>3521' => 'NS1',
          '<3522..>4202' => 'NS2a',
          '<4203..>4595' => 'NS2b',
          '<4596..>6449' => 'NS3',
          '<6450..>6827' => 'NS4a',
          '<6828..>6896' => '2k',
          '<6897..>7670' => 'NS4b',
           '7671..10388' => 'NS5',
        },

        # Hepatitis C virus genotype 1, used for all genotypes of Hepatitis C
        'NC_004102' => {
             '342..914' => 'C',
            '915..1490' => 'E1',
           '1491..2579' => 'E2',
           '2580..2768' => 'p7',
           '2769..3419' => 'NS2',
           '3420..5312' => 'NS3',
           '5313..5474' => 'NS4a',
           '5475..6257' => 'NS4b',
           '6258..7601' => 'NS5a',
           '7602..9374' => 'NS5b',
        },

        # Norwalk virus
        'NC_001959' => {
              '5..1198' => 'Nterm',
           '1199..2287' => 'NTPase',
           '2288..2890' => 'p22',
           '2891..3304' => 'VPg',
           '3305..3847' => 'Pro',
           '3848..5371' => 'Pol',
        },

    };

    my $gene_symbol = 'unk';
    my $acc = $reffeat->seq->accession_number;
    my $loc = $reffeat->location->to_FTstring;
    $gene_symbol = $gene_symbols->{$acc}->{$loc} ? $gene_symbols->{$acc}->{$loc} : $gene_symbol;
    $debug && print "get_gene_symbol: primary_tag=".($reffeat->seq->accession_number);
    $debug && print "\tlocation=".($reffeat->location->to_FTstring);
    $debug && print "\tsymbol=$gene_symbol\n";

    return $gene_symbol;
} # sub get_gene_symbol


## sub get_dna_byloc fetches the sequence directly from seq string according to
## the location range. Sequence is always extracted from $seq_obj_seq
sub get_dna_byloc {
    my ($feat_obj, $seq_obj_seq) = @_;

      my $debug = 0 && $debug_all;
      $debug && print STDERR "get_dna_byloc: \$feat_obj = \n".Dumper($feat_obj)."End of \$feat_obj\n\n";
      my ($s, $product);
      if ($feat_obj->location->isa('Bio::Location::Simple')) {
#           print "  Simple:DNA    :seq=". $feat_obj->seq->seq ."\n";
#           $s = substr($seq_obj_seq, $feat_obj->location->start-1, $feat_obj->location->end - $feat_obj->location->start + 1);
           my $start = $feat_obj->location->start-1;
           my $len = $feat_obj->location->end - $feat_obj->location->start + 1;
           $debug && print STDERR "get_dna_byloc: \$start = $start \$len=$len\n";
           if ($feat_obj->has_tag('codon_start')) {
               my $codon_start = [$feat_obj->get_tag_values('codon_start')];
               $codon_start = $codon_start->[0];
               $start = $start + $codon_start -1;
               $len = $len - $codon_start +1;
               $debug && print STDERR "get_dna_byloc: \$codon_start = $codon_start\n";
           }
           $debug && print STDERR "get_dna_byloc: \$start = $start \$len=$len\n";
           $s = substr($seq_obj_seq, $start, $len);
#           print "  Simple:protein:seq=". $feat_obj->seq->translate()->seq ."\n";
      } elsif ($feat_obj->location->isa('Bio::Location::Split')) {
           my $s0 = $seq_obj_seq;
#           print "location=". $feat_obj->location->to_FTstring . "\n";
           for my $loc ($feat_obj->location->sub_Location) {
#              $s .= substr($s0, $loc->start -1, $loc->end +1 - $loc->start);
               my $start = $loc->start -1;
               my $len = $loc->end +1 - $loc->start;
              if ($feat_obj->has_tag('codon_start')) {
                  my $codon_start = [$feat_obj->get_tag_values('codon_start')];
                  $codon_start = $codon_start->[0];
                  $start = $start + $codon_start -1;
                  $len = $len - $codon_start +1;
                  $debug && print STDERR "get_dna_byloc: \$codon_start = $codon_start\n";
              }
              $s .= substr($s0, $start, $len);
#              print "start=". $loc->start ."\tend=". $loc->end ."\n";
#              print "  \$s=". $s ."\n";
           }
      } elsif ($feat_obj->location->isa('Bio::Location::Fuzzy')) {
           my $s0 = $seq_obj_seq;
           $debug && print STDERR "get_dna_byloc: location=". $feat_obj->location->to_FTstring . "\n";
           for my $loc ($feat_obj->location->each_Location) {
#              $s .= substr($s0, $loc->start -1, $loc->end +1 - $loc->start);
               my $start = $loc->start -1;
               my $len = $loc->end +1 - $loc->start;
              if ($feat_obj->has_tag('codon_start')) {
                  my $codon_start = [$feat_obj->get_tag_values('codon_start')];
                  $codon_start = $codon_start->[0];
                  $start = $start + $codon_start -1;
                  $len = $len - $codon_start +1;
                  $debug && print STDERR "get_dna_byloc: \$codon_start = $codon_start\n";
              }
              $s .= substr($s0, $start, $len);
#              $debug && print STDERR "get_dna_byloc: start=". $loc->start ."\tend=". $loc->end ."\n";
#              $debug && print STDERR "get_dna_byloc:   \$s=". $s ."\n";
           }
      }

    return $s;
} # sub get_dna_byloc


=head2 check_refseq

Takes 1 refseq and 1 target genome. Check is the taxid for both are same

=cut

sub check_refseq {
    my ($refseq, $inseq) = @_;

    my $debug = 0 && $debug_all;
    my $refseq_fit = 0;
    return $refseq_fit if (!$refseq || !$inseq);

    my $taxid = $inseq->species->ncbi_taxid;
    my $ref_taxid = $refseq->species->ncbi_taxid;
    if ($taxid == $ref_taxid) {
      $refseq_fit = 1;
    }

    return $refseq_fit;
} # sub check_refseq


=head2 convert_range_protein2dna

Takes a CDS, and an array of gaps from alignment,
 return the gaps converted to DNA range, taking into account the start of CDS

=cut

sub convert_range_protein2dna {
    my ($f, $gaps) = @_;

    my $debug = 1 && $debug_all;
    my $subname = 'convert_range_protein2dna';
    my @gaps2;
    my $f_start = $f->location->start;
    my $under_water = 0;
    my $above_water = 0;
    $debug && print STDERR "$subname: gaps are converted to DNA coordinates from protein coordinate.\n";
    $debug && print STDERR "$subname: DNA coordinates include the start of CDS, also include the gaps.\n";
    for (my $i=0; $i<=$#{@$gaps}; $i++) {
        my $gap = $gaps->[$i];

        if ($gap eq '') {
                push @gaps2, $gap;
                next;
        }

        # This indicates the gap at the beginning of the CDS, special case
        if ($gap == $i) {
            my $sss = $f_start - ($gap*3)+0;
            $debug && print STDERR "$subname: \$i=$i \$gap=$gap start=$f_start new=". ($f_start-($gap*3)+2) ." \$under_water=$under_water \$sss=$sss\n";
            unshift @gaps2, (shift @gaps2, $f_start - ($gap*3)+2);
            if ($gaps2[1] <= 0) {
                $under_water++;
            } else {
                $above_water++;
            }
            unshift @gaps2, (shift @gaps2, $f_start - ($gap*3)+1);
            if ($gaps2[1] <= 0) {
                $under_water++;
            } else {
                $above_water++;
            }
            unshift @gaps2, (shift @gaps2, $f_start - ($gap*3)+0);
            if ($gaps2[1] <= 0) {
                $under_water++;
            } else {
                $above_water++;
            }
        } else {
            my $sss = $f_start -$under_water -$above_water + ($gap-1)*3;
            $debug && print STDERR "$subname: \$i=$i \$gap=$gap start=$f_start new=". ($f_start + ($gap-1)*3) ." \$under_water=$under_water \$sss=$sss\n";
            push @gaps2, $f_start -$under_water -$above_water + ($gap-1)*3;
            push @gaps2, $f_start -$under_water -$above_water + ($gap-1)*3+1;
            push @gaps2, $f_start -$under_water -$above_water + ($gap-1)*3+2;
        }
    }

    # check the gaps to ensure all gaps are included
    my $f_end = $f->location->end;
    for (my $i=1; $i<=$#gaps2; $i++) {
        $debug && print STDERR "$subname: \$i=$i \$gaps[$i]=$gaps2[$i] \$f_end=$f_end\n";
        if ($gaps2[$i] < $f_end) {
            $f_end++;
        } elsif ($gaps2[$i] == $f_end) {
            last;
        } elsif ($gaps2[$i] > $f_end) {
            @gaps2 = (@gaps2, $f_end+1 .. $gaps2[$i]-1);
            last;
        }
    }

    $debug && print STDERR "$subname: \$under_water=$under_water \$above_water=$above_water\n";
    my $gap1 = shift @gaps2;
    $debug && print STDERR "$subname: sorting \@gaps2=@gaps2\n";
    @gaps2 = sort {$a<=>$b} @gaps2;
    @gaps2 = ($gap1, @gaps2);
    $debug && print STDERR "$subname: sorted \@gaps2=@gaps2\n";
    return (\@gaps2);
} # sub convert_range_protein2dna


=head2 get_aln_seq

Takes an alignment object and accession, search for the alignment seq that contains the accession
 return id of the seq in alignment

=cut

sub get_aln_seq {
    my ($aln, $acc) = @_;

    my $debug = 0 && $debug_all;
    my $id;
    my $seqs = [ $aln->each_seq() ];
    foreach my $seq (@$seqs) {
        if ($seq->id =~ /$acc/i) {
            $id = $seq->id;
            last;
        }
    }
#    $refid = 'gi|22129792|ref|NC_004102';
    $debug && print "annotate_msa: \$id=$id.\n";
    my $refaln = [ $aln->each_seq_with_id($id) ];

    return ($id);
} # sub get_aln_seq


=head2 get_new_translation

Takes a new feature object, and a CDS
 return translation

=cut

sub get_new_translation {
    my ($feat, $cds) = @_;

    my $debug = 0 && $debug_all;
    my $subname = 'get_new_translation';
#    my $translation = '';
    my $translation = undef; # emptry string '' suppresses error msg, but cause other problems
    my $acc = $cds->seq->accession_number;

    $debug && print STDERR "$subname: \$feat = \n".Dumper($feat)."End of \$feat\n\n";
    {
        my $overhang1 = ($feat->location->start - $cds->location->start) % 3;
        my $overhang2 = ($feat->location->end   - $cds->location->start + 1) % 3;
        if ($overhang1 || $overhang2 ) {
            $debug && print STDERR "$subname: ERROR: \$acc=$acc \$feat = '".$feat->location->to_FTstring."' \$cds = '".$cds->location->to_FTstring."'\n";
            $debug && print STDERR "$subname: ERROR: \$acc=$acc \$overhang1=$overhang1 \$overhang2=$overhang2\n";
#            croak "$subname: ERROR: \$overhang1=$overhang1 \$overhang2=$overhang2";
        }
    }

    my $s = Annotate_Util::get_dna_byloc( $feat, $cds->{_gsf_seq}->seq);
#    $debug && print STDERR "$subname: \$s  ='$s'\n";
    my $f = Bio::PrimarySeq->new(
                         -seq      => $s,
                         -id       => '',	# id can't contain space
                         -alphabet => 'dna'
                                );
    $f->revcom() if ($feat->strand() == -1);

    $s = $f->translate()->seq;
    $debug && print STDERR "$subname: \$s  ='$s'\n";
    $debug && print STDERR "$subname: \$s  =".length($s)."\n";

    # Only keep the annotated mat_peptide that confirms to CDS. It should, just to double check
    my @parent_cds_seq = $cds->get_tag_values('translation');
    $debug && print STDERR "$subname: \$cds  =$parent_cds_seq[0]\n";
    $debug && print STDERR "$subname: \$cds  =".length($parent_cds_seq[0])."\n";

    if ($parent_cds_seq[0] =~ /$s/i) {
        $translation = $s;
    } else {
        print STDERR "$subname: ERROR: \$acc=$acc translation for feature doesn't match CDS.\n";
        print STDERR "$subname: \$s  ='$s'\n";
        print STDERR "$subname: \$cds='".$parent_cds_seq[0]."'\n";
        print STDERR "$subname: diff=".Annotate_Verify::diff_2str( $s, $parent_cds_seq[0])."\n";
#        return undef;
    }

    $debug && print STDERR "$subname: \$acc=$acc \$translation  ='$translation'\n";
    return ($translation);
} # sub get_new_translation


=head2 assemble_new_feature

Takes refcds, reffeat, loc2 of new feature, and target cds.
 Return a new Bio::SeqFeature::Generic that has:
 primary_tag from reffeat
 location from loc2
 translation from target cds based on loc2
 product from reffeat
 note="Annotated by VIPRBRC, MUSCLE, refseq=$reffeat->accession_number"
 note="Desc:CDS=GI:11559447|Loc=311..883|gene_symbol=unk|mat_peptide=core protein"

=cut

sub assemble_new_feature {
    my ($refcds, $reffeat, $cds, $aln, $aln_q, $aln_h, $note) = @_;

    my $debug = 0 && $debug_all;
    my $subname = 'assemble_new_feature';

    my ($loc2, $errcode) = Annotate_Util::msa_get_feature_loc(
                                 $reffeat->location,
                                 $refcds->location,
                                 $cds,
                                 $aln,
                                 $aln_q,
                                 $aln_h,
                                 );
    $debug && print STDERR "$subname: \$loc2 = \n".Dumper($loc2)."End of \$loc2\n\n";
    if (!$loc2) {
        $debug && print STDERR "$subname: \$loc2 is empty\n";
        return undef;
    }
    if ($errcode->{OutsideCDS}==1) {
        print STDERR "$subname: \$loc2=".$loc2->to_FTstring." doesn't match \$cds=".$cds->location->to_FTstring."\n";
        return undef;
    }
    if ($errcode->{long_internal_gap}==1) {
        print STDERR "$subname: \$loc2=".$loc2->to_FTstring." has too long internal gap in alignment, discard.\n";
        return undef;
    }

    my $feat;
    $feat = Bio::SeqFeature::Generic->new();
    $feat->primary_tag($reffeat->primary_tag); # Take the primary_tag from reffeat
    $debug && print STDERR "$subname: \$feat is ".ref($feat)." \n";

    $feat->location($loc2);
    $debug && print STDERR "$subname: \$loc2=\n".Dumper($loc2)."end of \$loc2\n\n";

#    $feat->attach_seq($cds->seq);

    $feat->add_tag_value('note', $note);
    $feat->strand($reffeat->strand());

    if ($reffeat->has_tag('product')) {
        my $tag = [$reffeat->get_tag_values('product')];
        $feat->add_tag_value('product', $tag->[0]);
    }

    if ($cds->has_tag('codon_start')) {
        my $codon_start = [$cds->get_tag_values('codon_start')];
        $codon_start = $codon_start->[0];
#        $feat->add_tag_value('codon_start', $codon_start);
    }

    # now get translation
    my $translation = '';
    if ($cds) {
        my $s = Annotate_Util::get_new_translation( $feat, $cds);
        if ($s) {
            $translation = $s;
        } else {
            print STDERR "$subname: \$feat = \n".Dumper($feat)."End of \$feat\n\n";
            print STDERR "$subname: WARNING: translation for mat_peptide doesn't match CDS.\n";
            print STDERR "$subname: \$s='$s'\n";
            print STDERR "$subname: \$cds='".$cds->seq->translate->seq."'\n";
            return undef;
        }

    }
    $feat->add_tag_value('translation', $translation);
#    print STDERR "$subname: \$feat=\n".Dumper($feat)."end of \$feat\n\n";

    if ($cds->has_tag('locus_tag')) {
        my $tag = [$cds->get_tag_values('locus_tag')];
        $feat->add_tag_value('locus_tag', $tag->[0]);
    }

    my $gene_symbol = Annotate_Util::get_gene_symbol( $reffeat);
    my ($id, $desc) = &get_feature_id_desc($feat, $errcode, $gene_symbol, $cds);
    $feat->add_tag_value('note', 'Desc:'.$desc);

    $debug && print STDERR "$subname: \$loc = \n".Dumper($reffeat->location)."End of \$loc\n\n";
    $debug && print STDERR "$subname: \$loc2 = \n".Dumper($loc2)."End of \$loc2\n\n";

    return ($feat);
} # sub assemble_new_feature


=head2 get_feature_id_desc

Takes $feat, $errcode, $refcds, $reffeat, $cds, $note
 Return a suitable id and desc for the new feature
 note="Annotated by VIPRBRC, MUSCLE, refseq=$reffeat->accession_number"
 note="Desc:CDS=GI:11559447|Loc=311..883|gene_symbol=unk|mat_peptide=core protein"

=cut

sub get_feature_id_desc {
    my ($feat, $errcode, $gene_symbol, $cds) = @_;

    my $debug = 0 && $debug_all;
    my $subname = 'get_feature_id_desc';

    my (@id);
    my $id = 0;
    if ($cds) {
        # get the id from parent CDS feature, 1) GI, 2) NP
        if ($cds->has_tag('db_xref')) {
            @id = $cds->get_tag_values('db_xref');
            for my $id1 (@id) {
               if ($id1 =~ /^GI:/i) {
                    $id = 'CDS='. $id1 ."|";
                    last;
               }
            }
        }
        if (($id !~ /^CDS=GI:/) && $cds->has_tag('protein_id')) {
            @id = $cds->get_tag_values('protein_id');
            $id = 'CDS='. $id[0] ."|";
        }

        # In case $id is still empty
        if (($id !~ /^CDS=/) ) {
            my @tags = ('db_xref', 'protein_id'); # per Chris request
            for my $tag (@tags) {
                if ($cds->has_tag($tag)) {
                   @id = $cds->get_tag_values($tag);
                   $id .= 'CDS='. $id[0] ."|";
                   last;
                }
            }
        }
    }
    $id = "CDS=unknown|" if (($id !~ /^CDS=/) ); # indicates a potential problem

    # add the accession
    if ( 1 ) {
        $id = 'ACC='. $cds->seq->accession_number .'|'. $id;
    }

    # add the range of mat_peptide
    my $allow_fuzzy_location = 0;
    if ($allow_fuzzy_location) {
        $id .= 'Loc='. $feat->location->to_FTstring .'|';
    } else {
        $id .= 'Loc='. $feat->location->start .'..'. $feat->location->end .'|';
    }

    # add gene symbol
    my $desc = $id;
    $desc .= 'gene_symbol='. $gene_symbol .'|';

    # error code: 0: none; 1: gaps around cleavage site; 2: partial mat_peptide; 3: outside of CDS
    if (exists($errcode->{GapAtCleavage}) && $errcode->{GapAtCleavage}==1) {
        $desc = $desc. 'GapAtCleavage=Y|'; # This could be removed in fix_cleavage_gaps
    }
    if (exists($errcode->{OutsideCDS}) && $errcode->{OutsideCDS}==1) {
       $desc = $desc. 'OutsideCDS=Y|';
    }
    if (exists($errcode->{partial_mat_peptide}) && $errcode->{partial_mat_peptide}==1) {
       $desc = $desc. 'Partial=Y|';
    } else {
#       $desc = $desc. 'Partial=N|';
    }

    # add the product of mat_peptide
    my @tags = ('product');
    for my $tag (@tags) {
        if ($feat->has_tag($tag)) {
           @id = $feat->get_tag_values($tag);
           if ( 1 ) {
               $desc = $desc. 'product='.$id[0];
#               $desc = $desc. 'mat_peptide='.$id[0];
           } else {
               $desc = $desc. $feat->primary_tag.'='. $id[0];
           }
           $id .= 'product='.$id[0];
           last;
        }
    }

    return ($id, $desc);
} # sub get_feature_id_desc


=head2 fix_cleavage_gaps

Takes the new annotation in an array, refcds, cds, and the gaps.
 Check if all new mat_peptides are connected
 if not, decide if the missing sequence at the cleavage site should go to the tail of previous mat_peptide,
 or head of the next mat_peptide

If the alignment is not definite, as in EF407458 EF407463 EF407467 when run together and separately, the annotation will change according to the MSA. This is a consequence of the shifting alignment.

=cut

sub fix_cleavage_gaps {
    my ($feats_all, $refset, $gaps_q, $gaps_h, $note) = @_;

    my $debug = 0 && $debug_all;
    my $subname = 'fix_cleavage_gaps';

    my $cds = $feats_all->[0];
    my $refcds = $refset->[0];
    my $max_length = 10;
    my $has_gap = 0;
    my $feats;
    $debug && print STDERR "$subname: \$feats_all=\n".Dumper($feats_all)."end of \$feats_all\n";
    $debug && print STDERR "$subname: \$refset=\n".Dumper($refset)."end of \$refset\n";

    my $feat1_n = 1;
    for (my $i=2; $i<=$#{@$feats_all}; $i++) {
        my ($feat1, $feat2);
        my ($str, $l, $n, $f);
        if ($refset->[$i]->location->end+1 == $refset->[$i+1]->location->start) {
#            $feat1_n = $i-1;
        }
        $debug && print STDERR "$subname: \$i=$i \$feat1_n=$feat1_n\n";
        $feat1 = $feats_all->[$feat1_n];
        $feat2 = $feats_all->[$i];

        # Skip if $feat2 starts before the end of $feat1, mostly indicating $feat2 is a part of $feat1
        next if ($feat1->location->end > $feat2->location->start);
        # Skip if there is no gap at cleavage site
        if ($feat1->location->end+1 == $feat2->location->start) {
            $feat1_n = $i;
            next;
        }

        $has_gap = 1;
        print STDERR "$subname: \$i=$i \$feat1_n=$feat1_n\n";
        $debug && print STDERR "$subname: \$feat1=\n".Dumper($feat1)."end of \$feat1\n";
        $debug && print STDERR "$subname: \$feat2=\n".Dumper($feat2)."end of \$feat2\n";

        if ($feat1->location->end+1 != $feat2->location->start) {
            # Update the dscription to include GapAtCleavage
            my @values = $feat2->get_tag_values('note');
            for (my $k=0; $k<=$#values; $k++) {
                my $value = $values[$k];
                next if ($value !~ /^Desc:(.+)$/i);
                my $pattn = 'gene_symbol=[a-zA-Z0-9]+';
                if ($value =~ s/($pattn)/$1\|GapAtCleavage=Y/i) {
                    $values[$k] = $value;
#                    $debug && print STDERR "$subname: \$value123 = '$value'\n";
                    $debug && print STDERR "$subname: \@values = \n".Dumper(@values)."End of \@values\n\n";
                }
            }
            $feat2->remove_tag('note');
            $feat2->add_tag_value('note', @values);
        }

        print STDERR "$subname: Found a gap after #$i: ".$feat1->location->end ."..". $feat2->location->start."\n";

        my $ts = ['', '', ''];
        $str = [ $feat1->get_tag_values('translation') ]->[0];
        $l = length($str);
        $n = ($l>=$max_length) ? $max_length : $l;
        $ts->[0] = substr($str, $l-$n, $n);

        $str = [ $feat2->get_tag_values('translation') ]->[0];
        $l = length($str);
        $n = ($l>=$max_length) ? $max_length : $l;
        $ts->[2] = substr($str, 0, $n);

        my $loc2 = Bio::Location::Simple->new();
        $loc2->start($feat1->location->end+1);
        $loc2->end  ($feat2->location->start-1);
        $debug && print STDERR "$subname: \$loc2=\n".Dumper($loc2)."end of \$loc2\n";

        my $feat;
        $feat = Bio::SeqFeature::Generic->new();
        $feat->primary_tag('miscel');
        $feat->location($loc2);
        $feat->strand($cds->strand());

        # now get translation
        my $translation = Annotate_Util::get_new_translation( $feat, $cds);
        $feat->add_tag_value('translation', $translation);
        $ts->[1] = $translation;

        print STDERR "$subname: \$ts = '@$ts'\n";


        my $rs = ['', '', ''];
        my $reffeat1 = $refset->[$feat1_n+1];
        my $reffeat2 = $refset->[$i+1];
        $debug && print STDERR "$subname: \$reffeat1=\n".Dumper($reffeat1)."end of \$reffeat1\n";
        $debug && print STDERR "$subname: \$reffeat2=\n".Dumper($reffeat2)."end of \$reffeat2\n";
        $str = $reffeat1->seq->seq;
        $f = Bio::PrimarySeq->new(-seq => $str, -id => '', -alphabet => 'dna');
        $f->revcom() if ($reffeat1->strand == -1);
#        $debug && print STDERR "$subname: \$f=\n".Dumper($f)."end of \$f\n";
        $str = $f->translate()->seq;
        $l = length($str);
        $n = ($l>=$max_length) ? $max_length : $l;
        $rs->[0] = substr($str, $l-$n, $n);

        $str = $reffeat2->seq->seq;
        $f = Bio::PrimarySeq->new(-seq => $str, -id => '', -alphabet => 'dna');
        $f->revcom() if ($reffeat1->strand == -1);
        $str = $f->translate()->seq;
        $l = length($str);
        $n = ($l>=$max_length) ? $max_length : $l;
        $rs->[2] = substr($str, 0, $n);

        print STDERR "$subname: \$ts = '@$ts'\n";
        print STDERR "$subname: \$rs = '@$rs'\n";

        # The gap would be lumped into the beginning of next mat_peptide provided that
        # 1) there are min. of 5 identical residues in the last 10 residues of the previous mat_peptide;
        # 2) the last residue in refseq and target are identical or from same class: (AG) (ST) (RK) (C);
        # 3) the gap is 5 or less residues long.
        my $identical_last10 = 0;
        my $id_last_residue = 0;
        my $length_gap = 0;

        $n = (length($rs->[0])<=length($ts->[0])) ? length($rs->[0]) : length($ts->[0]);
        foreach my $j (0 .. $n-1) {
            my $rc = uc substr($rs->[0], length($rs->[0])-$j, 1);
            my $tc = uc substr($ts->[0], length($ts->[0])-$j, 1);
            $identical_last10++ if ($rc eq $tc);
        }

        my $rc = uc substr($rs->[0], length($rs->[0]), 1);
        my $tc = uc substr($ts->[0], length($ts->[0]), 1);
        if ($rc eq $tc) {
            $id_last_residue = 1;
        } else {
            my $groups = {
                'A' => 'AG',
                'R' => 'KR',
                'N' => 'NQ',
                'D' => 'DE',
                'C' => 'C',
                'Q' => 'NQ',
                'E' => 'DE',
                'G' => 'AG',
                'H' => 'H',
                'I' => 'ILV',
                'L' => 'ILV',
                'K' => 'KR',
                'M' => 'M',
                'F' => 'FY',
                'P' => 'P',
                'S' => 'ST',
                'T' => 'ST',
                'W' => 'W',
                'Y' => 'FY',
                'V' => 'ILV',
                };
            $id_last_residue = 1 if ($groups->{$rc} eq $groups->{$tc});
        }

        $length_gap = length($ts->[1]);

        my $msg;
        if ($identical_last10>6 && $id_last_residue && $length_gap<=5) {
            my $loc2 = $feat2->location;
            if ($loc2->isa('Bio::Location::Simple') ) {

                my $start = $loc2->start - $length_gap*3;
                $loc2->start($start);

                $debug && print STDERR "$subname: \$loc2=\n".Dumper($loc2)."end of \$loc2\n";

            } elsif ($loc2->isa('Bio::Location::Fuzzy')) {

                my $min_start;
                $min_start = $loc2->min_start - $length_gap*3;
                $loc2->min_start($min_start);
                my $max_start;
                $max_start = $loc2->max_start - $length_gap*3;
                $loc2->max_start($max_start);

                $debug && print STDERR "$subname: \$loc2=\n".Dumper($loc2)."end of \$loc2\n";

            }
            my $translation = '';
            if ($cds) {
                my $s = Annotate_Util::get_new_translation( $feat2, $cds);
                if ($s) {
                    $translation = $s;
                } else {
                    print STDERR "$subname: \$feat = \n".Dumper($feat)."End of \$feat\n\n";
                    print STDERR "$subname: ERROR: translation for mat_peptide doesn't match CDS.\n";
                    print STDERR "$subname: \$s='$s'\n";
                    print STDERR "$subname: \$cds='".$cds->seq->translate->seq."'\n";
                    return undef;
                }

            }
            $feat2->remove_tag('translation');
            $feat2->add_tag_value('translation', $translation);

            # Update the dscription
            my @values = $feat2->get_tag_values('note');
            for (my $k=0; $k<=$#values; $k++) {
                my $value = $values[$k];
                next if ($value !~ /^Desc:(.+)$/i);
                if ($value =~ /(Loc=[0-9]+[.]{2,2}[0-9]+)/i) {
                    my $loc_old = $1;
                    my $loc_new = 'Loc='. $feat2->location->start .'..'. $feat2->location->end .'';
                    $debug && print STDERR "$subname: \$loc_old = '$loc_old'\n";
                    $debug && print STDERR "$subname: \$loc_new = '$loc_new'\n";
                    $value =~ s/$loc_old/$loc_new/;
                    $values[$k] = $value;
                    $debug && print STDERR "$subname: \$value = '$value'\n";
                    $debug && print STDERR "$subname: \@values = \n".Dumper(@values)."End of \@values\n\n";
                }
            }
            $feat2->remove_tag('note');
            $feat2->add_tag_value('note', @values);

            $msg =  "$subname: WARNING: There is a gap ($ts->[1]) between mat_peptides #$feat1_n and #$i:'$ts->[1]'. ";
            $msg .= " The gap has been lumped into the next mat_peptide\n";

        } else {
          $msg =  "$subname: ERROR: There is a gap between mat_peptides #".($feat1_n+1)." and #";
          $msg .=  ($i+1) .":$ts->[1]\n";
          if ($identical_last10<=6) {
            $msg .= "$subname: ERROR: Only $identical_last10 identical AAs amoung the last 10 of mat_peptide #$i\n";
            $msg .= "$subname: ERROR: This is less than required 7 or more\n";
          }

          if (!$id_last_residue) {
            $msg .= "$subname: ERROR: Last residue of mat_peptide #$i is $tc, differs from that in refseq $rc\n";
          }

          if ($length_gap>5) {
            $msg .= "$subname: ERROR: Length of gap between #$feat1_n-$i is $length_gap, more than required 5.\n";
          }

        }
        print STDERR $msg;
        print STDOUT $msg;

        if ($feat1->location->end+1 == $feat2->location->start) {
            # Update the dscription
            my @values = $feat2->get_tag_values('note');
            for (my $k=0; $k<=$#values; $k++) {
                my $value = $values[$k];
                next if ($value !~ /^Desc:(.+)$/i);
                my $pattn = 'GapAtCleavage=Y\|';
                if ($value =~ s/($pattn)//i) {
                    $values[$k] = $value;
                    $debug && print STDERR "$subname: \$value = '$value'\n";
                    $debug && print STDERR "$subname: \@values = \n".Dumper(@values)."End of \@values\n\n";
                }
            }
            $feat2->remove_tag('note');
            $feat2->add_tag_value('note', @values);

            $msg =  "$subname: WARNING: There is a gap ($ts->[1]) between mat_peptides #$feat1_n and #$i:$ts->[1]. ";
            $msg .= " The gap has been lumped into the next mat_peptide\n";

        }

        $feat1_n = $i; # Keep the current mat_peptide as starting point.
        $debug && print STDERR "$subname: \$feat1_n = $feat1_n\n";
    }

    $has_gap && print STDERR "$subname: \$feats_all = \n".Dumper($feats_all)."End of \$feats_all\n\n";

    return ($feats);
} # sub fix_cleavage_gaps


=head2 do_bl2seq_search

Takes 2 sequence strings to align against each other.
Runs a pairwise bl2seq query looking for a segment of the second sequence
which matches the first string.

Returns a object of Bio::Search::Result::BlastResult

=cut

sub run_bl2seq_search {
    # Get the shortstring and the query string
    my ($s1string, $s2string) = @_;

    my $debug = 0 && $debug_all;
    return unless ($s1string && $s2string);
    my $td = tempdir( CLEANUP => 1 );  # temp dir (threadsafe)
#    print "s1=$s1string\n";
#    print "s2=$s2string\n";

    # Create temp fasta strings
    my $s1 = ">s1\n$s1string\n";
    my $s2 = ">s2\n$s2string\n";
#    my $td = tempdir( CLEANUP => 1 );  # temp dir (threadsafe)
    my $s1temp = "$td/shortstring_input1.fasta";
    my $s2temp = "$td/shortstring_input2.fasta";
    `echo -e \"$s1\" > $s1temp`;
    `echo -e \"$s2\" > $s2temp`;
    my $blast_outfile = "$td/shortstring_blast.tmp";

    # Perform the 'bl2seq' Blast search, which will return the best local
    # alignment between the two sequences:
    #   -p blastp   : Run protein blast
    #   -i [file1]  : Fasta file of first sequence
    #   -i [file2]  : Fasta file of second sequence
    #   -o [output] : Output Blast report filename
    #   -g F        : Turn off gapped alignments
#    `bl2seq -p blastp -i $s1temp -j $s2temp -o $blast_outfile -F F`;
#    `$blast_path/bl2seq -p tblastn -i $s1temp -j $s2temp -o $blast_outfile -F F`;
    `bl2seq -p blastp -i $s1temp -j $s2temp -o $blast_outfile -F F`;

    if ( 1 && $debug ) {
        open my $debugfile, '<', $blast_outfile;
        print STDERR "\n";
#        print STDERR ">>>>>>>>>\$s1=\n$s1string\n";
#        print STDERR ">>>>>>>>>\$s2=\n$s2string\n";
        print STDERR ">>>>>>>>>do_bl2seq_search: $blast_outfile=\n";
        print STDERR <$debugfile>;
        print STDERR ">>>>>>>>>do_bl2seq_search: end of $blast_outfile\n\n";
        close $debugfile;
    }

    my $blast_report  = Bio::SearchIO->new(
                              -format => 'blast',
                              -file => $blast_outfile
                                         );
    my $single_report = $blast_report->next_result; # Bio::Search::Result::BlastResult object
#    print "do_bl2seq_search: \$single_report".Dumper($single_report)."\n";
    $debug && print STDERR "do_bl2seq_search: \$single_report ".ref($single_report)."\n";

    return $single_report;

} # sub run_bl2seq_search


=head2 is_new_annotation

Takes a new annotation,
 and an array of target genomes, in the form of [ accession, [ [CDS, mat_peptides, ...], []... ] ]
 check if the new annotation has same location with any of the existing mat_peptides,
 Return TRUE/FALSE

=cut

sub is_new_annotation {
    my ($feat, $inset) = @_;

    my $debug = 0 && $debug_all;

    $debug && print STDERR "\nis_new_annotation: \$feat = \n".Dumper($feat)."\$feat\n";
    my $new = 1;
    $debug && print STDERR "is_new_annotation: \$new = $new\n";
    for my $gbmatp (@$inset) {
        next if ($gbmatp->primary_tag eq 'CDS');
        my ($str1, $str2);
        if ( 0 ) {
            $str1 = $gbmatp->location->to_FTstring;
            $str2 = $feat->location->to_FTstring;
        } else {
            $str1 = $gbmatp->location->start .'..'. $gbmatp->location->end;
            $str2 = $feat->location->start .'..'. $feat->location->end;
        }
        $debug && print STDERR "is_new_annotation: \$str2='$str2' \$str1='$str1'\n";
        if ($str1 eq $str2) {
            $new = 0;
            last;
        }
    }

    return $new;
} # sub is_new_annotation


1;
