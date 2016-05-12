package Annotate_Math;

use strict;
use warnings;
use English;
use Carp;
use Data::Dumper;

use version; our $VERSION = qv('1.1.2'); # December 01 2010
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
# Annotate_Math contains the core math functions to perform annotation based on MUSCLE
#
#    Authors Chris Larsen, clarsen@vecna.com; Guangyu Sun, gsun@vecna.com
#    December 2010
#
##################

## //EXECUTE// ##


sub adjust_start {
    my ($qstart, $refcds_loc, $cds, $aln, $aln_q, $aln_h) = @_;

    return undef if (!$qstart);
    my $debug = 0 && $debug_all;
    my $subname = 'adjust_start';

        # Get any potential gap caused by split locations in refseq CDS
        my $qstart_gap = Annotate_Math::get_split_gap( $qstart, $refcds_loc);
        $debug && print STDERR "$subname: \$qstart_gap=$qstart_gap\n";
        $qstart -= $qstart_gap if ($qstart_gap % 3 ==0);

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

        # Get any potential gap caused by split locations in CDS
        $qstart_gap = Annotate_Math::get_split_gap( $hstart, $cds->location);
        $debug && print STDERR "$subname: \$qstart_gap=$qstart_gap\n";
        $hstart += $qstart_gap if ($qstart_gap % 3 ==0);

    return $hstart;
} # sub adjust_start

sub adjust_end {
    my ($qend, $refcds_loc, $cds, $aln, $aln_q, $aln_h) = @_;

    return undef if (!$qend);
    my $debug = 0 && $debug_all;
    my $subname = 'adjust_end';

        # Get any potential gap caused by split locations in CDS
        my $qend_gap = Annotate_Math::get_split_gap( $qend, $refcds_loc);
        $debug && print STDERR "$subname: \$qend_gap=$qend_gap\n";
        $qend -= $qend_gap if ($qend_gap % 3 ==0);

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
    $debug && print STDERR "$subname: \$cds_loc->start=".$cds_loc->start." \$qend=$hend\n";

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

        # Get any potential gap caused by split locations in CDS
        $qend_gap = Annotate_Math::get_split_gap( $hend, $cds->location);
        $debug && print STDERR "$subname: \$qend_gap=$qend_gap\n";
        $hend += $qend_gap if ($qend_gap % 3 ==0);

    $debug && print STDERR "$subname: \$qend=$hend\n";
    return $hend;

} # sub adjust_end


=head2 get_split_gap

Takes a position in the DNA sequence, and count how many gaps (U34999:join(44..5695,5699..7540)) prior to it
 Return the number of gaps in DNA unit

=cut

sub get_split_gap {
    my ($pos, $cds_loc) = @_;

    my $debug = 0 && $debug_all;
    my $subname = 'get_split_gap';

        # For split locations, count the gap between sublocations prior to the requested position
        my $split = 0;
        if ($cds_loc->isa('Bio::Location::Split')) {
            my $locs = [ $cds_loc->sub_Location ];
            for my $i (1 .. $#{$locs}) {
               $debug && print STDERR "$subname: \$pos=$pos \$locs->[$i]->start=".$locs->[$i]->start." \$split=$split\n";
               if ($pos < $locs->[$i]->start) {
                 $debug && print STDERR "$subname: 1\n";
                 next;
               } else {
                 $debug && print STDERR "$subname: 2\n";
                 $split = $locs->[$i]->start - $locs->[$i-1]->end -1;
               }
               $debug && print STDERR "$subname: \$pos=$pos \$locs->[$i]->start=".$locs->[$i]->start." \$split=$split\n";
            }
        }

    return $split;
} # sub get_split_gap

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
    if ($reffeat_loc->isa('Bio::Location::Simple') || $reffeat_loc->isa('Bio::Location::Split')) {
        my $qstart = $reffeat_loc->start;
        my $qend = $reffeat_loc->end;
        $debug && print STDERR "$subname: \$qstart=$qstart \$qend=$qend\n";

        $qstart = Annotate_Math::adjust_start( $qstart, $refcds_loc, $cds, $aln, $aln_q, $aln_h);
        $qend = Annotate_Math::adjust_end( $qend, $refcds_loc, $cds, $aln, $aln_q, $aln_h);
        my $codon_start = 1;
        if ($cds->has_tag('codon_start')) {
            $codon_start = [$cds->get_tag_values('codon_start')];
            $codon_start = $codon_start->[0];
        }

        $qstart = $qstart + $codon_start -1;
        $qend = $qend + $codon_start -1;

        # Get any potential gap caused by split locations in CDS
#        my $qstart_gap = Annotate_Math::get_split_gap( $qstart, $cds->location);
#        my $qend_gap = Annotate_Math::get_split_gap( $qend, $cds->location);
#        $debug && print STDERR "$subname: \$qstart_gap=$qstart_gap \$qend_gap=$qend_gap\n";
#        $qstart += $qstart_gap if ($qstart_gap % 3 ==0);
#        $qend += $qend_gap if ($qend_gap % 3 ==0);

        $debug && print STDERR "$subname: \$qstart=$qstart \$qend=$qend\n";

        if ($qstart<=$qend) {
          if ($cds_loc->isa('Bio::Location::Simple') || $cds_loc->isa('Bio::Location::Fuzzy')) {
            $loc2 = Bio::Location::Simple->new();
            $loc2->start($qstart);
            $loc2->end  ($qend);
          } elsif ( $cds_loc->isa('Bio::Location::Split')) {
            my $locs = [ $cds->location->sub_Location ];
            $debug && print STDERR "$subname: \$locs=\n".Dumper($cds->location)."end of \$locs\n\n";
            $debug && print STDERR "$subname: \$locs=\n".Dumper($locs)."end of \$locs\n\n";
            for (my $i=0; $i<=$#{$locs}; $i++ ) {
              my $loc = $locs->[$i];
              $debug && print STDERR "$subname: \$i=$i \$qstart=$qstart, \$loc->start=".$loc->start."\n";
              $debug && print STDERR "$subname: \$i=$i   \$qend=$qend,   \$loc->end=".$loc->end."\n";
              next if ($qstart>$loc->end || $qend<$loc->start); # Skip if the sublocation is outside of the mat_peptide
              if ($loc->start<=$qstart && $qend<=$loc->end) {
                $loc2 = Bio::Location::Simple->new();
                $loc2->start($qstart);
                $loc2->end  ($qend);
              } elsif ($loc->start<=$qstart && $loc->end < $qend) {
                $loc2 = Bio::Location::Split->new();
                my $qs1 = $qstart;
                my $qe1 = $loc->end;
                my $subloc = Bio::Location::Simple->new();
                $subloc->start($qs1);
                $subloc->end  ($qe1);
                $loc2->add_sub_Location($subloc);
              
                $loc = $locs->[++$i];
                $debug && print STDERR "$subname: \$loc=\n".Dumper($loc)."end of \$loc\n\n";
                while ($i<=$#{$locs} && $loc->start < $qend) {
                  $qs1 = $loc->start;
                  $qe1 = ($loc->end < $qend) ? $loc->end : $qend;
                
                  $subloc = Bio::Location::Simple->new();
                  $subloc->start($qs1);
                  $subloc->end  ($qe1);
                  $loc2->add_sub_Location($subloc);
              
                  $loc = $locs->[++$i];
                }
                $debug && print STDERR "$subname: \$loc2=\n".Dumper($loc2)."end of \$loc2\n\n";
                last;
              }
            }
          }
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
            $qstart_min = Annotate_Math::adjust_start( $qstart_min, $refcds_loc, $cds, $aln, $aln_q, $aln_h);
            $qstart_min = $qstart_min + $codon_start -1;
            $start = $qstart_min;
        }

        my $qend_min = $reffeat_loc->min_end;
        if (defined($qend_min)) {
            $debug && print STDERR "$subname: \$qend_min=$qend_min\n";
            $qend_min = Annotate_Math::adjust_end( $qend_min, $refcds_loc, $cds, $aln, $aln_q, $aln_h);
            $qend_min = $qend_min + $codon_start -1;
            $end = $qend_min;
        }

        my $qstart_max = $reffeat_loc->max_start;
        if (defined($qstart_max)) {
            $qstart_max = Annotate_Math::adjust_start( $qstart_max, $refcds_loc, $cds, $aln, $aln_q, $aln_h);
            $qstart_max = $qstart_max + $codon_start -1;
            $start = $qstart_max if (not $start<=$qstart_max);
        }

        my $qend_max = $reffeat_loc->max_end;
        if (defined($qend_max)) {
            $qend_max = Annotate_Math::adjust_end( $qend_max, $refcds_loc, $cds, $aln, $aln_q, $aln_h);
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
    my $errcode = {
                    OutsideCDS => 0,
                    partial_mat_peptide => 0,
                    gaps_start => 0,
                    gaps_start => 0,
                  };

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

    $errcode->{partial_mat_peptide} = Annotate_Math::msa_check_loc_start( $loc2, $reffeat_loc, $refcds_loc, $cds_loc, $aln, $aln_q, $aln_h);

    $errcode->{partial_mat_peptide} = $errcode->{partial_mat_peptide} || Annotate_Math::msa_check_loc_end( $loc2, $reffeat_loc, $refcds_loc, $cds_loc, $aln, $aln_q, $aln_h);

    $errcode->{long_internal_gap} = Annotate_Math::msa_check_loc_internal( $loc2, $reffeat_loc, $refcds_loc, $cds_loc, $aln, $aln_q, $aln_h);

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


=head2 convert_range_protein2dna

Takes a CDS, and an array of gaps from alignment,
 return the gaps converted to DNA range, taking into account the start of CDS

=cut

sub convert_range_protein2dna {
    my ($f, $gaps) = @_;

    my $debug = 0 && $debug_all;
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


## sub get_dna_byloc fetches the sequence directly from seq string according to
## the location range. Sequence is always extracted from $seq_obj_seq
sub get_dna_byloc {
    my ($feat_obj, $seq_obj_seq) = @_;

      my $debug = 0 && $debug_all;
      my $subname = 'get_dna_byloc';
      
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

    $debug && print STDERR "get_dna_byloc: \$s = $s\n";
    return $s;
} # sub get_dna_byloc


1;
