package Annotate_Math;

use strict;
use warnings;
use English;
use Carp;
use Data::Dumper;

use version; our $VERSION = qv('1.1.6'); # Feb 05 2013
use File::Temp qw/ tempfile tempdir /;
use Bio::Perl;
use Bio::SeqIO;
use Bio::Seq;
use Bio::AlignIO;
use Bio::Tools::Run::StandAloneBlast;
use IO::String;
use POSIX;

use Annotate_Verify;

my $debug_all = 0;

####//README//####
#
# Annotate_Math contains the core math functions to perform annotation based on MUSCLE
#
#    Authors Chris Larsen, clarsen@vecna.com; Guangyu Sun, gsun@vecna.com
#    December 2010
#
##################

## //EXECUTE// ##

sub setDebugAll {
    my ($debug) = @_;
    $debug_all = $debug;
} # sub setDebugAll


sub adjust_start {
    my ($qstart, $reffeat, $refcds, $cds, $aln, $aln_q, $aln_h, $errcode) = @_;

    return undef if (!$qstart);
    my $debug = 0 && $debug_all;
    my $subn = 'adjust_start';

    my $refcds_loc = $refcds->location;
    my $refcds_loc_str = $refcds->location->to_FTstring;
    my $cds_loc_str = $cds->location->to_FTstring;

    $debug && print STDERR "$subn: \$qstart=$qstart \$refcds=$refcds_loc_str \$cds=$cds_loc_str\n";
#    $debug && print STDERR "$subn: \$qstart=$qstart \$refcds=$refcds_loc_str \$refcds=\n".Dumper($refcds_loc)."\n";
    $errcode->{alnstart} = '' if (!$errcode || !exists($errcode->{alnstart}));

    # Get any codon_start in refseq CDS and reffeat
    my $codon_start = 1;
    if ($refcds->has_tag('codon_start')) {
        $codon_start = [$refcds->get_tag_values('codon_start')];
        $codon_start = $codon_start->[0];
        $debug && print STDERR "$subn: After adjusting refcds \$codon_start=$codon_start\n";
        $codon_start -= [$reffeat->get_tag_values('codon_start')]->[0]-1 if ($reffeat->has_tag('codon_start'));
        $debug && print STDERR "$subn: After adjusting reffeat \$codon_start=$codon_start\n";
    }
    # Get any potential gap caused by split locations in refseq CDS
    my $qstart_gap = Annotate_Math::get_split_gap( $qstart, $refcds_loc);
    $debug && print STDERR "$subn: refcds \$qstart_gap=$qstart_gap\n";
    # For SARS genome DQ497008, there is a stuter at CDS=join(265..13398,13398..21485)
    $debug && print STDERR "$subn: \$qstart=$qstart \$codon_start=$codon_start \$qstart_gap=$qstart_gap\n";

    my $cds_loc = $cds->location;
    # Convert to protain from DNA
    if ($refcds_loc->strand !=-1) {
        $qstart = ($qstart - $refcds_loc->start - $codon_start+1 - $qstart_gap)/3 +1;
    } else {
        $qstart = ($refcds_loc->end - $qstart - $codon_start+1 - $qstart_gap)/3 +1;
    }
    $debug && print STDERR "$subn: After DNA    ->protein: \$qstart=$qstart\n";
#    $debug && print STDERR "$subn: \$aln_q=\n".Dumper($aln_q)."end of $aln_q\n";
    my $qseq = $aln_q->seq;
#    $qseq = substr($qseq, 0 , $qstart-1); # Bad, as qstart can increase in for loop
    $debug && print STDERR "$subn: \$qseq='$qseq' \$qseq=".length($qseq)."\n";
    my $gap_char = $aln->gap_char;
    my $c;
    for (my $n=1; $n<=$qstart && $n<=length($aln_q->seq); $n++) {
        $c = substr($qseq, $n-1, 1);
        $qstart++ if ($c eq $gap_char);
        if (($n-1)%10 ==0) {
#            $debug && print STDERR "\n$subn: \$n=$n \$gap_char=$gap_char ";
        }
#        $debug && print STDERR "\$c=$c \$qstart=$qstart ";
    }
    $debug && print STDERR "\n";

    my ($hstart, $hend);
    $hstart = $qstart;
    my $hseq = $aln_h->seq;
    $hseq = substr($hseq, 0, $hstart-1); # OK, since hstart only decreases in for loop
    $debug && print STDERR "$subn: \$hseq='$hseq' \$hseq=".length($hseq)."\n";
    for (my $n=length($hseq);  $n>=1; $n--) {
        $c = substr($hseq, $n-1, 1);
        $hstart-- if ($c eq $gap_char);
        if (($n-0)%10 ==0 || $n==length($hseq)) {
#            $debug && print STDERR "\n$subn: \$n=$n \$gap_char=$gap_char ";
        }
#        $debug && print STDERR "\$c=$c \$hstart=$hstart ";
    }
    $debug && print STDERR "\n";

    # Get the start of the reffeat in the MSA, this is used later during checking
    $errcode->{alnstart} = $qstart;
    $debug && print STDERR "$subn: \$hstart=$hstart, \$errcode->{alnstart}=$errcode->{alnstart}\n";
    $hseq = $aln_h->seq;
    $hseq = substr($hseq, $qstart-1);
    $debug && print STDERR "$subn: \$hseq='$hseq' \$hseq=".length($hseq)."\n";
    for my $n (1 .. length($hseq)) {
        $c = substr($hseq, $n-1, 1);
        $debug && print STDERR "$subn: \$n='$n' \$c='$c'\n";
        if ($c eq $gap_char) {
            $errcode->{alnstart}++;
        } else {
            last;
        }
    }
    $debug && print STDERR "$subn: \$hstart=$hstart, \$errcode->{alnstart}=$errcode->{alnstart} \$cds=".$cds->strand."\n";

    # Convert back to DNA
    if ($cds->strand!=-1) {
        $hstart = $cds_loc->start + ($hstart-1)*3;
    } else {
        $hstart = $cds_loc->end - ($hstart-1)*3;
        my $codon_start = 1;
        if ($cds->has_tag('codon_start')) {
            $codon_start = [$cds->get_tag_values('codon_start')];
            $codon_start = $codon_start->[0];
            $debug && print STDERR "$subn: cds \$codon_start=$codon_start\n";
        }
        $hstart -= $codon_start-1;
    }
    $debug && print STDERR "$subn: \$hstart=$hstart, \$errcode->{alnstart}=$errcode->{alnstart} \$cds=".$cds->strand."\n";

        # Get any potential gap caused by split locations in refseq CDS
        my $hstart_gap = Annotate_Math::get_split_gap( $hstart, $cds->location);
        $debug && print STDERR "$subn: cds \$hstart_gap=$hstart_gap\n";
        # For SARS genome DQ497008, there is a stuter at CDS=join(265..13398,13398..21485)
#        $qstart += $qstart_gap;
#        $qstart += $qstart_gap if ($qstart_gap % 3 ==0); # Don't remember why used such condition in the beginning
        $debug && print STDERR "$subn: \$hstart=$hstart after checking CDS\n";

    $hstart = $hstart + $hstart_gap; # Include the stutering
        $debug && print STDERR "$subn: \$hstart=$hstart after checking CDS and stutering\n";

        # Get any potential gap caused by split locations in CDS
#        $qstart_gap = Annotate_Math::get_split_gap( $hstart, $cds->location);
#        $debug && print STDERR "$subn: \$qstart_gap=$qstart_gap\n";
#        $hstart += $qstart_gap if ($qstart_gap % 3 ==0);

    return $hstart;
} # sub adjust_start

sub adjust_end {
    my ($qend, $reffeat, $refcds, $cds, $aln, $aln_q, $aln_h, $errcode) = @_;

    return undef if (!$qend);
    my $debug = 0 && $debug_all;
    my $subn = 'adjust_end';

    my $refcds_loc = $refcds->location;
    my $cds_loc_str = $cds->location->to_FTstring;
    $debug && print STDERR "$subn: \$qend=$qend \$refcds_loc=".$refcds_loc->to_FTstring." \$cds=$cds_loc_str\n";
    $errcode->{alnend} = '' if (!$errcode || exists($errcode->{alnend}));

    # Get any codon_start in refseq CDS and reffeat
    my $codon_start = 1;
    if ($refcds->has_tag('codon_start')) {
        $codon_start = [$refcds->get_tag_values('codon_start')];
        $codon_start = $codon_start->[0];
        $debug && print STDERR "$subn: After adjusting refcds \$codon_start=$codon_start\n";
        # Asume the end of reffeat always ends with a complete codon
#        $codon_start -= [$reffeat->get_tag_values('codon_start')]->[0]-1 if ($reffeat->has_tag('codon_start'));
        $debug && print STDERR "$subn: After adjusting reffeat \$codon_start=$codon_start\n";
    }
    # Get any potential gap caused by split locations in refCDS
    my $qend_gap = Annotate_Math::get_split_gap( $qend, $refcds_loc);
    $debug && print STDERR "$subn: refcds \$qend_gap=$qend_gap\n";
    # For SARS genome DQ497008, there is a stuter at CDS=join(265..13398,13398..21485)
#    $qend -= $qend_gap;
#    $qend -= $qend_gap if ($qend_gap % 3 ==0); # Don't remember why used such condition in the beginning
    $debug && print STDERR "$subn: \$qend=$qend after checking refCDS\n";

    my $last_mat_peptide = 0;
    $last_mat_peptide = 1 if ($qend == $refcds_loc->end || $qend == $refcds_loc->end-3);
    $debug && print STDERR "$subn: \$last_mat_peptide=$last_mat_peptide\n";

    my $cds_loc = $cds->location;
    # Convert to protain from DNA
    if ($refcds_loc->strand!=-1) {
        $qend = ($qend+1 - $refcds_loc->start - $codon_start+1 - $qend_gap)/3 -1 +1;
    } else {
        $qend = ($refcds_loc->end - $qend+1 - $codon_start+1 - $qend_gap)/3 -1 +1;
    }
    $debug && print STDERR "$subn: \$qend=$qend\n";
#    $debug && print STDERR "$subn: \$aln_q=\n".Dumper($aln_q)."end of $aln_q\n";
    my $qseq = $aln_q->seq;
#    $qseq = substr($qseq, 0, $qend); # Bad, since qend can increase in for loop
    $debug && print STDERR "$subn: \$qseq='$qseq' \$qseq=".length($qseq)."\n";
    my $gap_char = $aln->gap_char;
    my $c;
    for (my $n=1; $n<=$qend && $n<=length($aln_q->seq); $n++) {
        $c = substr($qseq, $n-1, 1);
        $qend++ if ($c eq $gap_char);
        if (($n-1)%10 ==0) {
#            $debug && print STDERR "\n$subn: \$n=$n \$gap_char=$gap_char ";
        }
#        $debug && print STDERR "\$c=$c \$qend=$qend ";
    }
    $debug && print STDERR "\n";
    $debug && print STDERR "$subn: \$qend=$qend\n";

    # If this is the last mat_peptide within refseq CDS, see if there is any trailing short sequence in target.
    # If so, and the length is 5 AA or less, include it in the annotation.
    if ($last_mat_peptide) {
        if ($qend<length($qseq) && abs(length($qseq)-$qend)<=5) {
            my $old = $qend;
            $qend = length($qseq);
            print STDERR "$subn: \$qend=$old changed to $qend to include the trailing sequence(".substr($qseq,$qend-10).") after last mat_peptide\n";
        }
    }

    my ($hstart, $hend);
    $hend = $qend;
    my $hseq = $aln_h->seq;
    $hseq = $aln_h->seq;
    $hseq = substr($hseq, 0, $hend); # OK, as hend only decreases in for loop
    $debug && print STDERR "$subn: \$hseq='$hseq' \$hseq=".length($hseq)."\n";
    for (my $n=length($hseq);  $n>=1; $n--) {
        $c = substr($hseq, $n-1, 1);
        $hend-- if ($c eq $gap_char);
        if (($n-0)%10 ==0 || $n==length($hseq)) {
#            $debug && print STDERR "\n$subn: \$n=$n \$gap_char=$gap_char ";
        }
#        $debug && print STDERR "\$c=$c \$hend=$hend ";
    }
    $debug && print STDERR "\n";

    # Get the end of the reffeat in the MSA, this is used later during checking
    $errcode->{alnend} = $qend;
    $debug && print STDERR "$subn: \$hend=$hend, \$errcode->{alnend}=$errcode->{alnend}\n";
    $hseq = $aln_h->seq;
    $hseq = substr($hseq, 0, $qend);
    $debug && print STDERR "$subn: \$hseq='$hseq' \$hseq=".length($hseq)."\n";
    for (my $n=length($hseq);  $n>=1; $n--) {
        $c = substr($hseq, $n-1, 1);
        if ($c eq $gap_char) {
            $errcode->{alnend}--;
        } else {
            last;
        }
    }
    $debug && print STDERR "$subn: \$hend=$hend, \$errcode->{alnend}=$errcode->{alnend} after checking gap at end\n";

    # Convert back to DNA
    if ($cds->strand!=-1) {
        $hend = $cds_loc->start + ($hend-1)*3 +2;
    } else {
#        $hend = $cds_loc->end - ($hend-1)*3 +2;
        $hend = $cds_loc->end - ($hend-1)*3 -2;
    }
    $debug && print STDERR "$subn: \$cds_loc=".$cds_loc->to_FTstring." \$hend=$hend\n";

        # Get any potential gap caused by split locations in CDS
        my $hend_gap = Annotate_Math::get_split_gap( $hend, $cds->location);
        $debug && print STDERR "$subn: cds \$hend=$hend \$hend_gap=$hend_gap\n";
        # For SARS genome DQ497008, there is a stuter at CDS=join(265..13398,13398..21485)
#        $hend += $hend_gap;
#        $qend += $hend_gap if ($hend_gap % 3 ==0); # Don't remember why used such condition in the beginning
        $debug && print STDERR "$subn: \$hend=$hend after checking CDS\n";

    $hend = $hend + $hend_gap; # Include the stuttering
    $debug && print STDERR "$subn: \$cds_loc->start=".$cds_loc->start." \$hend=$hend\n";

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
    if ($hend-1 == $cds_loc->end) {
            my $tail = $cds->seq->seq;
            $tail = substr($tail, length($tail)-2, 2);
            $debug && print STDERR "$subn: \$tail=$tail\n";
            $hend-- if (exists($good_tail{$tail}) && $good_tail{$tail});
    }
    $codon_start = 1;
    if ($cds->has_tag('codon_start')) {
        $codon_start = [$cds->get_tag_values('codon_start')];
        $codon_start = $codon_start->[0];
    }
    $debug && print STDERR "$subn: \$codon_start=$codon_start \$hend=$hend\n";
    if ($cds->strand!=-1) {
        if ($codon_start>1 && ($cds_loc->end - $hend)<=2) {
            $hend -= $codon_start -1;
            $debug && print STDERR "$subn: \$codon_start=$codon_start \$hend=$hend\n";
        }
    } else {
            $hend -= $codon_start -1;
    }

        # Get any potential gap caused by split locations in CDS
#        $qend_gap = Annotate_Math::get_split_gap( $hend, $cds->location);
#        $debug && print STDERR "$subn: \$qend_gap=$qend_gap\n";
#        $hend += $qend_gap if ($qend_gap % 3 ==0);

    $debug && print STDERR "$subn: final value \$hend=$hend\n";
#    $hend += 1 if ($hend+1 == $cds->location->end); # Handles incomplete CDS w/ untranslated partial codon
#    $debug && print STDERR "$subn: \$hend=$hend\n";
    return $hend;

} # sub adjust_end


=head2 get_split_gap

Takes a position in the DNA sequence, and count how many gaps (U34999:join(44..5695,5699..7540)) prior to it
 Return the number of gaps in DNA unit

=cut

sub get_split_gap {
    my ($pos, $cds_loc) = @_;

    my $debug = 0 && $debug_all;
    my $subn = 'get_split_gap';

        # For split locations, count the gap between sublocations prior to the requested position
        my $split = 0;
        if ($cds_loc->isa('Bio::Location::Split')) {
            my $locs = [ $cds_loc->sub_Location ];
            for my $i (1 .. $#{$locs}) {
               $debug && print STDERR "$subn: \$pos=$pos \$locs->[$i]=".$locs->[$i]->to_FTstring." \$split=$split\n";
               if ($pos < $locs->[$i]->start) {
                 $debug && print STDERR "$subn: 1 pos=$pos start=".$locs->[$i]->to_FTstring."\n";
                 next;
               } else {
                 $debug && print STDERR "$subn: 2 pos=$pos start=".$locs->[$i]->to_FTstring."\n";
                 $split += $locs->[$i]->start - $locs->[$i-1]->end -1;
               }
               $debug && print STDERR "$subn: \$pos=$pos \$locs->[$i]=".$locs->[$i]->to_FTstring." \$split=$split\n";
            }
        }

    return $split;
} # sub get_split_gap

=head2 msa_get_feature_loc

Takes reffeat, gaps in the reffeat, and gaps (in DNA coordinate) in the target
 calculate the location of the feature in target genome
 Return a Bio::Location object for the new annotation
=cut

sub msa_get_feature_loc {
    my ($reffeat, $refcds, $cds, $aln, $aln_q, $aln_h) = @_;

    my $debug = 1 && $debug_all;
    my $subn = 'msa_get_feature_loc';

    my $reffeat_loc = $reffeat->location;
    my $cds_loc = $cds->location;
    my $cds_loc_str = $cds->location->to_FTstring;
    my $refcds_loc = $refcds->location;
    my $loc2 = undef;
    my $errcode = { alnstart=>'', alnend=>'',};

    $debug && print STDERR "$subn: \$aln=\n".Dumper($aln)."end of \$aln\n\n";
    $debug && print STDERR "$subn: \$reffeat_loc=\n".Dumper($reffeat_loc)."end of \$reffeat_loc\n\n";
    $debug && print STDERR "\n$subn: starting reffeat=".$reffeat_loc->to_FTstring." strand=".$reffeat_loc->strand."\n";
    if ($reffeat_loc->isa('Bio::Location::Simple') || $reffeat_loc->isa('Bio::Location::Split')) {
        my $qstart = $reffeat_loc->start;
        my $qend = $reffeat_loc->end;
        $debug && print STDERR "$subn: \$qstart=$qstart \$qend=$qend\n";

        my $codon_start = 1;
        if ($cds->has_tag('codon_start')) {
            $codon_start = [$cds->get_tag_values('codon_start')];
            $codon_start = $codon_start->[0];
        }
        if ($reffeat_loc->strand !=-1 && $cds->strand!=-1) {
            $qstart = Annotate_Math::adjust_start( $qstart, $reffeat, $refcds, $cds, $aln, $aln_q, $aln_h, $errcode);
            $debug && print STDERR "$subn: \$qstart=$qstart \$qend=$qend\n";
            $qend = Annotate_Math::adjust_end( $qend, $reffeat, $refcds, $cds, $aln, $aln_q, $aln_h, $errcode);
            $qstart = $qstart + $codon_start -1;
            $qend = $qend + $codon_start -1;
            $qend += 1 if ($qend+1 == $cds->location->end); # Handles incomplete CDS w/ untranslated partial codon
            $qend += 2 if ($qend+2 == $cds->location->end); # Handles incomplete CDS w/ untranslated partial codon
            $debug && print STDERR "$subn: After codon_start \$qstart=$qstart \$qend=$qend for \$cds=$cds_loc_str\n";
        } elsif ($reffeat_loc->strand !=-1 && $cds->strand==-1) {
            my $t = Annotate_Math::adjust_start( $qstart, $reffeat, $refcds, $cds, $aln, $aln_q, $aln_h, $errcode);
            $debug && print STDERR "$subn: \$reffeat_loc=".$reffeat_loc->to_FTstring." \$qstart=$qstart \$qend=$qend \$cds=".$cds->strand."\n";
            $qstart = Annotate_Math::adjust_end( $qend, $reffeat, $refcds, $cds, $aln, $aln_q, $aln_h, $errcode);
            $qend = $t;
            $debug && print STDERR "$subn: \$reffeat_loc=".$reffeat_loc->to_FTstring." \$qstart=$qstart \$qend=$qend \$cds=".$cds->strand."\n";
#            $qstart = $qstart + $codon_start -1;
#            $qend = $qend + $codon_start -1;
            $qstart -= 1 if ($qstart-1 == $cds->location->start); # Handles incomplete CDS w/ untranslated partial codon
            $qstart -= 2 if ($qstart-2 == $cds->location->start); # Handles incomplete CDS w/ untranslated partial codon
            $debug && print STDERR "$subn: After codon_start \$qstart=$qstart \$qend=$qend for \$cds=$cds_loc_str\n";
        } elsif ($reffeat_loc->strand ==-1 && $cds->strand==-1) {
            $qend = Annotate_Math::adjust_start( $qend, $reffeat, $refcds, $cds, $aln, $aln_q, $aln_h, $errcode);
            $debug && print STDERR "$subn: \$qstart=$qstart \$qend=$qend\n";
            $qstart = Annotate_Math::adjust_end( $qstart, $reffeat, $refcds, $cds, $aln, $aln_q, $aln_h, $errcode);
            $qstart = $qstart + $codon_start -1;
            $qend = $qend + $codon_start -1;
            $qstart -= 1 if ($qstart-1 == $cds->location->start); # Handles incomplete CDS w/ untranslated partial codon
            $qstart -= 2 if ($qstart-2 == $cds->location->start); # Handles incomplete CDS w/ untranslated partial codon
            $debug && print STDERR "$subn: After codon_start \$qstart=$qstart \$qend=$qend for \$cds=$cds_loc_str\n";
        } elsif ($reffeat_loc->strand ==-1 && $cds->strand!=-1) {
            my $t = Annotate_Math::adjust_start( $qend, $reffeat, $refcds, $cds, $aln, $aln_q, $aln_h, $errcode);
            $debug && print STDERR "$subn: \$qstart=$qstart \$qend=$qend\n";
            $qend = Annotate_Math::adjust_end( $qstart, $reffeat, $refcds, $cds, $aln, $aln_q, $aln_h, $errcode);
            $qstart = $t;
            $qstart = $qstart + $codon_start -1;
            $qend = $qend + $codon_start -1;
            $qend += 1 if ($qend+1 == $cds->location->end); # Handles incomplete CDS w/ untranslated partial codon
            $qend += 2 if ($qend+2 == $cds->location->end); # Handles incomplete CDS w/ untranslated partial codon
            $debug && print STDERR "$subn: After codon_start \$qstart=$qstart \$qend=$qend for \$cds=$cds_loc_str\n";
        } else {
            print STDERR "\n$subn: ERROR: Wired that we end up here, please contact author for any fix.\n";
            print STDERR "\n$subn: ERROR: acc=".$cds->seq->accession_number." refseq=".$refcds->seq->accession_number.".\n\n";
            exit(1);
        }
        $debug && print STDERR "$subn: After checking MSA \$qstart=$qstart \$qend=$qend\n";
        $debug && print STDERR "$subn: \$errcode=\n".Dumper($errcode)."end of \$errcode\n\n";


        # Get any potential gap caused by split locations in CDS
#        my $qstart_gap = Annotate_Math::get_split_gap( $qstart, $cds->location);
#        my $qend_gap = Annotate_Math::get_split_gap( $qend, $cds->location);
#        $debug && print STDERR "$subn: \$qstart_gap=$qstart_gap \$qend_gap=$qend_gap\n";
#        $qstart += $qstart_gap if ($qstart_gap % 3 ==0);
#        $qend += $qend_gap if ($qend_gap % 3 ==0);

        $debug && print STDERR "$subn: \$qstart=$qstart \$qend=$qend for \$cds=$cds_loc_str\n";

#        if ($qstart<=$qend) {
        if ($qstart<$qend) {
          if ($cds_loc->isa('Bio::Location::Simple') || $cds_loc->isa('Bio::Location::Fuzzy')) {
            $loc2 = Bio::Location::Simple->new();
            $loc2->start($qstart);
            $loc2->end  ($qend);
#            $loc2->strand(-1) if ($cds_loc->strand() == -1);
          } elsif ( $cds_loc->isa('Bio::Location::Split')) {
            my $locs = [ $cds->location->sub_Location ];
            $debug && print STDERR "$subn: \$locs=\n".Dumper($cds->location)."end of \$locs\n\n";
            $debug && print STDERR "$subn: \$locs=\n".Dumper($locs)."end of \$locs\n\n";
            for (my $i=0; $i<=$#{$locs}; $i++ ) {
              my $loc = $locs->[$i];
              $debug && print STDERR "$subn: \$i=$i \$qstart=$qstart, \$loc->start=".$loc->start."\n";
              $debug && print STDERR "$subn: \$i=$i   \$qend=$qend,   \$loc->end=".$loc->end."\n";
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
                $debug && print STDERR "$subn: \$loc=\n".Dumper($loc)."end of \$loc\n\n";
                while ($i<=$#{$locs} && $loc->start < $qend) {
                  $qs1 = $loc->start;
                  $qe1 = ($loc->end < $qend) ? $loc->end : $qend;
                
                  $subloc = Bio::Location::Simple->new();
                  $subloc->start($qs1);
                  $subloc->end  ($qe1);
                  $loc2->add_sub_Location($subloc);
              
                  $loc = $locs->[++$i];
                }
                $debug && print STDERR "$subn: \$i=$i \$loc2=\n".Dumper($loc2)."end of \$loc2\n\n";
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
            $debug && print STDERR "$subn: \$qstart_min=$qstart_min\n";
            $qstart_min = Annotate_Math::adjust_start( $qstart_min, $reffeat, $refcds, $cds, $aln, $aln_q, $aln_h, $errcode);
            $qstart_min = $qstart_min + $codon_start -1;
            $start = $qstart_min;
            $debug && print STDERR "$subn: Fuzzy \$qstart_min=$qstart_min \$start=$start \$end=$end\n";
            $debug && print STDERR "$subn: \$errcode=\n".Dumper($errcode)."end of \$errcode\n\n";
        }

        my $qend_min = $reffeat_loc->min_end;
        if (defined($qend_min)) {
            $debug && print STDERR "$subn: \$qend_min=$qend_min\n";
            $qend_min = Annotate_Math::adjust_end( $qend_min, $reffeat, $refcds, $cds, $aln, $aln_q, $aln_h, $errcode);
            $qend_min = $qend_min + $codon_start -1;
            $qend_min += 1 if ($qend_min+1 == $cds->location->end); # Handles incomplete CDS w/ untranslated partial codon
            $qend_min += 2 if ($qend_min+2 == $cds->location->end); # Handles incomplete CDS w/ untranslated partial codon
            $end = $qend_min;
            $debug && print STDERR "$subn: Fuzzy \$qend_min=$qend_min \$start=$start \$end=$end\n";
            $debug && print STDERR "$subn: \$errcode=\n".Dumper($errcode)."end of \$errcode\n\n";
        }

        my $qstart_max = $reffeat_loc->max_start;
        if (defined($qstart_max)) {
            $qstart_max = Annotate_Math::adjust_start( $qstart_max, $reffeat, $refcds, $cds, $aln, $aln_q, $aln_h, $errcode);
            $qstart_max = $qstart_max + $codon_start -1;
            $start = $qstart_max if (not $start<=$qstart_max);
            $debug && print STDERR "$subn: Fuzzy \$qstart_max=$qstart_max \$start=$start \$end=$end\n";
            $debug && print STDERR "$subn: \$errcode=\n".Dumper($errcode)."end of \$errcode\n\n";
        }

        my $qend_max = $reffeat_loc->max_end;
        if (defined($qend_max)) {
            $qend_max = Annotate_Math::adjust_end( $qend_max, $reffeat, $refcds, $cds, $aln, $aln_q, $aln_h, $errcode);
            $qend_max = $qend_max + $codon_start -1;
            $qend_max += 1 if ($qend_max+1 == $cds->location->end); # Handles incomplete CDS w/ untranslated partial codon
            $qend_max += 2 if ($qend_max+2 == $cds->location->end); # Handles incomplete CDS w/ untranslated partial codon
            $end = $qend_max if ($end<=$qend_max);
            $debug && print STDERR "$subn: Fuzzy \$qend_max=$qend_max \$start=$start \$end=$end\n";
            $debug && print STDERR "$subn: \$errcode=\n".Dumper($errcode)."end of \$errcode\n\n";
        }

#        if ($start<=$end) {
        if ($start<$end) {
            # start must be less than end
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
    if (!$loc2) {;
        $debug && print STDERR "$subn: \$loc2 is undef, return before checking\n";
        return ($loc2, $errcode);
    }
    # check the validity of $loc2, this ensures the new annotated feature lies within the cds,
    my $err1 = Annotate_Math::msa_check_loc( $loc2, $reffeat_loc, $refcds_loc, $cds_loc, $aln, $aln_q, $aln_h, $errcode);
    for my $k (keys %$err1) {
#        $errcode->{$k} = $err1->{$k};
    }

if ( 0 ) {
    # Found a new problem, eg AF126284 w/ refseq=NC_003900. --10/12/2012
    # The last mat_peptide=join(9884..10012,10012..10029) w/ frameshift doesn't have corresponding CDs
    # in AF126284. This causes CDS=7628..11362 to be annotated, leading to wrong result.
    # Added check to ensure resulting mat_peptide must be split if refseq is split when length is similar
    $debug && $loc2 && print STDERR "$subn: \$loc2=".$loc2->to_FTstring." compares to reffeat=".$reffeat_loc->to_FTstring."\n";
  if ( 1 ) {
    if ( $loc2 && (($loc2->end - $loc2->start)<90) ) {
        if ( $reffeat_loc->isa('Bio::Location::Split') && !$loc2->isa('Bio::Location::Split')) { # AY692454
            checkSplit($loc2, $reffeat_loc, $refcds_loc, $cds_loc, $aln, $aln_q, $aln_h, $errcode);
            print STDERR "$subn: Changing \$loc2=".$loc2->to_FTstring." to undef as reffeat=".$reffeat_loc->to_FTstring." is Split\n";
            $loc2 = undef;
#        exit;
        } elsif ( !$reffeat_loc->isa('Bio::Location::Split') && $loc2->isa('Bio::Location::Split')) {
        # Found a new problem, eg EU714028 w/ refseq=NC_001451. --10/17/2012
        # Added check to ensure resulting mat_peptide must be non-split if refseq is non-split when length is similar
            print STDERR "$subn: Changing \$loc2=".$loc2->to_FTstring." to undef as reffeat=".$reffeat_loc->to_FTstring." is non-Split\n";
            $loc2 = undef;
#        exit;
        }
    }
  } else {
  }
} else {
            Annotate_Math::checkSplitLocation( $loc2, $cds, $reffeat_loc, $refcds_loc, $cds_loc, $aln, $aln_q, $aln_h, $errcode);
}

    # Print out the alignment for $loc2
    my @qs = ('', '', '');
    my @hs = ('', '', '');
    my $tail_length = 5;
    if ($errcode->{alnstart}>1) {
        $qs[0] = substr($aln_q->seq, $errcode->{alnstart}-1-$tail_length, $tail_length);
        $hs[0] = substr($aln_h->seq, $errcode->{alnstart}-1-$tail_length, $tail_length);
    }
    $qs[1] = substr($aln_q->seq, $errcode->{alnstart}-1, $errcode->{alnend}-$errcode->{alnstart});
    $hs[1] = substr($aln_q->seq, $errcode->{alnstart}-1, $errcode->{alnend}-$errcode->{alnstart});
    if ($errcode->{alnend}<length($aln_q->seq)) {
        $qs[2] = substr($aln_q->seq, $errcode->{alnend}, $tail_length);
        $hs[2] = substr($aln_h->seq, $errcode->{alnend}, $tail_length);
    }
    $debug && print STDERR "$subn: \@qs='@qs'\n";
    $debug && print STDERR "$subn: \@hs='@hs'\n";

    $debug && print STDERR "$subn: \$loc2=\n".Dumper($loc2)."end of \$loc2\n";
    $debug && print STDERR "$subn: \$errcode=\n".Dumper($errcode)."end of \$errcode\n";

    return ($loc2, $errcode);

} # sub msa_get_feature_loc


=head2 checkSplitLocation
To match the new annotation with reference mat_peptide in their split-ness.
=cut

sub checkSplitLocation {
    my ($loc2, $cds, $reffeat_loc, $refcds_loc, $cds_loc, $aln, $aln_q, $aln_h, $errcode) = @_;

    my $debug = 0 && $debug_all;
    my $subn = 'checkSplitLocation';
    # error code: OutsideCDS; partial_mat_peptide; GapAtCleaveSite
#    my $errcode = {
#                    OutsideCDS => 0,
#                    partial_mat_peptide => 0,
#                    long_internal_gap => 0,   # used to throw out any new location w/ long gap
#                    has_internal_gap => 0,    # used to throw out any new location not confirming with reffeat in term of split location
#                    gaps_start => 0,
#                  };
    $errcode->{SplitMisatch} = 0 if (!$errcode || !exists($errcode->{SplitMisatch}));

    return ($errcode) if (!$loc2);

    # Found a new problem, eg AF126284 w/ refseq=NC_003900. --10/12/2012
    # The last mat_peptide=join(9884..10012,10012..10029) w/ frameshift doesn't have corresponding CDs
    # in AF126284. This causes CDS=7628..11362 to be annotated, leading to wrong result.
    # Added check to ensure resulting mat_peptide must be split if refseq is split when length is similar
    # test case: AF201929.1 (GI:6625759) shouldn't have nsp12=13424..16198
    # test case: AF220295.1 (GI:17529670) shouldn't have nsp9=13332..16100
    # test case: AY692454.1 (GI:56555212) should have nsp9=12459..15131
    # test case: EF424619.1 (GI:145208932) should have nsp9=13325..16084
    # test case: FJ884686.1 (GI:226693969) shouldn't have nsp12=13566..16340
    # test case: JF330899.1 (GI:324963481) should have nsp9=12331..15077
    # test case: AF126284.1 (GI:4240567) shouldn't have TF=9884..10030
    # test case: DQ241304.1 (GI:78499741) shouldn't have TF=9770..10069
    # test case: U34999 (GI:1144527) should have nsp3=join(4031..5695,5699..5716) as the gap at is stop codon and thus treated as read-through
    $debug && print STDERR "$subn: \$loc2=".$loc2->to_FTstring." compares to reffeat=".$reffeat_loc->to_FTstring."\n";
    my $mismatch = 0;
    if ( $reffeat_loc->isa('Bio::Location::Split') && !$loc2->isa('Bio::Location::Split')) { # AY692454
        $debug && print STDERR "$subn: \$loc2=".$loc2->to_FTstring." is '".ref($loc2)."'\n";
        $debug && print STDERR "$subn: reffeat=".$reffeat_loc->to_FTstring." is '".ref($reffeat_loc)."'\n";

        $mismatch = 1;
        my $len = 0;
        my $loc2_len = $loc2->end - $loc2->start +1;
        my $locs = [ $reffeat_loc->sub_Location ];
        for (my $i=0; $i<=$#{$locs}; $i++ ) {
            my $loc = $locs->[$i];
            if ($i==0) {
                $len = $loc->end - $loc->start +1;
            } else {
                $len = $loc->end - $loc->start +1;
                
            }
            $mismatch = 0 if ($loc2_len <= $len);
            $debug && print STDERR "$subn: \$i=$i \$loc2_len=$loc2_len \$len=$len \$mismatch=$mismatch\n";
            
        }
#        exit;

    } elsif ( !$reffeat_loc->isa('Bio::Location::Split') && $loc2->isa('Bio::Location::Split')) {
        # Found a new problem, eg EU714028 w/ refseq=NC_001451. --10/17/2012
        # Added check to ensure resulting mat_peptide must be non-split if refseq is non-split when length is similar
        $debug && print STDERR "$subn: \$loc2=".$loc2->to_FTstring." is '".ref($loc2)."'\n";
        $debug && print STDERR "$subn: reffeat=".$reffeat_loc->to_FTstring." is '".ref($reffeat_loc)."'\n";
        $mismatch = 0;
        my $len = $reffeat_loc->end - $reffeat_loc->start +1;
        my $loc2_len = 0;
        my $locs = [ $loc2->sub_Location ];
        for (my $i=0; $i<$#{$locs}; $i++ ) {
            my $loc = $locs->[$i];
if (($locs->[$i+1]->start-1 - $loc->end)!=3) {
            $debug && print STDERR "$subn: Changing \$loc2=".$loc2->to_FTstring." to undef as reffeat=".$reffeat_loc->to_FTstring." is non-Split\n";
            $mismatch = 1;
} else {
            my $loc1 = Bio::Location::Simple->new(
                                                   -start => $loc->end+1,
                                                   -end   => $locs->[$i+1]->start-1,
                                                 );
#            $debug && print STDERR "$subn: \$loc1=".Dumper($loc1)."'\n";
            my $feat;
            $feat = Bio::SeqFeature::Generic->new( -location => $loc1, );
            my $s = Annotate_Math::get_dna_byloc( $feat, $cds->{_gsf_seq}->seq);
#            $debug && print STDERR "$subn: \$s='$s' \$loc1=".Dumper($loc1)."'\n";
            # because the annotation $loc2 is from #cds, and by reference mat_peptide, the gap will be treated as read-through
            if ($s !~ /(TAG|TGA|TGG)/i) {
                print STDERR "$subn: Gap between ".$loc->to_FTstring."-".$locs->[$i+1]->to_FTstring." is '$s', not a stop codon\n";
                $mismatch = 1;
            } else {
                print STDERR "$subn: Gap between ".$loc->to_FTstring."-".$locs->[$i+1]->to_FTstring." is '$s', is a stop codon, treated as read-through\n";
            }
}
        }
    }

    $debug && print STDERR "$subn: \$mismatch=$mismatch \$loc2=".$loc2->to_FTstring." as reffeat=".$reffeat_loc->to_FTstring."\n";
    if ($mismatch) {
        print STDERR "$subn: \$loc2=".$loc2->to_FTstring." is '".ref($loc2)."'\n";
        print STDERR "$subn: reffeat=".$reffeat_loc->to_FTstring." is '".ref($reffeat_loc)."'\n";
        print STDERR "$subn: Changing \$loc2=".$loc2->to_FTstring." to undef as reffeat=".$reffeat_loc->to_FTstring."\n";
        #$loc2 = undef;
        $_[0] = undef;
        $errcode->{SplitMisatch} = $mismatch;
    }

    return ($errcode->{SplitMisatch});
} # sub checkSplitLocation


=head2 msa_check_loc
Takes reffeat, gaps in the reffeat, and gaps (in DNA coordinate) in the target
 calculate the location of the feature in target genome
 Return a Bio::Location object for the new annotation
=cut

sub msa_check_loc {
    my ($loc2, $reffeat_loc, $refcds_loc, $cds_loc, $aln, $aln_q, $aln_h, $errcode) = @_;

    my $debug = 0 || $debug_all;
    my $subn = 'msa_check_loc';
    # error code: OutsideCDS; partial_mat_peptide; GapAtCleaveSite
#    my $errcode = {
#                    OutsideCDS => 0,
#                    partial_mat_peptide => 0,
#                    long_internal_gap => 0,   # used to throw out any new location w/ long gap
#                    has_internal_gap => 0,    # used to throw out any new location not confirming with reffeat in term of split location
#                    gaps_start => 0,
#                  };
    $errcode->{OutsideCDS}          = 0 if (!$errcode || !exists($errcode->{OutsideCDS}));
    $errcode->{partial_mat_peptide} = 0 if (!$errcode || !exists($errcode->{partial_mat_peptide}));
    $errcode->{long_internal_gap}   = 0 if (!$errcode || !exists($errcode->{long_internal_gap}));
    $errcode->{has_internal_gap}    = 0 if (!$errcode || !exists($errcode->{has_internal_gap}));

    return ($errcode) if (!$loc2);

#    my ($qstart, $qend);
    # take qstart & qend from $errcode
#    if ( 1 ) {
#        $qstart = $errcode->{alnstart};
#        $qend   = $errcode->{alnend  };
#    }
#    my $qseq = $aln_q->seq;
#    $qseq = substr($qseq, $qstart-1, $qend-$qstart+1); # Any gap at the beginning of CDS, plus the 1st AA of current feat.
#    $debug && print STDERR "$subn: \$loc2=$loc2_str \$qseq='$qseq' \$qseq=".length($qseq)."\n";
#    $qseq =~ s/^[$gap_char]*//; # remove any gap at start
#    $qseq =~ s/[$gap_char]*$//; # remove any gap at end
#

    my $msg;
    my $max_gap_end = 3;
    my $max_gap_middle = 3;

    # See if $loc2 is outside of $cds_loc
    if ($loc2->start < $cds_loc->start) {  # if $loc2 is outside of $cds_loc
            $debug && print STDERR "$subn: a \$loc2=".$loc2->start.' CDS='.$cds_loc->start."\n";
    }
    if ($loc2->end <= $cds_loc->start) {  # if $loc2 is outside of $cds_loc
            $debug && print STDERR "$subn: b \$loc2=".$loc2->end.' CDS='.$cds_loc->start."\n";
    }
    if ($loc2->start >= $cds_loc->end) {  # if $loc2 is outside of $cds_loc
            $debug && print STDERR "$subn: c \$loc2=".$loc2->start.' CDS='.$cds_loc->end."\n";
    }
    if ($loc2->end > $cds_loc->end) {  # if $loc2 is outside of $cds_loc
            $debug && print STDERR "$subn: d \$loc2=".$loc2->end.' CDS='.$cds_loc->end."\n";
    }
    if ($loc2->start < $cds_loc->start
         || $loc2->end <= $cds_loc->start
         || $loc2->start >= $cds_loc->end
         || $loc2->end > $cds_loc->end
       ) {  # if $loc2 is outside of $cds_loc
            print STDERR "$subn: e \$loc2=".$loc2->to_FTstring." doesn't match \$cds=".$cds_loc->to_FTstring."\n";
    # error code: OutsideCDS; partial_mat_peptide; GapAtCleaveSite
            $errcode->{OutsideCDS} = 1;
            return ($errcode);
    }

#    $errcode->{partial_mat_peptide} = Annotate_Math::msa_check_loc_start( $loc2, $reffeat_loc, $refcds_loc, $cds_loc, $aln, $aln_q, $aln_h);
#    $errcode->{partial_mat_peptide} = $errcode->{partial_mat_peptide} || Annotate_Math::msa_check_loc_end( $loc2, $reffeat_loc, $refcds_loc, $cds_loc, $aln, $aln_q, $aln_h);
    Annotate_Math::msa_check_loc_start( $loc2, $reffeat_loc, $refcds_loc, $cds_loc, $aln, $aln_q, $aln_h, $errcode);

    Annotate_Math::msa_check_loc_end( $loc2, $reffeat_loc, $refcds_loc, $cds_loc, $aln, $aln_q, $aln_h, $errcode);

#    $errcode->{long_internal_gap} = Annotate_Math::msa_check_loc_internal( $loc2, $reffeat_loc, $refcds_loc, $cds_loc, $aln, $aln_q, $aln_h);
    Annotate_Math::msa_check_internal( $loc2, $reffeat_loc, $refcds_loc, $cds_loc, $aln, $aln_q, $aln_h, $errcode);

    $debug && print STDERR "$subn: \$errcode=\n".Dumper($errcode)."\n";

    $debug && print STDERR "$subn: subroutine finished\n";
    return ($errcode);

} # sub msa_check_loc

# V1.1.2
sub msa_check_internal {
#    my ($loc2, $reffeat_loc, $refcds_loc, $cds_loc, $aln, $aln_q, $aln_h) = @_;
    my ($loc2, $reffeat_loc, $refcds_loc, $cds_loc, $aln, $aln_q, $aln_h, $errcode) = @_;

    my $debug = 0 || $debug_all;
    my $subname = 'msa_check_internal';
#    my $errcode = { long_internal_gap => 0 };
    $errcode->{long_internal_gap} = 0 if (!$errcode || !exists($errcode->{long_internal_gap}));
    $errcode->{has_internal_gap} = 0 if (!$errcode || !exists($errcode->{has_internal_gap}));

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
        $debug && print STDERR "$subname: \$qseq='$qseq' \$qseq=".length($qseq)."\n";
        $qseq = substr($qseq, $qstart, $qend-$qstart+1); # Any gap at the beginning of CDS, plus the 1st AA of current feat.
        $debug && print STDERR "$subname: \$qseq='$qseq' \$qseq=".length($qseq)."\n";
        $qseq =~ s/^[$gap_char]*//; # remove any gap at start
        $qseq =~ s/[$gap_char]*$//; # remove any gap at end
        my $hseq = $aln_h->seq;
        $debug && print STDERR "$subname: \$hseq='$hseq' \$hseq=".length($hseq)."\n";
        $hseq = substr($hseq, $qstart, $qend-$qstart+1); # Any gap at the beginning of CDS, plus the 1st AA of current feat.
        $debug && print STDERR "$subname: \$hseq='$hseq' \$hseq=".length($hseq)."\n";
        $hseq =~ s/^[$gap_char]*//; # remove any gap at start
        $hseq =~ s/[$gap_char]*$//; # remove any gap at end
        my $count = 0;
        my @chars = split(//, $qseq);
        # Count the gaps in querry
        foreach my $c (@chars) {
            $count++ if ($c eq $gap_char);
        }
        $debug && print STDERR "$subname: \$hseq='$hseq' \$hseq=".length($hseq)." \$count=$count\n";
        if ($count > 0.333*length($qseq)) {
            print STDERR "$subname: ERROR: refcds=".$refcds_loc->to_FTstring." refmat=".$reffeat_loc->to_FTstring." cds=".$cds_loc->to_FTstring." \$qstart=$qstart \$qend=$qend \$gap_char=$gap_char\n";
            print STDERR "$subname: ERROR: \$qseq='$qseq' \$qseq=".length($qseq)." \$count=$count\n";
            $errcode->{long_internal_gap} = 1;
        }
        $count = 0;
        @chars = split(//, $hseq);
        # Count the gaps in hit
        foreach my $c (@chars) {
            $count++ if ($c eq $gap_char);
        }
        $debug && print STDERR "$subname: \$hseq='$hseq' \$hseq=".length($hseq)." \$count=$count\n";
        if ($count > 0.333*length($hseq)) {
            print STDERR "$subname: ERROR: refcds=".$refcds_loc->to_FTstring." refmat=".$reffeat_loc->to_FTstring." cds=".$cds_loc->to_FTstring." \$qstart=$qstart \$qend=$qend \$gap_char=$gap_char\n";
            print STDERR "$subname: ERROR: \$hseq='$hseq' \$hseq=".length($hseq)." \$count=$count\n";
            $errcode->{long_internal_gap} = 1;
        }
    $debug && print STDERR "$subname: \$errcode=\n".Dumper($errcode)."\n";

    return ($errcode->{long_internal_gap});

} # sub msa_check_internal

=head2
 sub msa_check_internal2 checks if a new $loc2 has long internal gap. --11/07/2012
 After much consideration, the rule for sequence with internal gap is to include as mat_peptide following:
 1. ungapped run of alignment has to be at least 10AA (eg. AB182800 has short sequence without gap), and
 2. length of gaped sequence is less than 33% of total length. Gapped sequence is defined as shorter than 10AA
=cut
# V1.1.7
sub msa_check_internal2 {
    my ($loc2, $reffeat_loc, $refcds_loc, $cds_loc, $aln, $aln_q, $aln_h, $errcode) = @_;

    my $debug = 0 || $debug_all;
    my $subn = 'msa_check_internal2';
    $errcode->{long_internal_gap} = 0 if (!$errcode || !exists($errcode->{long_internal_gap}));
    $errcode->{has_internal_gap} = 0 if (!$errcode || !exists($errcode->{has_internal_gap}));

    $debug && print STDERR "$subn: \$loc2=$loc2 \$errcode=\n".Dumper($errcode)."\n";
    my ($qstart, $qend);
    my $gap_char = $aln->gap_char;
    $qstart = $reffeat_loc->start; # get start of reffeat
    $qstart = ($qstart - $refcds_loc->start)/3 +1; # convert to AA coordinate in CDS
    my $qseq = $aln_q->seq;
#    $debug && print STDERR "$subn: \$qseq='$qseq' \$qseq=".length($qseq)."\n";
    my $char;
    for (my $n=1; $n<=$qstart && $n<=length($qseq); $n++) {
            $char = substr($qseq, $n-1, 1);
            $qstart++ if ($char eq $gap_char);
#        $debug && print STDERR "$subn: \$n=$n \$char=$char \$gap_char=$gap_char \$qstart=$qstart\n";
    }

    $qend = $reffeat_loc->end; # get end of reffeat
    $qend = ($qend - $refcds_loc->start +1)/3; # convert to AA coordinate in CDS
    $qseq = $aln_q->seq;
#    $debug && print STDERR "$subn: \$qseq='$qseq' \$qseq=".length($qseq)."\n";
    for (my $n=1; $n<=$qend && $n<=length($qseq); $n++) {
        $char = substr($qseq, $n-1, 1);
        $qend++ if ($char eq $gap_char);
#        $debug && print STDERR "$subn: \$n=$n \$char=$char \$gap_char=$gap_char \$qend=$qend\n";
    }

    # take qstart & qend from $errcode
    if ( 1 ) {
        $qstart = $errcode->{alnstart};
        $qend   = $errcode->{alnend  };
    }

    # Count the number of gaps in the hit sequence
    my $loc2_str = $loc2->to_FTstring;
    $debug && print STDERR "$subn: \$loc2=$loc2_str \$qseq='$qseq' \$qseq=".length($qseq)."\n";
    $qseq = substr($qseq, $qstart-1, $qend-$qstart+1); # Any gap at the beginning of CDS, plus the 1st AA of current feat.
    $debug && print STDERR "$subn: \$loc2=$loc2_str \$qseq='$qseq' \$qseq=".length($qseq)."\n";
    $qseq =~ s/^[$gap_char]*//; # remove any gap at start
    $qseq =~ s/[$gap_char]*$//; # remove any gap at end

    my $hseq = $aln_h->seq;
    $debug && print STDERR "$subn: \$loc2=$loc2_str \$hseq='$hseq' \$hseq=".length($hseq)."\n";
    $hseq = substr($hseq, $qstart-1, $qend-$qstart+1); # Any gap at the beginning of CDS, plus the 1st AA of current feat.
    $debug && print STDERR "$subn: \$loc2=$loc2_str \$hseq='$hseq' \$hseq=".length($hseq)."\n";
    $hseq =~ s/^[$gap_char]*//; # remove any gap at start
    $hseq =~ s/[$gap_char]*$//; # remove any gap at end
    my $gapct = 0;
    my @chars = split(//, $qseq);
    # Count the gaps in querry
    foreach my $c (@chars) {
            $gapct++ if ($c eq $gap_char);
    }
    $debug && print STDERR "$subn: \$loc2=$loc2_str \$hseq='$hseq' \$hseq=".length($hseq)." \$gapct=$gapct\n";

    # Check if refseq has any long gap. If so, discard the annotation
    # The ultimate test is AF126284, in which the TF protein has only 6/45 different from non-frameshifted version.
    my $GAP_THRESH = (length($qseq)>30) ? 0.2 *length($qseq) : 5;
    if ($gapct > $GAP_THRESH) {
        if ( 1 ) {
            print STDERR "$subn: ERROR: refcds=".$refcds_loc->to_FTstring." refmat=".$reffeat_loc->to_FTstring." cds=".$cds_loc->to_FTstring." \$qstart=$qstart \$qend=$qend \$gap_char='$gap_char'\n";
            print STDERR "$subn: ERROR: \$qseq='$qseq' \$qseq=".length($qseq)." \$gapct=$gapct\n";
            $errcode->{long_internal_gap} = 1;
        }
    }
    $gapct = 0;
    @chars = split(//, $hseq);
    # Count the gaps in hit
    my $sequngapped = 0;
    my $gapped = 0;
    my $temp = 0;
    my $i = 0;
    while ($i<=$#chars) {
        my $c = $chars[$i];
        while ( $i<=$#chars && $c ne $gap_char ) {
            $temp++;
            $debug && print STDERR "$subn: \$i=$i \$c=$c \$temp=$temp \$sequngapped=$sequngapped \$gapped=$gapped \$gapct=$gapct\n";
            $c = $chars[++$i];
        }
        if ($temp<10) {
            $gapped += $temp;
        } else {
            $sequngapped += $temp;
        }
#        if ($temp>$sequngapped) {
#            $gapped += $sequngapped;
#            $sequngapped = $temp;
#        }
        $temp = 0;
        while ( $i<=$#chars && $c eq $gap_char) {
            $gapct++;
            $debug && print STDERR "$subn: \$i=$i \$c=$c \$temp=$temp \$sequngapped=$sequngapped \$gapped=$gapped \$gapct=$gapct\n";
            $c = $chars[++$i];
        }
        $debug && print STDERR "$subn: \$i=$i \$temp=$temp \$sequngapped=$sequngapped \$gapped=$gapped \$gapct=$gapct\n";
    }
    $debug && print STDERR "$subn: Summary: \$loc2=$loc2_str \$hseq='$hseq' \$hseq=".length($hseq)." \$gapct=$gapct\n";

  if ( 0 ) {
    if ($gapct > 0.333*length($hseq)) {
            print STDERR "$subn: ERROR: \$loc2=$loc2_str \$qstart=$qstart \$qend=$qend \$gap_char='$gap_char'\n";
            print STDERR "$subn: ERROR: \$loc2=$loc2_str \$hseq='$hseq' \$hseq=".length($hseq)." \$gapct=$gapct\n";
            $errcode->{long_internal_gap} = 1;
    }
  } else {
    # The criteria for a good mat_peptide is: --Nov 2012
    # 1. if the run of good alignment is longer than 10AA (eg AB182800), and
    # 2. total of gapped alignment is less than 33% of total alignment. Gapped sequence is defined as shorter than 10AA
    return 0 if (($gapped + $sequngapped)<1);
    my $pct = $gapped/($gapped + $sequngapped);
    $debug && printf STDERR "$subn: \$loc2=$loc2_str \$sequngapped=$sequngapped \$gapped=$gapped percent=%4.2f%%\n", $pct*100;
    # Added '$gapct<5 || ' to fish out those with only very short gap - Dec 20, 2012
    if ($gapct<5 || $sequngapped>=10 && $gapped/($gapped + $sequngapped) <0.333) {
#    if ($sequngapped>=10 && $gapped/($gapped + $sequngapped) <0.333) {
    } else {
        print STDERR "$subn: Summary: \$loc2=$loc2_str \$hseq='$hseq' \$hseq=".length($hseq)." \$gapct=$gapct\n";
        printf STDERR "$subn: \$loc2=$loc2_str \$sequngapped=$sequngapped \$gapped=$gapped percent=%4.2f%%\n", $pct*100;
        $errcode->{long_internal_gap} = 1;
    }
  }
    $errcode->{has_internal_gap} = 1 if ($gapct>0);

    $debug && print STDERR "$subn: \$loc2=$loc2_str \$errcode=\n".Dumper($errcode)."\n";

    return ($errcode->{long_internal_gap});

} # sub msa_check_internal2


sub msa_check_loc_start {
    my ($loc2, $reffeat_loc, $refcds_loc, $cds_loc, $aln, $aln_q, $aln_h, $errcode) = @_;

    my $debug = 0 || $debug_all;
    my $subn = 'msa_check_loc_start';
    # error code: OutsideCDS; partial_mat_peptide; GapAtCleaveSite
    $errcode->{ partial_mat_peptide } = 0 if (!$errcode || !exists($errcode->{partial_mat_peptide}));

        my $qstart = $reffeat_loc->start; # get start of reffeat
        $qstart = ($qstart - $refcds_loc->start)/3 +1; # convert to AA coordinate in CDS
        my $qseq = $aln_q->seq;
        my $hseq = $aln_h->seq;
        $debug && print STDERR "$subn: \$qseq='$qseq' \$qseq=".length($qseq)."\n";
        my $gap_char = $aln->gap_char;
        my ($char, $char1);
        my ($s, $s1);
        my $n1 = 1;
        for (my $n=1; $n<=$qstart && $n<=length($qseq); $n++) {
            $char = substr($qseq, $n-1, 1);
            $s .= $char;
            $s .= ' ' if (10* floor(($n-1)/10) == 10* floor($n/10-1));
            $char1 = substr($hseq, $n-1, 1);
            $s1 .= $char1;
            $s1 .= ' ' if (10* floor(($n-1)/10) == 10* floor($n/10-1));
            $qstart++ if ($char eq $gap_char);
            my $v = 50* floor(($n-1)/50);
            my $v1 = 50* floor($n/50-1);
            if ($v ==$v1 || $n==$qstart) {
                $debug && printf STDERR "$subn: \$n=%3d..%3d \$s =$s \n", $n1, $n;
                $debug && printf STDERR "$subn: \$n=%3d..%3d \$s1=$s1 \$gap_char=$gap_char \$qstart=$qstart\n", $n1, $n;
                $s = undef;
                $s1 = undef;
                $n1 = $n+1;
            }
        }

        # Any gap in hseq at the cleavage site makes the mat_peptide partial
        #my $hseq = $aln_h->seq;
        $debug && print STDERR "$subn: \$hseq='$hseq' \$hseq=".length($hseq)."\n";
        $char = substr($hseq, 0, $qstart); # Any gap at the beginning of CDS, plus the 1st AA of current feat.
        $debug && print STDERR "$subn: \$char='$char' \$char=".length($char)."\n";
        if ($char =~ /[$gap_char]+$/) {
            $debug && print STDERR "$subn: \$qstart=$qstart \$char=$char \$gap_char=$gap_char\n";
            $errcode->{partial_mat_peptide} = 1;
        }
    $debug && print STDERR "$subn: \$errcode=\n".Dumper($errcode)."\n";

    return ($errcode->{partial_mat_peptide});

} # sub msa_check_loc_start

sub msa_check_loc_end {
    my ($loc2, $reffeat_loc, $refcds_loc, $cds_loc, $aln, $aln_q, $aln_h, $errcode) = @_;

    my $debug = 0 || $debug_all;
    my $subn = 'msa_check_loc_end';
    # error code: OutsideCDS; partial_mat_peptide; GapAtCleaveSite
    $errcode->{ partial_mat_peptide } = 0 if (!$errcode || !exists($errcode->{partial_mat_peptide}));

        my $qend = $reffeat_loc->end; # get end of reffeat
        $qend = ($qend - $refcds_loc->start +1)/3; # convert to AA coordinate in CDS
        my $qseq = $aln_q->seq;
        my $hseq = $aln_h->seq;
        $debug && print STDERR "$subn: \$qseq='$qseq' \$qseq=".length($qseq)."\n";
        my $gap_char = $aln->gap_char;
        my ($char, $char1);
        my ($s, $s1);
        my $n1 = 1;
        for (my $n=1; $n<=$qend && $n<=length($qseq); $n++) {
            $char = substr($qseq, $n-1, 1);
            $s .= $char;
            $s .= ' ' if (10* floor(($n-1)/10) == 10* floor($n/10-1));
            $char1 = substr($hseq, $n-1, 1);
            $s1 .= $char1;
            $s1 .= ' ' if (10* floor(($n-1)/10) == 10* floor($n/10-1));
            $qend++ if ($char eq $gap_char);
            #$debug && print STDERR "$subn: \$n=$n \$char=$char \$gap_char=$gap_char \$qend=$qend\n";
            my $v = 50* floor(($n-1)/50);
            my $v1 = 50* floor($n/50-1);
            if ($v ==$v1 || $n==$qend) {
                $debug && printf STDERR "$subn: \$n=%3d..%3d \$s =$s \n", $n1, $n;
                $debug && printf STDERR "$subn: \$n=%3d..%3d \$s1=$s1 \$gap_char=$gap_char \$qend=$qend\n", $n1, $n;
                $s = undef;
                $s1 = undef;
                $n1 = $n+1;
            }
        }
        # Any gap in hseq at the cleavage site makes the mat_peptide partial
        #my $hseq = $aln_h->seq;
        $debug && print STDERR "$subn: \$hseq='$hseq' \$hseq=".length($hseq)."\n";
        $char = substr($hseq, $qend-1, length($hseq)-$qend+1); # Any gap at the end of CDS, plus the last AA of current feat.
        $debug && print STDERR "$subn: \$char='$char' \$char=".length($char)."\n";
        if ($char =~ /^[$gap_char]+/) {
            $debug && print STDERR "$subn: \$qend=$qend \$char=$char \$gap_char=$gap_char\n";
            $errcode->{partial_mat_peptide} = 1;
        }

    $debug && print STDERR "$subn: \$errcode=\n".Dumper($errcode)."\n";

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

    my $debug = 0 && $debug_all;
    my $subn = 'msa_check_loc_fuzzy';
    my $errcode = {};
    my $msg;
    my $max_gap_end = 3;
    my $max_gap_middle = 3;

    if ($loc2->start < $cds_loc->start
       ) {  # if $feat starts before CDS starts
            print STDERR "$subn: 1 ".$loc2->start.' '.$cds_loc->start."\n";
    }
    if ($loc2->end <= $cds_loc->start
       ) {  # if $feat end before CDS starts
            print STDERR "$subn: 2 ".$loc2->end.' '.$cds_loc->start."\n";
    }
    if ($loc2->start >= $cds_loc->end
       ) {  # if $feat starts after CDS ends
            print STDERR "$subn: 3 ".$loc2->start.' '.$cds_loc->end."\n";
    }
    if ($loc2->end > $cds_loc->end
       ) {  # if $feat ends after CDS ends
            print STDERR "$subn: 4 ".$loc2->end.' '.$cds_loc->end."\n";
    }
    if ($loc2->start < $cds_loc->start
         || $loc2->end <= $cds_loc->start
         || $loc2->start >= $cds_loc->end
         || $loc2->end > $cds_loc->end
       ) {  # if $feat is outside of CDS
            print STDERR "$subn: \$loc2=".$loc2->to_FTstring." doesn't match \$cds=".$cds_loc->to_FTstring."\n";
    # error code: 0: none; 1: gaps around cleavage site; 2: partial mat_peptide; 3: outside of CDS
            $errcode->{OutsideCDS} = 1;
            return ($errcode);
    }

    $debug && print STDERR "$subn: \n";


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
        foreach my $i (1 .. $#{$gaps_q}) {
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
    my $subn = 'convert_range_protein2dna';
    my @gaps2;
    my $f_start = $f->location->start;
    my $under_water = 0;
    my $above_water = 0;
    $debug && print STDERR "$subn: gaps are converted to DNA coordinates from protein coordinate.\n";
    $debug && print STDERR "$subn: DNA coordinates include the start of CDS, also include the gaps.\n";
    for (my $i=0; $i<=$#{$gaps}; $i++) {
        my $gap = $gaps->[$i];

        if ($gap eq '') {
                push @gaps2, $gap;
                next;
        }

        # This indicates the gap at the beginning of the CDS, special case
        if ($gap == $i) {
            my $sss = $f_start - ($gap*3)+0;
            $debug && print STDERR "$subn: \$i=$i \$gap=$gap start=$f_start new=". ($f_start-($gap*3)+2) ." \$under_water=$under_water \$sss=$sss\n";
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
            $debug && print STDERR "$subn: \$i=$i \$gap=$gap start=$f_start new=". ($f_start + ($gap-1)*3) ." \$under_water=$under_water \$sss=$sss\n";
            push @gaps2, $f_start -$under_water -$above_water + ($gap-1)*3;
            push @gaps2, $f_start -$under_water -$above_water + ($gap-1)*3+1;
            push @gaps2, $f_start -$under_water -$above_water + ($gap-1)*3+2;
        }
    }

    # check the gaps to ensure all gaps are included
    my $f_end = $f->location->end;
    for (my $i=1; $i<=$#gaps2; $i++) {
        $debug && print STDERR "$subn: \$i=$i \$gaps[$i]=$gaps2[$i] \$f_end=$f_end\n";
        if ($gaps2[$i] < $f_end) {
            $f_end++;
        } elsif ($gaps2[$i] == $f_end) {
            last;
        } elsif ($gaps2[$i] > $f_end) {
            @gaps2 = (@gaps2, $f_end+1 .. $gaps2[$i]-1);
            last;
        }
    }

    $debug && print STDERR "$subn: \$under_water=$under_water \$above_water=$above_water\n";
    my $gap1 = shift @gaps2;
    $debug && print STDERR "$subn: sorting \@gaps2=@gaps2\n";
    @gaps2 = sort {$a<=>$b} @gaps2;
    @gaps2 = ($gap1, @gaps2);
    $debug && print STDERR "$subn: sorted \@gaps2=@gaps2\n";
    return (\@gaps2);
} # sub convert_range_protein2dna


## sub get_dna_byloc fetches the sequence directly from seq string according to
## the location range. Sequence is always extracted from $seq_obj_seq
## Codon_start is taken care of, meaning the result excludes overhang
sub get_dna_byloc {
    my ($feat_obj, $seq) = @_;

    my $debug = 0 && $debug_all;
    my $subn = 'Annotate_Math::get_dna_byloc';
      
    $debug && print STDERR "$subn: \$feat_obj=\n".Dumper($feat_obj)."End of \$feat_obj\n\n";
    $debug && print STDERR "$subn: original \$seq=\n".Dumper($seq)."End of \$seq\n\n";
    my ($s, $product);
    if ($feat_obj->location->isa('Bio::Location::Simple')) {
        $debug && print STDERR "$subn: \$loc=\n".Dumper($feat_obj->location)."End of \$loc\n\n";
        $debug && print STDERR "$subn: \$loc=\n".Dumper($feat_obj)."End of \$loc\n\n";
#        print "  Simple:DNA    :seq=". $feat_obj->seq->seq ."\n";
#        $s = substr($seq, $feat_obj->location->start-1, $feat_obj->location->end - $feat_obj->location->start + 1);
        my $start = $feat_obj->location->start-1;
        my $len = $feat_obj->location->end - $feat_obj->location->start + 1;
        $debug && print STDERR "$subn: \$start = $start \$len=$len\n";
        $s = substr($seq, $start, $len);
        $debug && print STDERR "$subn: \$s=\n".Dumper($s)."End of \$s\n\n";
#        print "  Simple:protein:seq=". $feat_obj->seq->translate()->seq ."\n";
    } elsif ( $feat_obj->location->isa('Bio::Location::Split') ) {
        my $s0 = $seq;
#        print "location=". $feat_obj->location->to_FTstring . "\n";
        for my $loc ($feat_obj->location->sub_Location) {
            $debug && print STDERR "$subn: \$loc=\n".Dumper($loc)."End of \$loc\n\n";
            $debug && print STDERR "$subn: \$loc->strand()=".$loc->strand()."\n";
            $debug && print STDERR "$subn: \$loc->strand()=".$feat_obj->location->strand()."\n";
#            $s .= substr($s0, $loc->start -1, $loc->end +1 - $loc->start);
            my $start = $loc->start -1;
            my $len = $loc->end +1 - $loc->start;
            my $stemp = substr($s0, $start, $len);
              $s .= $stemp;
#           print "start=". $loc->start ."\tend=". $loc->end ."\n";
#           print "  \$s=". $s ."\n";
        }
    } elsif ($feat_obj->location->isa('Bio::Location::Fuzzy')) {
        my $s0 = $seq;
        $debug && print STDERR "$subn: location=". $feat_obj->location->to_FTstring . "\n";
        for my $loc ($feat_obj->location->each_Location) {
#            $s .= substr($s0, $loc->start -1, $loc->end +1 - $loc->start);
            my $start = $loc->start -1;
            my $len = $loc->end +1 - $loc->start;
            $s .= substr($s0, $start, $len);
            $debug && print STDERR "get_dna_byloc: start=". $loc->start ."\tend=". $loc->end ."\n";
            $debug && print STDERR "get_dna_byloc:   \$s=". $s ."\n";
        }
    }
    $debug && print STDERR "$subn: len=".length($s)." \$s = $s\n";

    # Take care of the positive/negative strand
    if ($feat_obj->strand() && $feat_obj->strand() == -1) {
        $s = revcom($s); # reverse the sequence, however this creates a 'Bio::PrimarySeq'
        $debug && print STDERR "$subn: after revcom \$s=\n".Dumper($s)."End of \$s\n\n";
        $s = $s->seq;
        $debug && print STDERR "$subn: after revcom \$s=\n".Dumper($s)."End of \$s\n\n";
    }

    # Take care of codon_start
    if ($feat_obj->has_tag('codon_start')) {
        my $codon_start = [$feat_obj->get_tag_values('codon_start')];
        $codon_start = $codon_start->[0];
        $s = substr($s, $codon_start-1);
    }
    $debug && print STDERR "$subn: len=".length($s)." \$s = $s\n";

    return $s;
} # sub get_dna_byloc

1;
