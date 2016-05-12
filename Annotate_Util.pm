package Annotate_Util;

use strict;
use warnings;
use English;
use Carp;
use Data::Dumper;

use version; our $VERSION = qv('1.1.3'); # February 09 2011
use File::Temp qw/ tempfile tempdir /;
use Bio::SeqIO;
use Bio::Seq;
use Bio::AlignIO;
use Bio::Tools::Run::StandAloneBlast;
use IO::String;

use Annotate_Math;
use Annotate_Verify;

my $debug_all = 1;

####//README//####
#
# Annotate_Util contains the core functions to perform annotation based on both MUSCLE and bl2seq alignment
#
#    Authors Chris Larsen, clarsen@vecna.com; Guangyu Sun, gsun@vecna.com
#    December 2010
#
##################

## //EXECUTE// ##

my $taxon = {
               taxon_loaded => 0,
               taxon_fn     => "Annotate_taxon_records.txt",
            };
my $gene_symbol2 = {
               symbol_loaded => 0,
               symbol_fn     => "Annotate_symbol_records.txt",
            };


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

    if (!$refseq || !$refseq->isa('Bio::Seq::RichSeq') || !check_refseq($refseq, $inseq)) {

        $refseq = Annotate_Util::get_refseq( $inseq, $exe_dir);
        if (!$refseq) {
            my $acc = $inseq->accession;
            my $id  = $inseq->species->ncbi_taxid;
            print STDERR "$subname: ERROR genome=".$acc." w/ taxid=".$id." is not covered in V$VERSION.";
            print STDERR " Please contact script author for any update.\n";
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

    # For a new refseq, go through all features in refseq, look for polyproptein CDS
    my $reffeats = [ $refseq->get_SeqFeatures ];

    for (my $refct = 0; $refct<=$#{@$reffeats}; $refct++) {
        my $reffeat = $reffeats->[$refct];
        $debug && print STDERR "$subname: \$refct=$refct \$reffeat=".$reffeat->primary_tag."\n";
        next if ($reffeat->primary_tag ne 'CDS'); # Looking for the first CDS

        my $acc = $reffeat->seq->accession_number;
        my $is_poly = 0;
        $is_poly = Annotate_Util::is_polyprotein( $reffeat, ['product', 'note'], $reffeats->[$refct+1]);
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
        my %allowed_tags = (
                             'mat_peptide' => 1,
                             'sig_peptide' => 1,
                           );
        while ($refct<$#{@$reffeats} && $allowed_tags{$reffeats->[$refct+1]->primary_tag}) {

            my $reffeat = $reffeats->[$refct+1];
            # Special cases to correct refseq location of preM of NC_001809
            if ($acc eq "NC_001809" && $reffeat->location->to_FTstring eq "466..972") {
                $reffeat->location->end(969); # The annotation in gbk has 1 extra codon at the end
                print STDERR "$subname: \$reffeat=\n".Dumper($reffeat)."End of $reffeat\n";
            }
            # Special cases to correct refseq location of preM of NC_003635
            if ($acc eq "NC_003635" && $reffeat->location->to_FTstring eq "<440..>700") {
                $reffeat->location->end(925); # The annotation stops at the pre part
                print STDERR "$subname: \$reffeat=\n".Dumper($reffeat)."End of $reffeat\n";
            }
            
            $debug && print STDERR "$subname: \$refct=$refct \$reffeat=".$reffeats->[$refct+1]->primary_tag.' '.$reffeats->[$refct+1]->location->to_FTstring."\n";
            if ($reffeats->[$refct+1]->primary_tag eq 'sig_peptide') {
                ++$refct;
                next;
            };
            
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


        # check if this is a new CDS, e.g. NC_001547.1 has a 2nd CDS that's part of 1st one
        # If so, choose the longer CDS, and combine the list of mat_peptides
        my $seen=0;
#        foreach my $j (keys %$polyprots) {
        foreach my $k (0 .. $#{@$refpolyprots}) {
            my $old_cds = $refpolyprots->[$k]->[1];
            $debug && print STDERR "$subname: \$k=$k \$old_cds=\n".Dumper($old_cds)."end of \$old_cds\n\n";
            $debug && print STDERR "$subname: \$k=$k \$refcds=\n".Dumper($refcds)."end of \$refcds\n\n";
            my $s1 = Annotate_Util::get_new_translation( $old_cds, $old_cds);
            my $s2 = Annotate_Util::get_new_translation( $refcds, $refcds);
            $debug && print STDERR "$subname: \$s1=$s1\n";
            $debug && print STDERR "$subname: \$s2=$s2\n";
            if ($s1 =~ /$s2/i || $s2 =~ /$s1/i) { # $s1 is related to $s2
              $seen=1;
              printf STDERR "$subname: old CDS=%d..%d is related to CDS=%d..%d, skip the shorter one\n", $old_cds->location->start, $old_cds->location->end, $refcds->location->start, $refcds->location->end;
              $old_cds = $refcds if ($s2 =~ /$s1/i);
              # See if each mat_peptide is a duplicate, leave out duplicates. E.g. NC_003899
              foreach my $refmatp (@$refmatps) {
                my $seen1 = 0;
                for my $n (2 .. $#{$refpolyprots->[$k]}) {
                  my $refpolyprot = $refpolyprots->[$k]->[$n]->location;
                  if ($refmatp->location->start == $refpolyprot->start || $refmatp->location->end == $refpolyprot->end) {
                     $seen1 = 1;
                  }
                }
                push @{$refpolyprots->[$k]}, $refmatp if (!$seen1);
              }
              last;
            } else {
              # if $s1 is not related to $s2
              $debug && print STDERR "$subname: \$s1 is not related to \$s2\n";
            }
        }
#        }

        if (!$seen && $#{@$refmatps}>=0) {
            $debug && print STDERR "$subname: $#{@$refmatps}+1 mat_peptides in ".$refseq->accession_number."|$cds_id\n";
            push @$refpolyprots, [
                               $cds_id,
                               $refcds,
                               @$refmatps,
                                 ];
        }
        $debug && print STDERR "$subname: \$#refpolyprots=$#{@$refpolyprots}+1\n";
        $debug && print STDERR "$subname: \$refpolyprots=\n".Dumper($refpolyprots)."end of \$refpolyprots\n\n";
    }

    return $refpolyprots;
} # sub get_refpolyprots


sub is_polyprotein {
    my ($feat, $allowed_tags, $feat2) = @_;
#        $is_poly = Annotate_Util::is_polyprotein( $reffeat, ['product', 'note']);

    my $debug = 0 && $debug_all;
    my $subname = 'is_polyprotein';

    my $is_poly = 0;
    foreach my $tag (@$allowed_tags) {
#        $debug && print STDERR "$subname: \$refct=$refct \$tag=$tag\n";
        next if (!$feat->has_tag($tag));
#        $debug && print STDERR "$subname: \$refct=$refct Found \$tag=$tag\n";
        my @tag = $feat->get_tag_values($tag);
        foreach my $t (@tag) {
#            $debug && print STDERR "$subname: \$refct=$refct \$t=$t\n";
            if ($t =~ /polyprotein/i) {
                $is_poly = 1;
                last;
            }
        }
        last if $is_poly;
    }

    if (!$is_poly && $feat2) { {
        $debug && print STDERR "$subname: \$is_poly=$is_poly \$feat2=$feat2\n";
        last if ($feat2->primary_tag ne 'mat_peptide');
        last if ($feat2->location->start < $feat->location->start);
        last if ($feat->location->end   < $feat2->location->end);
        $is_poly = 1;
    } }

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
    my $num_cds = -1;
    my $comment = 'No refpolyprots provided';
    if ($#{@$refpolyprots} < 0) {
            return ($polyprots, $num_cds, $comment);
    }
    $comment = 'No CDS defined in genome file';

    my $acc = $inseq->accession_number;
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
        my $found_polyprot = 0;
        for (my $ct = 0; $ct<=$#{@$feats}; $ct++) {
            my $feat = $feats->[$ct];
#            $debug && print STDERR "$subname: \$ct=$ct \$feat=".$feat->primary_tag."\n";
            if ($feat->primary_tag ne 'CDS') {
                # Skip those features such as source, gene, 5'URT
                $debug && print STDERR "$subname: \$ct=$ct \$feat=".$feat->primary_tag." is not CDS. Skip.\n";
                next;
            }
            $comment = 'Found CDS, but doesn\'t aligns with refseq';
            $debug && print STDERR "$subname: \$ct=$ct \$feat=".$feat->primary_tag."\n";

            my $cds = $feat;
            # polyprotein is determined by bl2seq similarity search, not by label of polyprotein
            $debug && print STDERR "$subname: Found polyprotein at \$ct=$ct\n";

            # double check to ensure the translation of CDS is right. For instance, DQ430819
            my $tr2 = &get_new_translation($cds, $cds);
            if (!$tr2) {
                print STDERR "$subname: get_new_translation returned translation='undef'. Skip.\n";
                $comment = 'Problem with translation of CDS';
                next; # Skip if no translation can be obtained
            }
            $debug && print STDERR "$subname: translation='$tr2'\n";

            my $match_refcds = 0;
            $match_refcds = &match_cds_bl2seq($cds, $refcds);
            if (!$match_refcds) {
                $debug && print STDERR "$subname: \$ct=$ct \$feat=".$feat->location->to_FTstring." doesn't match refcds".$refcds->location->to_FTstring."\n";
                $comment = "CDS sequence doesn't match refseq";
                next;
            }
            $found_polyprot = 1;

            # Found a polyprotein CDS, now collect all mat_peptide and sig_peptide following it
            my $matps = [];
            my %allowed_tags = ( 'mat_peptide' => 1, 'sig_peptide' => 1, 'misc_feature' => 1,);
            while ($ct<$#{@$feats} && $allowed_tags{$feats->[$ct+1]->primary_tag}) {
                ++$ct;
                next if ($feats->[$ct]->primary_tag eq 'misc_feature');
                my $feat = $feats->[$ct];
                $debug && print STDERR "$subname: \$ct=$ct \$feat=".$feat->primary_tag.' '.$feat->location->to_FTstring."\n";
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
                my $s1 = Annotate_Util::get_new_translation( $old_cds, $old_cds);
                my $s2 = Annotate_Util::get_new_translation( $cds, $cds);
                $debug && print STDERR "$subname: \$s1=$s1\n";
                $debug && print STDERR "$subname: \$s2=$s2\n";
                my $loc = $cds->location;
                my $old_loc = $old_cds->location;
                if ($s1 =~ /$s2/i) { # $s1 is same to or longer than $s2
                    printf STDERR "$subname: CDS=%d..%d is same to or longer than CDS=%d..%d, skip 2nd CDS\n", $old_loc->start, $old_loc->end, $loc->start, $loc->end;
                    $seen=1;
                    push @{$polyprots->{$j}->[$k]}, @$matps;
                    last;
                } elsif ($s2 =~ /$s1/i) {
                    # if $s1 is shorter than $s2, need to update $old_cds, and append any mat_peptide to list
                    printf STDERR "$subname: CDS=%d..%d is shorter than CDS=%d..%d, skip 1st CDS\n", $old_loc->start, $old_loc->end, $loc->start, $loc->end;
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
        if (!$found_polyprot) {
            print STDERR "$subname: No corresponding polyprotein for \$refseq=".$reffeat->seq->accession_number."|CDS=".$reffeat->location->to_FTstring." in \$inseq=$acc\n";
        }
        print STDERR "$subname: Found polyprotein for \$refseq=".$reffeat->seq->accession_number." \$i=$i in  \$inseq=$acc  \$num_cds=$num_cds\n";

    } #     for (my $i = 0; $i<=$#{@$refpolyprots}; $i++)

    my $n = [keys %$polyprots];
    $debug && print STDERR "$subname: polyprotein CDS for acc=$acc \$n=$#{@$n}\n";
    $debug && print STDERR "$subname: \$polyprots=\n".Dumper($polyprots)."end of \$polyprots\n\n";
    $comment = "Found $#{@$n} CDS for acc=$acc suitable for annotation";

    return ($polyprots, $num_cds, $comment);
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
    my $emsgs = '';

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
    return $match_cds if (!defined($hit_cds));
    while (my $hsp_cds = $hit_cds->next_hsp) {
        $debug && print STDERR "$subname: \$hsp_cds=\n".Dumper($hsp_cds)."end of \$hsp_cds\n\n";

        $s1 = $s1<$s2 ? $s1 : $s2;
        my $hsp_cds_conserved = (exists($hsp_cds->{'CONSERVED'})) ? $hsp_cds->{'CONSERVED'} : $hsp_cds->{'num_conserved'};
        $debug && print STDERR "$subname: \$s1=$s1 length of conserved=$hsp_cds_conserved\n";
        my $conserved_residues_required = 0.66667; # The target must have 66.667% residues conserved wrt refseq
        $conserved_residues_required = 0.5 if ($debug);
        if ($hsp_cds_conserved >= $s1 * $conserved_residues_required) {
            $match_cds = 1;
            $debug && print STDERR "$subname: CDS length=$s1 conserved=$hsp_cds_conserved.";
            $debug && print STDERR " Length of conserved meets required $conserved_residues_required\n";
            $debug && print STDERR "$subname: \$match_cds=$match_cds\n";
            last;
        } else {
            my $cs = $cds->location->start;
            my $rs = $refcds->location->start;
            $emsgs .= "$subname: CDS=$cs refcds=$rs length=$s1 conserved=$hsp_cds_conserved.";
            $emsgs .= " Length of conserved doesn't meet required $conserved_residues_required\n";
            $emsgs .= "$subname: \$match_cds=$match_cds\n";
            $emsgs .= "$subname: QUERY_SEQ='".$hsp_cds->query_string."'\n";
            $emsgs .= "$subname:           '".$hsp_cds->homology_string."'\n";
            $emsgs .= "$subname:   HIT_SEQ='".$hsp_cds->hit_string."'\n";
        }
    }
    $debug && (!$match_cds) && print STDERR "$emsgs";

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
            $refacc = Annotate_Util::get_refseq_acc( $inseq, $exe_dir);
            return undef if (!$refacc);
            if (-e "${exe_dir}refseq/${refacc}_matpeptide.gb") {
                $refseq_fn = $refacc.'_matpeptide.gb';
            } elsif (-e "${exe_dir}refseq/${refacc}_msaa.gb") {
                $refseq_fn = $refacc.'_msaa.gb';
            } elsif (-e "${exe_dir}refseq/${refacc}.gb") {
                $refseq_fn = $refacc.'.gb';
            } elsif (-e "${exe_dir}refseq/${refacc}.gbk") {
                $refseq_fn = $refacc.'.gbk';
            } elsif (-e "${exe_dir}refseq/${refacc}.genbank") {
                $refseq_fn = $refacc.'.genbank';
            } else {
                $refseq_fn = $refacc.'.gb';
            }
            print STDERR "get_refseq: REFSEQ is ${exe_dir}refseq/$refseq_fn\n";

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
            print STDERR "get_refseq: ERROR: Can't find refseq gb file: $refseq_fn\n";
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
    my ($inseq, $exe_dir) = @_;

    my $debug = 0 && $debug_all;
    my $subname = 'get_refseq_acc';

#    my $refseq_acc = undef;
#    return $refseq_acc if (!$inseq);

        # list of refseqs
        my $refseq_list = {
           # Flavivirus
           11079 => 'NC_000943', # Murray Valley encephalitis virus; Ver1.1.3
           11072 => 'NC_001437', # Japanese encephalitis virus
##           11073 => 'NC_001437', # Japanese encephalitis virus. Refseq has 13 mat_peptides, w/ fuzzy ranges
##           11076 => 'NC_001437', # Japanese encephalitis virus. Refseq has 13 mat_peptides, w/ fuzzy ranges
           11060 => 'NC_001474', # Dengue virus 2
           11069 => 'NC_001475', # Dengue virus 3
           11053 => 'NC_001477', # Dengue virus 1
           12637 => 'NC_001477', # Dengue virus
           11082 => 'NC_001563', # West Nile virus (lineage II strain 956)
           31658 => 'NC_001564', # Cell fusing agent virus; Ver1.1.3
           11084 => 'NC_001672', # Tick-borne encephalitis virus. Refseq has 13 mat_peptides, w/ fuzzy ranges
##           11092 => 'NC_001672', # Tick-borne encephalitis virus. Refseq has 13 mat_peptides, w/ fuzzy ranges
##           11094 => 'NC_001672', # Tick-borne encephalitis virus. Refseq has 13 mat_peptides, w/ fuzzy ranges
##           47300 => 'NC_001672', # Tick-borne encephalitis virus. Refseq has 13 mat_peptides, w/ fuzzy ranges
##          638787 => 'NC_001672', # Tick-borne encephalitis virus. Refseq has 13 mat_peptides, w/ fuzzy ranges
           11086 => 'NC_001809', # Louping ill virus; Ver1.1.3; Special treatment to correct end point of preM
           11089 => 'NC_002031', # Yellow fever virus (YFV). This refseq has 13 mat_peptides, w/ fuzzy ranges
           11070 => 'NC_002640', # Dengue virus 4
           64300 => 'NC_003635', # Modoc virus; Ver1.1.3
           64285 => 'NC_003675', # Rio Bravo virus; Ver1.1.3
           64280 => 'NC_003676', # Apoi virus; Ver1.1.3
           11083 => 'NC_003687', # Powassan virus; Ver1.1.3
##           11085 => 'NC_003690', # Langat virus
##          161675 => 'NC_003996', # Tamana bat virus
           11103 => 'NC_004102', # Hepatitis C virus genotype 1
           31646 => 'NC_004102', # Hepatitis C virus genotype 1a
           31647 => 'NC_004102', # Hepatitis C virus genotype 1b
           31649 => 'NC_004102', # Hepatitis C virus genotype 2a
           31650 => 'NC_004102', # Hepatitis C virus genotype 2b
           31655 => 'NC_004102', # Hepatitis C virus genotype 6a
          356426 => 'NC_004102', # Hepatitis C virus genotype 3a
           64312 => 'NC_004119', # Montana myotis leukoencephalitis virus; Ver1.1.3
          172148 => 'NC_004355', # Alkhurma hemorrhagic fever virus; Ver1.1.3
##           64294 => 'NC_005039', # Yokose virus
##           12542 => 'NC_005062', # Omsk hemorrhagic fever virus
##          218849 => 'NC_005064', # Kamiti River virus
           64286 => 'NC_006551', # Usutu virus; Ver1.1.3
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

           # Coronaviridae
#  11128 => 'NC_003045', #     11128 | Bovine coronavirus
#  11142 => 'NC_001846', #     11138 | Murine hepatitis virus
#  11144 => 'NC_006852', #     11138 | Murine hepatitis virus
#  28295 => 'NC_003436', #     28295 | Porcine epidemic diarrhea virus
# 290028 => 'NC_006577', #    290028 | Human coronavirus HKU1
# 694009 => 'NC_004718', #    694009 | Severe acute respiratory syndrome-related coronavirus
#  11152 => 'NC_010800', #    694014 | Avian coronavirus
#  11120 => 'NC_001451', #    694014 | Avian coronavirus

           # Caliciviridae
           11983 => 'NC_001959', # Norwalk virus
           95340 => 'NC_001959', # Norwalk virus
#          150080 => 'NC_001959', # Norwalk virus

           # Togaviridae
           # strain_id => 'accession_refseq',
           # Alphavirus
           11036 => 'NC_001449', # 11036 |     11036 | Venezuelan equine encephalitis virus; Ver1.1.3
           11027 => 'NC_001512', # 11027 |     11027 | O'nyong-nyong virus; Ver1.1.3
           11029 => 'NC_001544', # 11029 |     11029 | Ross River virus; Ver1.1.3
           11034 => 'NC_001547', # 11034 |     11034 | Sindbis virus; Ver1.1.3
           11020 => 'NC_001786', # 11020 |     11020 | Barmah Forest virus; Ver1.1.3
           11033 => 'NC_003215', # 11033 |     11033 | Semliki forest virus; Ver1.1.3
           59301 => 'NC_003417', # 59301 |     59301 | Mayaro virus; Ver1.1.3
           78540 => 'NC_003433', # 78540 |     84589 | Salmon pancreas disease virus; Ver1.1.3
           11021 => 'NC_003899', # 11021 |     11021 | Eastern equine encephalitis virus; Ver1.1.3
           44158 => 'NC_003900', # 44158 |     44158 | Aura virus; Ver1.1.3
           11039 => 'NC_003908', # 11039 |     11039 | Western equine encephalomyelitis virus; Ver1.1.3
           84589 => 'NC_003930', # 84589 |     84589 | Salmon pancreas disease virus; Ver1.1.3
#           37124 => 'NC_001512', # 37124 |     37124 | Chikungunya virus, refseq=O'nyong-nyong virus
           37124 => 'NC_004162', # 37124 |     37124 | Chikungunya virus; Ver1.1.3
#           59300 => 'NC_001544', # 59300 |     59300 | Getah virus, refseq=Ross River virus
           59300 => 'NC_006558', # 59300 |     59300 | Getah virus; Ver1.1.3
           11024 => 'NC_012561', # 11024 |     11024 | Highlands J virus; Ver1.1.3
           48544 => 'NC_013528', # 48544 |     48544 | Fort Morgan virus; Ver1.1.3
           # Rubivirus
           11041 => 'NC_001545', # 11041 |     11041 | Rubella virus; Ver1.1.3
           };
        $debug && print STDERR "$subname: \$refseq_list=\n".Dumper($refseq_list)."end of \$refseq_list\n\n";

#    my $refseq_acc = undef;
    my $refseq_acc = '';
    return $refseq_list if (!$inseq);

        $debug && print STDERR "$subname: Taxon=$taxon->{taxon_loaded}\n";
        if (!$taxon->{taxon_loaded} && -f "$exe_dir/$taxon->{taxon_fn}" ) {
            $debug && print STDERR "$subname: found taxon_fn=$exe_dir/$taxon->{taxon_fn}\n";

            open my $taxon_file, '<', "$exe_dir/$taxon->{taxon_fn}"
               or croak("$subname: '$exe_dir/$taxon->{taxon_fn}', but couldn't open: $OS_ERROR");
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
            close $taxon_file or croak "$subname: Couldn't close $exe_dir/$taxon->{taxon_fn}: $OS_ERROR";
            my @keys = keys %{$taxon->{speciesid}};
            $taxon->{taxon_loaded} = 1 if ($#keys>1);

            $debug && print STDERR "$subname: finished reading list file: '$exe_dir/$taxon->{taxon_fn}'.\n";
            $debug && print STDERR "$subname: \$taxon = \n".Dumper($taxon)."End of \$taxon\n\n";


        }
        $debug && print STDERR "$subname: \$taxon->{taxon_loaded}=$taxon->{taxon_loaded}\n";

        my $taxid = $inseq->species->ncbi_taxid;
        $debug && print STDERR "$subname: \$taxid=$taxid \$inseq=$inseq is a ".ref($inseq)."\n";


        # Determine the refseq by a 2-step process
        if (exists($refseq_list->{$taxid})) {
            $refseq_acc = $refseq_list->{$taxid};
        } elsif($taxon->{taxon_loaded}) {
            my @speciesids = keys %{$taxon->{'speciesid'}};
            foreach my $speciesid (@speciesids) {
                next if (!exists($taxon->{'speciesid'}->{$speciesid}->{$taxid}));
                # Could use 0 to signal the strains not covered
#                next if ($taxon->{'speciesid'}->{$speciesid}->{$taxid}!=1);
                my $refacc;
                $refacc = $taxon->{'speciesid'}->{$speciesid}->{$speciesid};
                $debug && print "$subname: \$taxid=$taxid \$speciesid=$speciesid \$refacc=$refacc\n";
                $refseq_acc = $refacc if ($refacc && $refacc =~ /^NC_0/i);
                print STDERR "$subname: \$taxid=$taxid \$speciesid=$speciesid \$refacc=$refacc \$refseq_acc=$refseq_acc\n";
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
    my ($reffeat,$exe_dir) = @_;

    my $debug = 0 && $debug_all;
    my $subname = 'get_gene_symbol';
#    print STDERR "get_gene_symbol: \$reffeat=\n".Dumper($reffeat)."end of \$reffeat\n\n";

    # These gene symbols are defined by CLarsen, after considering refseqs, e.g. NC_001477 & NC_009942
    my $gene_symbols;

#my $gene_symbol2 = {
#               symbol_loaded => 0,
#               symbol_fn     => "Annotate_symbol_records.txt",
#            };
    $debug && print STDERR "$subname: \$exe_dir=$exe_dir\n";
    $debug && print STDERR "$subname: \$gene_symbol2->{symbol_fn}=$gene_symbol2->{symbol_fn}\n";
    if (!$gene_symbol2->{symbol_loaded} && -f "$exe_dir/$gene_symbol2->{symbol_fn}" ) {
        $debug && print STDERR "$subname: \$gene_symbol2->{symbol_loaded}=$gene_symbol2->{symbol_loaded}\n";
        $debug && print STDERR "$subname: found symbol_fn=$exe_dir/$gene_symbol2->{symbol_fn}\n";

        open my $symbol_file, '<', "$exe_dir/$gene_symbol2->{symbol_fn}"
           or croak("$subname: found '$exe_dir/$gene_symbol2->{symbol_fn}', but couldn't open: $OS_ERROR");
        while (<$symbol_file>) {
                chomp;
                $debug && print STDERR "$subname: \$_=$_\n";
                next if (m/^[#]/x);  # Skip the comment
                next if (m/^\s*$/x); # Skip empty lines
                s/'//g; # Remove single quotes
                my $words = [split(/\s*;\s*/)];
                $debug && print STDERR "$subname: \$_='$_'\n";
                $debug && print STDERR "$subname: \$words($#{@$words})='@$words'\n";
                if ($#{@$words}<2) { # skip the lines without enough fields
                    $debug && print STDERR "$subname: not enough data in gene_symbol: '@$words'\n";
                    next;
                }

                my ($acc, $loc, $sym) = @{$words}[0..3];
                $loc = '0..0' if (!$loc);
                $debug && print STDERR "$subname: \$acc=$acc \$loc=$loc \$sym=$sym\n";
                if (exists($gene_symbol2->{$acc}->{$loc}) && $sym != $gene_symbol2->{$acc}->{$loc}) {
                  print STDERR "$subname: conflicting data, \$gene_symbol2->{$acc}->{$loc}=$gene_symbol2->{$acc}->{$loc}\n";
                  print STDERR "$subname: conflicting data, new \$sym=$sym\n";
                } else {
                    $gene_symbol2->{$acc}->{$loc} = $sym;
                }
        }
        close $symbol_file or croak "$subname: Couldn't close $exe_dir/$gene_symbol2->{symbol_fn}: $OS_ERROR";
        $gene_symbol2->{symbol_loaded} = 1;

        $debug && print STDERR "$subname: finished reading list file: '$exe_dir/$gene_symbol2->{symbol_fn}'.\n";
        $debug && print STDERR "$subname: \$gene_symbol2->{symbol_loaded}=$gene_symbol2->{symbol_loaded}\n";
        $debug && print STDERR "$subname: \$gene_symbol2 = \n".Dumper($gene_symbol2)."End of \$gene_symbol2\n\n";
    }

    $gene_symbols = $gene_symbol2;
    my $gene_symbol = 'unk';
    my $acc = $reffeat->seq->accession_number;
    my $loc = $reffeat->location->to_FTstring;
    $loc = $reffeat->location->start .'..'. $reffeat->location->end;
    $gene_symbol = $gene_symbols->{$acc}->{$loc} ? $gene_symbols->{$acc}->{$loc} : $gene_symbol;
    $debug && print "get_gene_symbol: primary_tag=".($reffeat->seq->accession_number);
    $debug && print "\tlocation=".($reffeat->location->to_FTstring);
    $debug && print "\tsymbol=$gene_symbol\n";

    return $gene_symbol;
} # sub get_gene_symbol


=head2 check_refseq

Takes 1 refseq and 1 target genome. Check if their taxids are same

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

    my $s = Annotate_Math::get_dna_byloc( $feat, $cds->{_gsf_seq}->seq);
    $debug && print STDERR "$subname: \$s  ='$s'\n";
    my $f = Bio::PrimarySeq->new(
                         -seq      => $s,
                         -id       => '',	# id can't contain space
                         -alphabet => 'dna'
                                );
    $f->revcom() if ($feat->strand() == -1);

    $s = $f->translate()->seq;
    $s =~ s/[*]$// if ($s =~ /[*]$/);
    $s =~ s/[*]/./g if ($s =~ /[*]/);
    $s =~ s/[X]/./ig if ($s =~ /[X]/i);
    $debug && print STDERR "$subname: \$s  ='$s'\n";
    $debug && print STDERR "$subname: \$s  =".length($s)."\n";

    # Only keep the annotated mat_peptide that confirms to CDS. It should, just to double check
    my @parent_cds_seq = $cds->get_tag_values('translation');
    $debug && print STDERR "$subname: \$cds  =$parent_cds_seq[0]\n";
    $debug && print STDERR "$subname: \$cds  =".length($parent_cds_seq[0])."\n";

    if ($parent_cds_seq[0] =~ /($s)/i) {
        $translation = $1;
        $debug && print STDERR "$subname: \$s  ='$s'\n";
        $debug && print STDERR "$subname: \$cds='".$parent_cds_seq[0]."'\n";
        $debug && print STDERR "$subname: diff='".Annotate_Verify::diff_2str( $s, $parent_cds_seq[0])."'\n";
    } else {
        print STDERR "$subname: ERROR: \$acc=$acc translation for feature=".$feat->location->to_FTstring." doesn't match CDS.\n";
        print STDERR "$subname: \$s  ='$s'\n";
        print STDERR "$subname: \$cds='".$parent_cds_seq[0]."'\n";
        print STDERR "$subname: diff='".Annotate_Verify::diff_2str( $s, $parent_cds_seq[0])."'\n";
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
    my ($refcds, $reffeat, $cds, $aln, $refcds_id, $cds_id, $note,$exe_dir) = @_;

    my $debug = 0 && $debug_all;
    my $subname = 'assemble_new_feature';

    my $aln_q = Annotate_Util::msa_get_aln( $aln, $refcds_id);
    my $aln_h = Annotate_Util::msa_get_aln( $aln, $cds_id);
    my ($loc2, $errcode) = Annotate_Math::msa_get_feature_loc(
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
        # codon_start is taken care of during calculation of start/end of new feature. No need here
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

    my $gene_symbol = Annotate_Util::get_gene_symbol( $reffeat,$exe_dir);
    my ($id, $desc) = Annotate_Util::get_feature_id_desc( $feat, $errcode, $gene_symbol, $cds, $refcds);
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
    my ($feat, $errcode, $gene_symbol, $cds, $refcds) = @_;

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

    # add the accession of target genome
    if ( 1 ) {
        my $desc = '';
        $desc .= 'src=MSA|'; # Add method to description
        $desc .= 'ACC='. $cds->seq->accession_number; # Add accession to description
        $desc .= '|Ver='. $cds->seq->accession_number; # Add accession to description
        $desc = $desc .'.'.$cds->entire_seq->{'_version'}; # Add version of accession to description
        $debug && print STDERR "$subname: \$cds=\n".Dumper($cds)."End of \$cds\n\n";
        $debug && print STDERR "$subname: \$seq=\n".Dumper($cds->entire_seq)."End of \$seq\n\n";
        $id = $desc .'|'. $id;

        # add accession.version of refseq
        $desc = 'ref='. $refcds->seq->accession_number; # Add accession to description
        $desc .= '.'.$refcds->entire_seq->{'_version'}; # Add version of accession to description
        $debug && print STDERR "$subname: \$refcds=\n".Dumper($refcds)."End of \$refcds\n\n";
        $debug && print STDERR "$subname: \$refseq=\n".Dumper($refcds->entire_seq)."End of \$refseq\n\n";
        $id .= $desc .'|';
    }

    # add the range of mat_peptide
    $id .= 'Loc='. &get_DNA_loc($feat) .'|';

    # Add range in polyprotein. This is simply calculated from existing DNA coordinates.
    # One caveat is, in some cases, the last AA of polyprotein may have 2 coden, need special handling.
    my $s = '';
    $s = &get_AA_loc($feat, $cds);
    $id .= 'AA='. $s .'|';

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
               $desc = $desc. $feat->primary_tag."=". $id[0];
           }
           $id .= 'product='.$id[0];
           last;
        }
    }
    $debug && print STDERR "$subname: \$desc='$desc'\n";

    return ($id, $desc);
} # sub get_feature_id_desc


=head2 get_DNA_loc

Takes $feat
 Returns the location range in a string
 This basically is $feat->location->to_FTstring, but we are turning off the <> notation
=cut

sub get_DNA_loc {
    my ($feat) = @_;

    my $debug = 0 && $debug_all;
    my $subname = 'get_DNA_loc';

    # get the range of mat_peptide
    my $location_allow_split = 1;
    my $location_allow_fuzzy = 0;
    my $s;
    if ($location_allow_split) {
        $s = $feat->location->to_FTstring;
        $debug && print STDERR "$subname: may contain '<>' in \$s=$s.\n";
    } else {
        $s = $feat->location->start .'..'. $feat->location->end;
        $debug && print STDERR "$subname: no '<>' in \$s=$s.\n";
    }

    $s =~ s/[<>]//g if (!$location_allow_fuzzy);

    return $s;
} # sub get_DNA_loc


=head2 get_AA_loc

Takes mat_peptide and the parent CDS
 Return the AA location

=cut

sub get_AA_loc {
    my ($feat, $cds) = @_;

    my $debug = 0 && $debug_all;
    my $subname = 'Annotate_Util::get_AA_loc';

    # Add range in polyprotein. This is simply calculated by comparing the sequence of mat_peptide with
    # that of the CDS, and finding the index.
    my $s = '0..0';
    if (1) {
        $debug && print STDERR "$subname: \$feat=".$feat->location->to_FTstring." \$cds=".$cds->location->to_FTstring."\n";
        my $s_temp = -1;
        my $s1 = Annotate_Util::get_new_translation( $feat, $cds);
        $s1 = '' if (!$s1);
        my $s0 = '';
        $s0 = [ $cds->get_tag_values('translation') ]->[0];
        $s_temp = $feat->location->start - $cds->location->start;
        $s_temp = $s_temp / 3 - 10;
        $s_temp = index($s0, $s1, $s_temp)+1;
        $debug && print STDERR "$subname: \$s_temp=$s_temp \$s1=$s1\n";
        $debug && print STDERR "$subname: \$s_temp=$s_temp \$s0=$s0\n";
        $s = sprintf("%d..%d", $s_temp, $s_temp+length($s1)-1) if ($s_temp>0);
        $debug && print STDERR "$subname: \$s=$s \$s1=$s1\n";
        $debug && print STDERR "$subname: \$s=$s \$s0=$s0\n";
    } else {
        my $codon_start_feat = 0;
        $codon_start_feat = [$feat->get_tag_values('codon_start')]->[0] -1 if ($feat->has_tag('codon_start'));
        my $codon_start_cds = 0;
        $codon_start_cds = [$cds->get_tag_values('codon_start')]->[0] -1 if ($cds->has_tag('codon_start'));

        my $b0 = $feat->location->start +$codon_start_feat - $cds->location->start -$codon_start_cds;
        $debug && print STDERR "$subname: \$b0=$b0 \$b0/3=". $b0 % 3 ."\n";
        $b0 = (($b0 % 3) ==0) ? $b0 /3 +1 : 0;
        $debug && print STDERR "$subname: \$b0=$b0\n";

        my $e0 = $feat->location->end +$codon_start_feat + 1 - $cds->location->start -$codon_start_cds;
        $debug && print STDERR "$subname: \$e0=$e0\n";
        $e0 +=1 if (($e0 +1) % 3 ==0); # some of the polyproteins have 2 nucleotides for the last codon
        $debug && print STDERR "$subname: \$e0=$e0\n";
        $e0 = ($e0 % 3 ==0) ? $e0 /3 : 0;
        $debug && print STDERR "$subname: \$e0=$e0\n";
        $debug && print STDERR "$subname: \$b0=$b0 \$e0=$e0\n";
        $s = $e0 - $b0 +1; # This is the length calculated from the DNA location

        # For split locations, any gap between sublocations needs to be accounted for
        if ($feat->location->isa('Bio::Location::Split')) {
            my $locs = [ $feat->location->sub_Location ];
            my $gap = 0;
            for my $i (1 .. $#{$locs}) {
               $gap = $locs->[$i]->start - $locs->[$i-1]->end -1; 
               $debug && print STDERR "$subname: \$gap=$gap\n";
            }
            $debug && print STDERR "$subname: \$s=$s \$gap=$gap\n";
            $s -= $gap/3 if (($gap%3)==0);
        }

        my $b = $b0;
        my $e = $e0;
        $debug && print STDERR "$subname: \$feat=\n".Dumper($feat)."End of \$feat\n\n";
        # check the length of existing translation against length from coordinates, just in case
        if ($feat->has_tag('translation')) {
            my $t = [ $feat->get_tag_values('translation') ];
            $t = $t->[0];
            $debug && print STDERR "$subname: \$seq=$t\n";
            if ($s!=length($t)) {
                print STDERR  "$subname: ".$cds->seq->accession_number ." has problem between existing translation and coordinate \$s=$s length=".length($t)."\n";
                print STDERR "$subname: \$s=$s \$t=".length($t)."\n";
                print STDERR "$subname: \$feat=\n".Dumper($feat)."End of \$feat\n\n";
                print STDERR "$subname: \$cds=\n".Dumper($cds)."End of \$cds\n\n";
                croak "$subname: ".$cds->seq->accession_number ." has problem between existing translation and coordinate \$s=$s length=".length($t)."\n";
            }
        }
        $s = $b .'..'. $e;
    }

    return $s;
} # sub get_AA_loc


=head2 project_matpept

Takes references to refset, $inset, alignment, ids of refcds and cds, and a note
 returns the list of features of CDS, mat_peptide, sig_peptide.

=cut

sub project_matpept {
    my ($refset, $inset, $aln, $refcds_id, $cds_id, $note,$exe_dir) = @_;

    my $debug = 0 && $debug_all;
    my $subname = 'project_matpept';

    my $feats_all = []; # Holds all new features for the target genome
#    $debug && print STDERR "$subname: \$refset = \n".Dumper($refset)."End of \$refset\n";
    # Go through all features in refseq, map the corresponding features in targe genome,
    # first, check CDS that is labeled as polyprotein, then all mat_pepride (and sig_peptide) after such CDS
    
    $debug && print STDERR "$subname: \$refset = \n".Dumper($refset)."End of \$refset\n\n";
    for (my $ct = 1; $ct <= $#{@$refset}; $ct++) {
        # First get CDS from refseq & target
        my $reffeat = $refset->[$ct];
        $debug && print STDERR "$subname: #$ct is ".$reffeat->primary_tag." \t$reffeat\n";
        if ($reffeat->primary_tag ne "CDS") {
            $debug && print STDERR "$subname: ERROR: refset($ct) is not CDS\n";
            next;
        }
        my @prod = $reffeat->get_tag_values('product');
#        $debug && print STDERR "project_matpept: ct=$ct \$reffeat is ".$reffeat->primary_tag." prod=$prod[0]\n";
        # Skip any CDS not labeled as "polyprotein". Potential problem as some polyproteins are not so labeled
#        next if ($prod[0] !~ /polyprotein/i);
        next if (!Annotate_Util::is_polyprotein( $reffeat, ['product', 'note'], $refset->[$ct+1]));
        $debug && print STDERR "$subname: \$reffeat = \n".Dumper($reffeat)."End of \$reffeat\n\n";

        my $refcds = $reffeat;
        my $cds = $inset->[0];
        $debug && print STDERR "$subname: \$cds=".$cds->seq->accession_number."\n";
        my %allowed_feats = ('mat_peptide' => 1, 'sig_peptide' => 1);
        while (($refset->[$ct+1]) && ($allowed_feats{$refset->[$ct+1]->primary_tag})) {
            $reffeat = $refset->[++$ct];
            $debug && print STDERR "$subname: \$reffeat #$ct is ".$reffeat->primary_tag."=\n".Dumper($reffeat)."End of \$reffeat\n\n";

            $debug && print STDERR "$subname: \$refcds_id=$refcds_id \$cds_id=$cds_id\n";

            # gets a Bio::SeqFeature::Generic object
            my $feat = Annotate_Util::assemble_new_feature(
                                 $refcds,
                                 $reffeat,
                                 $cds,
                                 $aln,
                                 $refcds_id,
                                 $cds_id,
                                 $note,$exe_dir,
                                 );
            if (!defined($feat)) {
                $debug && print STDERR "$subname: #$ct \$feat is undef, skip\n";
                next;
            }
            $debug && print STDERR "$subname: #$ct \$feat = \n".Dumper($feat)."End of \$feat\n\n";

            # See if this is entirely new annotation wrt genbank
            my $new = 1;
            my $feat_gbk;
            ($new, $feat_gbk) = Annotate_Util::is_new_annotation( $feat, $inset);
            if ($new) {
            	# add *new* to the description to indicate that GBK doesn't have this annotation
                my @tags = $feat->get_tag_values('note');
                $feat->remove_tag('note');
                for (my $i = 0; $i<=$#tags; $i++) {
                    $debug && print STDERR "$subname: \$tags[$i] = $tags[$i]\n";
                    if ($tags[$i] =~ /^Desc:/i) {
                        $tags[$i] = $tags[$i] .'|*new*';
                        $feat->add_tag_value('note', $tags[$i]);
                        $debug && print STDERR "$subname: \$tags[$i] = $tags[$i]\n";
                    } else {
                        $feat->add_tag_value('note', $tags[$i]);
                    }
                }
            } elsif($feat_gbk) {
                # add source GBK and other identifiers (GI) to description
                my $id_gbk = '';
                for my $tag ('db_xref', 'protein_id', 'product') { # per Chris request
                    if ($feat_gbk->has_tag($tag)) {
                       $id_gbk = [ $feat_gbk->get_tag_values($tag) ];
                       $id_gbk = $id_gbk->[0];
                       last;
                    }
                }

            }

            push @$feats_all, $feat;
            $debug && print STDERR "$subname: \$feat = \n".Dumper($feat)."End of \$feat\n\n";
            $debug && print STDERR "$subname: \$feats_all = \n".Dumper($feats_all)."End of \$feats_all\n\n";
        }
        unshift @$feats_all, $cds; # prepend CDS to the array

        # check if there is any gaps in new annotation, any gaps is shown as >=?=<
        my ($has_gap, $str) = Annotate_Verify::check_ranges( $feats_all);
        $debug && print STDERR "$subname: \$has_gap = $has_gap\n";
        if ($has_gap) {
            my $feats = Annotate_Util::fix_cleavage_gaps( $feats_all, $refset, $inset,);
            ($has_gap, $str) = Annotate_Verify::check_ranges( $feats_all);
            $debug && print STDERR "$subname: \$has_gap = $has_gap\n";
        }

        # check if there is any partial mat_peptide in the middle of polyprotein
        ($has_gap, $str) = Annotate_Verify::check_partial( $feats_all);
        $debug && print STDERR "$subname: \$has_gap = $has_gap\n";
        if ($has_gap) {
#            my $feats = Annotate_Util::fix_partial( $feats_all, $refset, $inset,);
#            ($has_gap, $str) = Annotate_Verify::check_partial( $feats_all);
#            $debug && print STDERR "$subname: \$has_gap = $has_gap\n";
        }

    } # while ($ct <= $#{@$reffeats})

    $debug && print STDERR "$subname: \$feats_all = \n".Dumper($feats_all)."End of \$feats_all\n\n";

    return $feats_all;
} # sub project_matpept


=head2 fix_partial

Takes the new annotation in an array, refcds, cds, and the gaps.
 Check if all new mat_peptides are connected
 if not, decide if the missing sequence at the cleavage site should go to the tail of previous mat_peptide,
 or head of the next mat_peptide

If the alignment is not definite, as in EF407458 EF407463 EF407467 when run together and separately, the annotation will change according to the MSA. This is a consequence of the shifting alignment.

=cut

sub fix_partial {
#    my ($feats_all, $refset, $gaps_q, $gaps_h, $note) = @_;
    my ($feats_all, $refset, $inset) = @_;

    my $debug = 0 && $debug_all;
    my $subname = 'fix_partial';

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
                my $new = &is_new_annotation( $feat2, $inset);
                $pattn = '\|[*]new[*]';
                if (!$new && $value =~ s/($pattn)//i) {
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
} # sub fix_partial


=head2 fix_cleavage_gaps

Takes the new annotation in an array, refcds, cds, and the gaps.
 Check if all new mat_peptides are connected
 if not, decide if the missing sequence at the cleavage site should go to the tail of previous mat_peptide,
 or head of the next mat_peptide

If the alignment is not definite, as in EF407458 EF407463 EF407467 when run together and separately, the annotation will change according to the MSA. This is a consequence of the shifting alignment.

=cut

sub fix_cleavage_gaps {
#    my ($feats_all, $refset, $gaps_q, $gaps_h, $note) = @_;
    my ($feats_all, $refset, $inset) = @_;

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
#        print STDOUT $msg;

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
                my $new = &is_new_annotation( $feat2, $inset);
                $pattn = '\|[*]new[*]';
                if (!$new && $value =~ s/($pattn)//i) {
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

    $has_gap && $debug && print STDERR "$subname: \$feats_all = \n".Dumper($feats_all)."End of \$feats_all\n\n";

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
    my $subname = 'is_new_annotation';

    $debug && print STDERR "\n$subname: \$feat = \n".Dumper($feat)."\$feat\n";
    $debug && print STDERR "\n$subname: \$inset = \n".Dumper($inset)."\$inset\n";
    my $new = 1;
    my $feat_gbk;
    $debug && print STDERR "$subname: \$new = $new\n";
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
        $debug && print STDERR "$subname: \$str2='$str2' \$str1='$str1'\n";
        if ($str1 eq $str2) {
            $new = 0;
            $feat_gbk = $gbmatp;
            last;
        }
    }

    return ($new, $feat_gbk);
} # sub is_new_annotation


1;
