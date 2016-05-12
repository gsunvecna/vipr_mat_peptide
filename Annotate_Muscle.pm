package Annotate_Muscle;

use strict;
use warnings;
use English;
use Carp;
use Data::Dumper;

use version; our $VERSION = qv('1.1.0');
use Bio::SeqIO;
use Bio::Seq;
use Bio::AlignIO;
use Bio::Tools::Run::Alignment::Muscle;
use IO::String;
use Annotate_Util;
use Annotate_Verify;

my $debug_all = 1;

####//README//####
#
# Annotate_Muscle contains the core functions to perform annotation based on MUSCLE alignment.
# The alignment is done with AA sequences of ref CDS and target CDS
#
#    Authors Chris Larsen, clarsen@vecna.com; Guangyu Sun, gsun@vecna.com
#    May 2010
#
##################


=head2 annotate_1gbk

Takes a genbank file, finds hash of refseqs, and an array of target genomes, in the form of
 [ accession, [CDS, mat_peptides] ... ], []... ]
 run muscle on each CDS,
 Return an array of alignments CDS1, CDS2, CDS3 ...

=cut

sub annotate_1gbk {
    my ($gbk, $exe_dir, $aln_fn) = @_;
    $aln_fn = [] if (!defined($aln_fn));

    my $debug = 0 && $debug_all;
    my $subname = 'annotate_1gbk';

    my $in_file2 = IO::String->new($gbk);
    my $in  = Bio::SeqIO->new( -fh => $in_file2, -format => 'genbank' );
    # Only take 1st sequence from each gbk (Note: gbk can hold multiple sequences, we ignore all after 1st)
    my $inseq = $in->next_seq();

    my $refseqs = {};
    my $inseqs  = [];
    my $refpolyprots = [];
    my $polyprots    = [];
    # determine the refseq, and get the CDS/mat_peptides in refseq
    $refpolyprots = Annotate_Util::get_refpolyprots( $refseqs, $inseq, $exe_dir);
    $debug && print STDERR "$subname: \$refpolyprots = $#{@$refpolyprots}\n";
#    $debug && print STDERR "$subname: \$refpolyprots = \n".Dumper($refpolyprots)."End of \$refpolyprots\n\n";
#    $debug && print STDERR "$subname: \$refpolyprots = \n".Dumper($refpolyprots->[0]->[0]->seq)."End of \$refpolyprots\n\n";

    # According to refseq, get the CDS/mat_peptides in inseq, use bl2seq to determine if the CDS matches
    $polyprots = Annotate_Util::get_polyprots( $inseq, $refpolyprots);
    $debug && print STDERR "$subname:    \$polyprots = $#{@$polyprots}\n";
    $debug && print STDERR "$subname: \$polyprots = \n".Dumper($polyprots)."End of \$polyprots\n\n";
    $debug && print STDERR "$subname: \$polyprots = \n".Dumper($polyprots->[0]->[0])."End of \$polyprots\n\n";

    # add refseq to hash, add inseq to array
    if ($#{@$refpolyprots} >=0) {
        if (!exists($refseqs->{$refpolyprots->[0]->[0]->seq->accession_number})) {
            $refseqs->{$refpolyprots->[0]->[0]->seq->accession_number} = $refpolyprots;
        }
    } else {
        return undef;
    }
    if ($#{@$polyprots} <0) {
        print STDERR "$subname: no polyprotein was found in ".$inseq->accession_number."\n";
        print STDERR "$subname: \$polyprots = $#{@$polyprots}\n";
        return undef;
    }
    push @$inseqs, [$polyprots->[0]->[0]->seq->accession_number, $polyprots];

    my $feats_all;
    $feats_all = Annotate_Muscle::muscle_profile( $refseqs, $inseqs, $aln_fn);

    return $feats_all->[0];

} # sub annotate_1gbk


=head2 muscle_profile

Takes a hash of refseqs, and an array of target genomes, in the form of
 [ accession, [CDS, mat_peptides] ... ], []... ]
 run muscle on each CDS,
 Return an array of alignments CDS1, CDS2, CDS3 ...

=cut

sub muscle_profile {
    my ($refseqs, $inseqs, $aln_fn) = @_;

    my $debug = 0 && $debug_all;
    my $subname = 'muscle_profile';

    my $feats_all;

    my $feats;
    my $alns = [];
    my $cds_all = [];
    my $vfile1 = IO::String->new('');
    my $seq_out = Bio::SeqIO->new(
                          '-fh'     => $vfile1,
                          '-format' => 'fasta'
                          );

    my $n_sets; # holds the number of [CDS, mat_peptides] ... ] in each genome. Since each CDS is unique, separate MUSCLE run has to be carried out
    foreach my $set (@$inseqs) {
        $n_sets = $set->[1];
        next if (!@$n_sets);
        $n_sets = $#{@$n_sets};
    }
    $debug && print STDERR "$subname: \$n_sets=$n_sets\n";
    if ($aln_fn->[0] && $#{@$aln_fn}!=$n_sets) {
        $debug && print STDERR "$subname: \$n_sets=$n_sets not equal to # of aln files $#{@$aln_fn}\n";
        $debug && print STDERR "$subname: please check aln files, and re-try\n";
        return undef;
    }

    # for each unique CDS
    for my $n_set (0 .. $n_sets) {
        my $aln;
        my @param;
        my $factory;
        @param = (
                   '-stable' => '',
                   '-outfile_name' => "test.afa",
                 );
        $factory = Bio::Tools::Run::Alignment::Muscle->new(@param);

        # First, get the CDS from refseq
        foreach my $key (keys %$refseqs) {
            my $refcds = $refseqs->{$key}->[$n_set]->[0];
            $debug && print STDERR "$subname: \$refcds=\n". Dumper($refcds) . "End of \$refcds\n";

            my @values = $refcds->get_tag_values('translation');
            my $s1 = $values[0];
            $s1 = $1 if ($s1 =~ /([^*]+)[*]$/); # to remove any trailing * in CDS
#            my $str1 = $refcds->seq->accession_number .'|';
#            $str1 .= $refcds->primary_tag.'='.$refcds->location->to_FTstring;
            my $str1 = 'refseq='.$refcds->seq->accession_number .'|';
            $str1 .= $refcds->primary_tag.'='.$refcds->location->to_FTstring;
#            $debug && print "$subname: \$str=$str\n";
            my $f1 = Bio::PrimarySeq->new(
                         -seq      => $s1,
                         -id       => $str1,	# id can't contain space
                         -alphabet => 'protein'
                                    );

            push @$cds_all, $f1;
        }
#        $debug && print STDERR "$subname: \$cds_all=\n". Dumper($cds_all) . "End of \$cds_all\n";

        # Second, add CDS from target genomes
        for (my $i = 0; $i<=$#{@$inseqs}; $i++) {
            my $cds = $inseqs->[$i]->[$n_set+1]->[0]->[0];
            next if (!$cds);
#            $debug && print STDERR "$subname: \$cds=\n". Dumper($cds) . "End of \$cds\n";
#            my $s3 = $cds->seq->translate->seq;
            my $s3 = Annotate_Util::get_new_translation( $cds, $cds);
            $s3 = $1 if ($s3 =~ /([^*]+)[*]$/); # to remove any trailing * in CDS
            my @values = $cds->get_tag_values('translation');
            my $s2 = $values[0];
            $s2 = $1 if ($s2 =~ /([^*]+)[*]$/); # to remove any trailing * in CDS
            if ($s2 ne $s3) {
                print STDERR "$subname: \$s2='$s2'\n";
                print STDERR "$subname: \$s3='$s3'\n";
                print STDERR "$subname: translation tag and translate don't match in CDS of ".$cds->seq->accession_number.".\n";
                print STDERR "$subname: Skip.\n\n";
                next;
            }
#            my $str2 = $cds->seq->accession_number .'|';
#            $str2 .= $cds->primary_tag.'='.$cds->location->to_FTstring;
            my $str2 = 'ACC='.$cds->seq->accession_number .'|';
            $str2 .= $cds->primary_tag.'='.$cds->location->to_FTstring;
            my $f2 = Bio::PrimarySeq->new(
                         -seq      => $s2,
                         -id       => $str2,	# id can't contain space
                         -alphabet => 'protein'
                                    );
            push @$cds_all, $f2;
        }
#        $debug && print STDERR "$subname: \$cds_all=\n". Dumper($cds_all) . "End of \$cds_all\n";

        # Run MUSCLE for each CDS in refseq
        # Returns a SimpleAlign object
        $aln = $factory->align(
                         [@$cds_all]
                         );
#        $debug && print STDERR "$subname: \$aln=\n". Dumper($aln) . "End of \$aln\n\n";

        push @$alns, $aln;

        $feats = &MSA_annotate($refseqs, $inseqs, $n_set, $aln);

        push @$feats_all, $feats;

    } # for my $n_set (0 .. $n_sets)

  return $feats_all;

} # sub muscle_profile


=head2 MSA_annotate

Takes a hash of refseqs, in the form of { accession => [ [CDS, mat_peptides] ... ], []... }
 and an array of target genomes, in the form of [ accession, [ [CDS, mat_peptides, ...], []... ] ]
 run muscle on each CDS,
 Return an array of alignments CDS1, CDS2, CDS3 ...

=cut

sub MSA_annotate {
    my ($refseqs, $inseqs, $n_set, $aln) = @_;

    my $debug = 0 && $debug_all;
    my $subname = 'MSA_annotate';

    my $refset = [];    # [CDS, mat_peptides] ... ]

    # get $refset
    foreach my $key (sort keys %$refseqs) {
        $debug && print STDERR "$subname: \$key=$key\n";
        $refset = $refseqs->{$key}->[$n_set];
#        $debug && print STDERR "$subname: \$refset=\n". Dumper($refset) . "End of \$refset\n";
        my $refcds = $refseqs->{$key}->[$n_set]->[0];
        $debug && print STDERR "$subname: \$refcds=". $refcds->seq->accession_number ."\n";
#        $debug && print STDERR "$subname: \$refcds=\n". Dumper($refcds) ."End of \$refcds\n";
        last;
    }

    # get $inset
    my $feats_all;
    my $inset = [];    # [CDS, mat_peptides] ... ]
    for (my $i = 0; $i<=$#{@$inseqs}; $i++) {
        my $cds = $inseqs->[$i]->[$n_set+1]->[0]->[0]; # [$n_set+1] b/c the [accession, ...] in $inseqs
        next if (!$cds);
        $debug && print STDERR "$subname: \$cds=\n". Dumper($cds) . "End of \$cds\n";
        $debug && print STDERR "$subname: \$i=$i \$refset=\n". Dumper($refset) . "End of \$refset\n";

        $inset  = $inseqs->[$i]->[$n_set+1]->[0];
        $debug && print STDERR "$subname: \$i=$i \$inset=\n". Dumper($inset) . "End of \$inset\n";

        # find the new annotations
        my $feats_new = &MSA_annotate_1cds($refset, $inset, $aln);

        push @$feats_all, $feats_new;

#        $debug && print STDERR "$subname: \$feats_new=\n". Dumper($feats_new) . "End of \$feats_new\n";

    }

    return $feats_all;
} # sub MSA_annotate


=head2 MSA_annotate_1cds

Takes a hash of refseqs, in the form of { accession => [ [CDS, mat_peptides] ... ], []... }
 and an array of target genomes, in the form of [ accession, [ [CDS, mat_peptides, ...], []... ] ]
 run muscle on each CDS,
 Return an array of alignments CDS1, CDS2, CDS3 ...

=cut

sub MSA_annotate_1cds {
    my ($refset, $inset, $aln) = @_;

    my $debug = 0 && $debug_all;
    my $subname = 'MSA_annotate_1cds';

    my $refcds = $refset->[0];
    $debug && print STDERR "$subname: \$refcds='".$refcds->seq->accession_number."'\n";
#    $debug && print STDERR "MSA_annotate_1cds: \$refcds=\n". Dumper($refcds) . "End of \$refcds\n";

    my @values = $refcds->get_tag_values('translation');
    my $s1 = $values[0];
    $s1 = $1 if ($s1 =~ /([^*]+)[*]$/); # to remove any trailing * in CDS
    my $str1 = 'refseq='.$refcds->seq->accession_number .'|';
    $str1 .= $refcds->primary_tag.'='.$refcds->location->to_FTstring;
    my $refcds_id = $str1;
    my $aln_q = Annotate_Util::msa_get_aln( $aln, $refcds_id);


    my $cds = $inset->[0];
    $debug && print STDERR "$subname: \$cds='".$cds->seq->accession_number."'\n";
    next if (!$cds);
#    $debug && print STDERR "$subname: \$cds=\n". Dumper($cds) . "End of \$cds\n";

#    $debug && print STDERR "$subname: \$inset=\n". Dumper($inset) . "End of \$inset\n";

#    my $s3 = $cds->seq->translate->seq;
    my $s3 = Annotate_Util::get_new_translation( $cds, $cds);
    $s3 = $1 if ($s3 =~ /([^*]+)[*]$/); # to remove any trailing * in CDS
    @values = $cds->get_tag_values('translation');
    my $s2 = $values[0];
    $s2 = $1 if ($s2 =~ /([^*]+)[*]$/); # to remove any trailing * in CDS
    if ($s2 ne $s3) {
        print STDERR "$subname: \$s2='$s2'\n";
        print STDERR "$subname: \$s3='$s3'\n";
        print STDERR "$subname: translation tag and translate don't match in CDS of ".$cds->seq->accession_number.".\n";
        print STDERR "$subname: Skip.\n\n";
        next;
    }
    my $str2 = 'ACC='.$cds->seq->accession_number .'|';
    $str2 .= $cds->primary_tag.'='.$cds->location->to_FTstring;
    my $cds_id = $str2;
#    $debug && print STDERR "$subname: \$cds_all=\n". Dumper($cds_all) . "End of \$cds_all\n";
    my $aln_h = Annotate_Util::msa_get_aln( $aln, $cds_id);
    $debug && print STDERR "$subname: \$aln_q=\n".Dumper($aln_q)."End of \$aln_q\n";
    $debug && print STDERR "$subname: \$aln_h=\n".Dumper($aln_h)."End of \$aln_h\n";

    my $feats_all;

        $feats_all = &project_matpept(
                         $refset,
                         $inset,
                         $aln,
                         $aln_q,
                         $aln_h,
                         );

    $debug && print STDERR "$subname: \$feats_all=\n". Dumper($feats_all) . "End of \$feats_all\n";

    return $feats_all;
} # sub MSA_annotate_1cds


=head2 project_matpept

Takes references to refset, $inset, gaps in query & hit, refaln
 returns the list of features of CDS, mat_peptide, sig_peptide.

=cut

sub project_matpept {
    my ($refset, $inset, $aln, $aln_q, $aln_h) = @_;

    my $debug = 0 && $debug_all;
    my $subname = 'project_matpept';

    my @feats_all = (); # Holds all new features for the target genome
#    $debug && print STDERR "$subname: \$refset = \n".Dumper($refset)."End of \$refset\n";
    # Go through all features in refseq, map the corresponding features in targe genome,
    # first, check CDS that is labeled as polyprotein, then all mat_pepride (and sig_peptide) after such CDS
    
    for (my $ct = 0; $ct <= $#{@$refset}; $ct++) {
        # First get CDS from refseq & target
        my $reffeat = $refset->[$ct];
        $debug && print STDERR "$subname: #$ct is ".$reffeat->primary_tag." \t$reffeat\n";
        next if ($reffeat->primary_tag ne "CDS");
        my @prod = $reffeat->get_tag_values('product');
#        $debug && print STDERR "project_matpept: ct=$ct \$reffeat is ".$reffeat->primary_tag." prod=$prod[0]\n";
        # Skip any CDS not labeled as "polyprotein". Potential problem as some polyproteins are not so labeled
#        next if ($prod[0] !~ /polyprotein/i);
        next if (!Annotate_Util::is_polyprotein( $reffeat, ['product', 'note']));
        $debug && print STDERR "$subname: \$reffeat = \n".Dumper($reffeat)."End of \$reffeat\n\n";

        my $refcds = $reffeat;
        my $cds = $inset->[0];
        my %allowed_feats = ('mat_peptide' => 1, 'sig_peptide' => 1);
        while (($refset->[$ct+1]) && ($allowed_feats{$refset->[$ct+1]->primary_tag})) {
            $reffeat = $refset->[++$ct];
            $debug && print STDERR "$subname: #$ct is ".$reffeat->primary_tag." \t$reffeat\n";

            $debug && print STDERR "$subname: \$aln_q=$aln_q\n";
            $debug && print STDERR "$subname: \$aln_h=$aln_h\n";

            my $note = "Annotated by VIPRBRC, MSA, refseq=".$refcds->seq->accession_number();
            my $feat = Annotate_Util::assemble_new_feature(
                                 $refcds,
                                 $reffeat,
                                 $cds,
                                 $aln,
                                 $aln_q,
                                 $aln_h,
                                 $note,
                                 );
            if (!defined($feat)) {
                print STDERR "$subname: \$feat is undef, skip\n";
                next;
            }
            $debug && print STDERR "$subname: \$feat = \n".Dumper($feat)."End of \$feat\n\n";

            # See if this is entirely new annotation wrt genbank
            my $new = 1;
            $new = Annotate_Util::is_new_annotation($feat, $inset);
            if ($new) {
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
            }

            push @feats_all, $feat;
            $debug && print STDERR "$subname: \$feat = \n".Dumper($feat)."End of \$feat\n\n";
        }

        # check if there is any gaps in new annotation, any gaps is shown as >=?=<
        my ($has_gap, $str) = Annotate_Verify::check_ranges( $cds, \@feats_all);
        $debug && print STDERR "$subname: \$has_gap = $has_gap\n";

        if ($has_gap) {
            my $note = "Annotated by VIPRBRC, MSA, refseq=".$refcds->seq->accession_number();
            my $feats = Annotate_Util::fix_cleavage_gaps(
                                 \@feats_all,
                                 $cds,
                                 $refcds,
                                 $refset,
                                 $aln_q,
                                 $aln_h,
                                 $note,
                                 );
#            @feats_all = @$feats;
            ($has_gap, $str) = Annotate_Verify::check_ranges( $cds, \@feats_all);
            $debug && print STDERR "$subname: \$has_gap = $has_gap\n";
        }
        unshift @feats_all, $cds; # prepend CDS to the array

    } # while ($ct <= $#{@$reffeats})

    $debug && print STDERR "$subname: \@feats_all = \n".Dumper(@feats_all)."End of \@feats_all\n\n";

    return \@feats_all;
} # sub project_matpept


1;
