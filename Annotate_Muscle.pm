package Annotate_Muscle;

use strict;
use warnings;
use English;
use Carp;
use Data::Dumper;

use version; our $VERSION = qv('1.1.2'); # December 01 2010
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
#    December 2010
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

    my $comment = '';
    my $acc = $inseq->accession_number;
    my $refseqs = {};
    my $inseqs  = [];
    my $refpolyprots = [];
    my $polyprots    = [];
    # determine the refseq, and get the CDS/mat_peptides in refseq
    $refpolyprots = Annotate_Util::get_refpolyprots( $refseqs, $inseq, $exe_dir);
    $debug && print STDERR "$subname: \$refpolyprots = $#{@$refpolyprots}\n";
    $debug && print STDERR "$subname: \$refpolyprots = \n".Dumper($refpolyprots)."End of \$refpolyprots\n\n";
#    $debug && print STDERR "$subname: \$refpolyprots = \n".Dumper($refpolyprots->[0]->[0]->seq)."End of \$refpolyprots\n\n";

    # According to refseq, get the CDS/mat_peptides in inseq, use bl2seq to determine if the CDS matches
    my $num_cds;
    ($polyprots, $num_cds, $comment) = Annotate_Util::get_polyprots( $inseq, $refpolyprots);
    $debug && print STDERR "$subname:    \$num_cds = $num_cds\n";
    $debug && print STDERR "$subname: \$polyprots = \n".Dumper($polyprots)."End of \$polyprots\n\n";
    # Skip the refseqs that has mat_peptide annotation from genbank
    if ($acc =~ /^NC_\d+$/i) {
    	if ( 1 ) {
            $comment = 'Refseq with mat_peptide annotation from NCBI, skip';
            return (undef, $comment);
    	}
        my $has_mat_peptide = 0;
        CHECK: for my $key (keys %$polyprots) {
            my $sets = $polyprots->{$key};
            $debug && print STDERR "$subname: \$sets=\n".Dumper($sets)."End of \$sets\n\n";
            for my $j (0 .. $#{@$sets}) {
                my $set = $sets->[$j];
                $debug && print STDERR "$subname: \$set=\n".Dumper($set)."End of \$set\n\n";
                for my $k (1 .. $#{@$set}) {
                    next if ($set->[$k]->primary_tag ne 'mat_peptide');
                    $has_mat_peptide = 1;
                    last CHECK;
                }
            }
        }
    }

    # add refseq to hash, add inseq to array
    if ($#{@$refpolyprots} >=0) {
        if (!exists($refseqs->{$refpolyprots->[0]->[1]->seq->accession_number})) {
            $refseqs->{$refpolyprots->[0]->[1]->seq->accession_number} = $refpolyprots;
        }
    } else {
    	$comment = 'No suitable refseq was found for taxid='. $inseq->species->ncbi_taxid;
    	return (undef, $comment);
    }
    if ($num_cds<0) {
        print STDERR "$subname: no polyprotein (suitable CDS) found in ".$inseq->accession_number." \$polyprots=$num_cds\n";
    	$comment = "No polyprotein: ". $comment;
    	return (undef, $comment);
    }
    push @$inseqs, [$inseq->accession_number, $polyprots];

    my $feats_all;
    $feats_all = Annotate_Muscle::muscle_profile( $refseqs, $inseqs, $aln_fn,$exe_dir);

    my $fs = $feats_all->[0]->[0];
    my $n = $#{@$fs};
  	$comment = "MSA returned $n mat_peptides";
  	
    return ($feats_all->[0], $comment);

} # sub annotate_1gbk


=head2 muscle_profile

Takes a hash of refseqs, and an array of target genomes, in the form of
 [ accession, [CDS, mat_peptides] ... ], []... ]
 run muscle on each CDS,
 Return an array of alignments CDS1, CDS2, CDS3 ...

=cut

sub muscle_profile {
    my ($refseqs, $inseqs, $aln_fn,$exe_dir) = @_;

    my $debug = 0 && $debug_all;
    my $subname = 'muscle_profile';

    my $feats_all;

    $debug && print STDERR "$subname: \$refseqs=\n". Dumper($refseqs) . "End of \$refseqs\n";
    $debug && print STDERR "$subname: \$inseqs=\n". Dumper($inseqs) . "End of \$inseqs\n";
    my $feats;
    my $alns = [];
    my $cds_all = [];
    my $vfile1 = IO::String->new('');
    my $seq_out = Bio::SeqIO->new(
                          '-fh'     => $vfile1,
                          '-format' => 'fasta'
                          );

    my $n_sets; # holds the number of [CDS, mat_peptides] ... ] in each genome. Since each CDS is unique, separate MUSCLE run has to be carried out
    foreach my $key (keys %$refseqs) {
        my $refseq = $refseqs->{$key};
        $debug && print STDERR "$subname: \$refseq=\n".Dumper($refseq)."End of \$refseq\n";
        $n_sets = $#{@$refseq};
        $debug && print STDERR "$subname: \$n_sets=\n".Dumper($n_sets)."End of \$n_sets\n";
    }
    $debug && print STDERR "$subname: \$n_sets=$n_sets\n";
    if ($aln_fn->[0] && $#{@$aln_fn}!=$n_sets) {
        $debug && print STDERR "$subname: \$n_sets=$n_sets not equal to # of aln files $#{@$aln_fn}\n";
        $debug && print STDERR "$subname: please check aln files, and re-try\n";
        return undef;
    }

    # for each unique CDS
    for my $n_set (0 .. $n_sets) {
      $debug && print STDERR "$subname: \$n_set=$n_set\n";
      my $aln;
      if ($aln_fn->[$n_set] && -f $aln_fn->[$n_set] && $debug) {
        print STDERR "$subname: \$n_set=$n_set alignment read from $aln_fn->[$n_set]\n";
        my $str = Bio::AlignIO->new('-file' => $aln_fn->[$n_set]);
        $debug && print STDERR "$subname: \$str=\n".Dumper($str)."End of \$str\n";
        $aln = $str->next_aln();
        $debug && print STDERR "$subname: \$aln=\n".Dumper($aln)."End of \$aln\n";

      } else {
        my @param;
        my $factory;
        @param = (
                   '-stable' => '',
                   '-outfile_name' => "test.afa",
                 );
        $factory = Bio::Tools::Run::Alignment::Muscle->new(@param);

        # First, get the CDS from refseq
        foreach my $key (keys %$refseqs) {
            my $refcds = $refseqs->{$key}->[$n_set]->[1];
            $debug && print STDERR "$subname: \$n_set=$n_set \$refcds=\n". Dumper($refcds) . "End of \$refcds\n";

            my @values = $refcds->get_tag_values('translation');
            my $s1 = $values[0];
            $s1 = $1 if ($s1 =~ /([^*]+)[*]$/); # to remove any trailing * in CDS
            my $str1 = 'refseq='.$refcds->seq->accession_number;
            $str1 .= '|'. $refcds->primary_tag.'='.$refcds->location->to_FTstring;
#            $debug && print "$subname: \$str=$str\n";
            my $f1 = Bio::PrimarySeq->new(
                         -seq      => $s1,
                         -id       => $str1,	# id can't contain space
                         -alphabet => 'protein'
                                    );

            push @$cds_all, $f1;

          # Second, add CDS from target genomes
          for (my $i = 0; $i<=$#{@$inseqs}; $i++) {
#            my $cds_set = $inseqs->[$i]->[$n_set+1]->{'GI:22129793'};
            my $cds_set = $inseqs->[$i]->[$n_set+1]->{$refseqs->{$key}->[$n_set]->[0]};
            $debug && print STDERR "$subname: \$cds_set=\n". Dumper($cds_set) . "End of \$cds_set\n";
            foreach my $j (0 .. $#{@$cds_set}) {
              my $cds = $cds_set->[$j]->[0];
              next if (!$cds);
              $debug && print STDERR "$subname: \$j=$j \$cds=\n". Dumper($cds) . "End of \$cds\n";
              my $acc = $cds->seq->accession_number;
#              my $s3 = $cds->seq->translate->seq;
              my $s3 = Annotate_Util::get_new_translation( $cds, $cds);
              if (!defined $s3) {
                  print STDERR "$subname: sub get_new_translation returned undef result. Skip \$j=$j accession=$acc\n";
                  print STDOUT "$subname: sub get_new_translation returned undef result. Skip \$j=$j accession=$acc\n";
                  next;
              }
              $s3 = $1 if ($s3 =~ /([^*]+)[*]$/); # to remove any trailing * in CDS
              my @values = $cds->get_tag_values('translation');
              my $s2 = $values[0];
              $s2 = $1 if ($s2 =~ /([^*]+)[*]$/); # to remove any trailing * in CDS
              if ($s2 ne $s3) {
                print STDERR "$subname: \$s2='$s2'\n";
                print STDERR "$subname: \$s3='$s3'\n";
                print STDERR "$subname: translation tag and translate don't match in CDS of ".$cds->seq->accession_number.". Skip.\n";
                print STDOUT "$subname: translation tag and translate don't match in CDS of ".$cds->seq->accession_number.". Skip.\n";
                next;
              }
              my $str2 = 'ACC='.$cds->seq->accession_number;
              $str2 .= '|'. $cds->primary_tag.'='.$cds->location->to_FTstring;
              my $f2 = Bio::PrimarySeq->new(
                         -seq      => $s2,
                         -id       => $str2,	# id can't contain space
                         -alphabet => 'protein'
                                    );
              push @$cds_all, $f2;
            }
          }

        }
        $debug && print STDERR "$subname: Finished reading CDS for \$n_set=$n_set\n";
        $debug && print STDERR "$subname: \$cds_all=\n". Dumper($cds_all) . "End of \$cds_all\n";

        # Run MUSCLE for each CDS in refseq
        # Returns a SimpleAlign object
        print STDERR "$subname: \$cds_all has $#{@$cds_all}+1 genomes for \$n_set=$n_set\n";
        $debug && print STDOUT "\n$subname: \$cds_all has $#{@$cds_all}+1 genomes for \$n_set=$n_set\n";
        if ($#{@$cds_all}<1) {
            print STDERR "$subname: Not enough genomes in \$cds_all $#{@$cds_all}+1. Skip \$n_set=$n_set\n";
            print STDERR "$subname: \$cds_all=\n". Dumper($cds_all) . "End of \$cds_all\n";
            next;
        }
        #  Returns : Reference to a SimpleAlign object containing the sequence alignment
        $aln = $factory->align(
                         [@$cds_all]
                         );
#        $debug && print STDERR "$subname: \$aln=\n". Dumper($aln) . "End of \$aln\n\n";

      }

      next if (!$aln);
      push @$alns, $aln;

      $feats = &MSA_annotate($refseqs, $inseqs, $n_set, $aln,$exe_dir);

      push @$feats_all, $feats;

    } # for my $n_set (0 .. $n_sets)

  return $feats_all;

} # sub muscle_profile


=head2 MSA_annotate

Takes a hash of refseqs, in the form of { accession => [ [CDS, mat_peptides] ... ], []... }
 and an array of target genomes, in the form of [ accession, [ [CDS, mat_peptides, ...], []... ] ]
 and alignment from MUSCLE
 Return an array of alignments CDS1, CDS2, CDS3 ...

=cut

sub MSA_annotate {
    my ($refseqs, $inseqs, $n_set, $aln,$exe_dir) = @_;

    my $debug = 0 && $debug_all;
    my $subname = 'MSA_annotate';

    my $refset = [];    # [CDS, mat_peptides] ... ]

    # get $refset
    foreach my $key (sort keys %$refseqs) {
        $debug && print STDERR "$subname: \$key=$key\n";
        $refset = $refseqs->{$key}->[$n_set];
#        $debug && print STDERR "$subname: \$refset=\n". Dumper($refset) . "End of \$refset\n";
#        my $refcds = $refseqs->{$key}->[$n_set]->[0];
#        $debug && print STDERR "$subname: \$refcds=". $refcds->seq->accession_number ."\n";
#        $debug && print STDERR "$subname: \$refcds=\n". Dumper($refcds) ."End of \$refcds\n";
        last;
    }

    # get $inset
    my $feats_all;
    my $inset = [];    # [CDS, mat_peptides] ... ]
=head1
          for (my $i = 0; $i<=$#{@$inseqs}; $i++) {
#            my $cds_set = $inseqs->[$i]->[$n_set+1]->{'GI:22129793'};
            my $cds_set = $inseqs->[$i]->[$n_set+1]->{$refseqs->{$key}->[$n_set]->[0]};
            $debug && print STDERR "$subname: \$cds_set=\n". Dumper($cds_set) . "End of \$cds_set\n";
            foreach my $j (0 .. $#{@$cds_set}) {
              my $cds = $cds_set->[$j]->[0];
              next if (!$cds);
=cut
    for (my $i = 0; $i<=$#{@$inseqs}; $i++) {
      my $inseq = $inseqs->[$i];
      for my $key (sort keys %$refseqs) {
        my $refcds_id = $refseqs->{$key}->[$n_set]->[0];
        $debug && print STDERR "$subname: \$refcds_id=$refcds_id\n";
        next if (!exists($inseq->[$n_set+1]->{$refcds_id}));

        $refset = $refseqs->{$key}->[$n_set];
        $debug && print STDERR "$subname: \$i=$i \$refset=\n". Dumper($refset) . "End of \$refset\n";

        my $cds_sets = $inseq->[$n_set+1]->{$refcds_id};
        $debug && print STDERR "$subname: \$i=$i \$cds_sets=\n". Dumper($cds_sets) . "End of \$cds_sets\n";
        foreach $inset (@$cds_sets) {
#          $inset  = $inseqs->[$i]->[$n_set+1]->[0];
          $debug && print STDERR "$subname: \$i=$i \$inset=\n". Dumper($inset) . "End of \$inset\n";

          # find the new annotations
          my $feats_new = Annotate_Muscle::MSA_annotate_1cds( $refset, $inset, $aln,$exe_dir);

          push @$feats_all, $feats_new;

#        $debug && print STDERR "$subname: \$feats_new=\n". Dumper($feats_new) . "End of \$feats_new\n";
        }
      }
    }
    $debug && print STDERR "$subname: \$feats_all=\n". Dumper($feats_all) . "End of \$feats_all\n";

    return $feats_all;
} # sub MSA_annotate


=head2 MSA_annotate_1cds

Takes a hash of refseqs, in the form of { accession => [ [CDS, mat_peptides] ... ], []... }
 and an array of target genomes, in the form of [ accession, [ [CDS, mat_peptides, ...], []... ] ]
 and alignment from MUSCLE
 Return an array of alignments CDS1, CDS2, CDS3 ...

=cut

sub MSA_annotate_1cds {
    my ($refset, $inset, $aln,$exe_dir) = @_;

    my $debug = 0 && $debug_all;
    my $subname = 'MSA_annotate_1cds';

    my $refcds = $refset->[1];
    $debug && print STDERR "$subname: \$refcds='".$refcds->seq->accession_number."'\n";
#    $debug && print STDERR "MSA_annotate_1cds: \$refcds=\n". Dumper($refcds) . "End of \$refcds\n";

    my @values = $refcds->get_tag_values('translation');
    my $s1 = $values[0];
    $s1 = $1 if ($s1 =~ /([^*]+)[*]$/); # to remove any trailing * in CDS
    my $refcds_id = 'refseq='.$refcds->seq->accession_number;
    $refcds_id .= '|'. $refcds->primary_tag.'='.$refcds->location->to_FTstring;


    my $cds = $inset->[0];
    next if (!$cds);
    my $acc = $cds->seq->accession_number;
#    $debug && print STDERR "$subname: \$cds=\n". Dumper($cds) . "End of \$cds\n";
#    $debug && print STDERR "$subname: \$inset=\n". Dumper($inset) . "End of \$inset\n";

    my $s3 = Annotate_Util::get_new_translation( $cds, $cds);
    if (!defined $s3) {
       print STDERR "$subname: sub get_new_translation returned undef result. Skip accession=$acc\n";
       print STDOUT "$subname: sub get_new_translation returned undef result. Skip accession=$acc\n";
       return;
    }
    $s3 = $1 if ($s3 =~ /([^*]+)[*]$/); # to remove any trailing * in CDS
    @values = $cds->get_tag_values('translation');
    my $s2 = $values[0];
    $s2 = $1 if ($s2 =~ /([^*]+)[*]$/); # to remove any trailing * in CDS
    if ($s2 ne $s3) {
        print STDERR "$subname: \$s2='$s2'\n";
        print STDERR "$subname: \$s3='$s3'\n";
        print STDERR "$subname: translation tag and translate don't match in CDS of '$acc'.\n";
        print STDERR "$subname: Skip.\n\n";
        next;
    }
    my $cds_id = "ACC=$acc";
    $cds_id .= '|'. $cds->primary_tag.'='.$cds->location->to_FTstring;
    $debug && print STDERR "$subname: \$refcds_id='$refcds_id' \$cds_id='$cds_id'\n";

    my $note = "Annotated by VIPRBRC, MSA, refseq=".$refcds->seq->accession_number();
    my $feats_all;
    $feats_all = Annotate_Util::project_matpept(
                         $refset,
                         $inset,
                         $aln,
                         $refcds_id,
                         $cds_id,
                         $note,$exe_dir,
                         );

    $debug && print STDERR "$subname: \$feats_all=\n". Dumper($feats_all) . "End of \$feats_all\n";

    return $feats_all;
} # sub MSA_annotate_1cds


1;
