package Annotate_Muscle;

use strict;
use warnings;
use English;
use Carp;
use Data::Dumper;

use version; our $VERSION = qv('1.1.4'); # May 09 2012
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
    my ($gbk, $exe_dir, $aln_fn, $dir_path) = @_;

    my $debug = 0 && $debug_all;
    my $subname = 'annotate_1gbk';

    $aln_fn = [] if (!defined($aln_fn));
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
    if ( !defined($refpolyprots) ) {
            $comment = $acc." w/ taxid=".$inseq->species->ncbi_taxid." is not covered in V$VERSION";
            return (undef, $comment);
    }

    $debug && print STDERR "$subname: \$refpolyprots = $#{$refpolyprots}\n";
    $debug && print STDERR "$subname: \$refpolyprots = \n".Dumper($refpolyprots)."End of \$refpolyprots\n\n";

    # According to refseq, get the CDS/mat_peptides in inseq, use bl2seq to determine if the CDS matches
    my $num_cds;
    ($polyprots, $num_cds, $comment) = Annotate_Util::get_polyprots( $inseq, $refpolyprots);
    $debug && print STDERR "$subname:    \$num_cds = $num_cds\n";
    $debug && print STDERR "$subname: \$polyprots = \n".Dumper($polyprots)."End of \$polyprots\n\n";
    # Skip the refseqs that has mat_peptide annotation from genbank
    if ($acc =~ /^NC_\d+$/i) {
        my $has_mat_peptide = 0;
        CHECK: for my $key (keys %$polyprots) {
            my $sets = $polyprots->{$key};
#            $debug && print STDERR "$subname: \$sets=\n".Dumper($sets)."End of \$sets\n\n";
            for my $j (0 .. $#{$sets}) {
                my $set = $sets->[$j];
#                $debug && print STDERR "$subname: \$set=\n".Dumper($set)."End of \$set\n\n";
                for my $k (1 .. $#{$set}) {
                    next if ($set->[$k]->primary_tag ne 'mat_peptide');
                    $has_mat_peptide = 1;
                    last CHECK;
                }
            }
        }
        if ( 0 && $has_mat_peptide && !$debug ) {
            $comment = 'Refseq with mat_peptide annotation from NCBI, skip';
            return (undef, $comment);
        }
    }

    # add refseq to hash, add inseq to array
    if ($#{$refpolyprots} >=0) {
        if (!exists($refseqs->{$refpolyprots->[0]->[1]->seq->accession_number})) {
            $refseqs->{$refpolyprots->[0]->[1]->seq->accession_number} = $refpolyprots;
        }
    } else {
        $comment = 'No suitable polyprotein was found for taxid='. $inseq->species->ncbi_taxid;
        return (undef, $comment);
    }
    if ($num_cds<0) {
        print STDERR "$subname: no polyprotein (suitable CDS) found in ".$inseq->accession_number." \$polyprots=$num_cds\n";
        $comment = "No polyprotein: ". $comment;
        return (undef, $comment);
    }
    push @$inseqs, [$inseq->accession_number, $polyprots];
    $debug && print STDERR "$subname: \$inseqs=\n". Dumper($inseqs) . "End of \$inseqs\n";

    my $feats_all;
    $feats_all = Annotate_Muscle::muscle_profile( $refseqs, $inseqs, $aln_fn,$exe_dir, $dir_path);
    $debug && print STDERR "$subname: \$feats_all=\n". Dumper($feats_all) . "End of \$feats_all\n";

    my $n = 0;
    my $key1;
    for my $key (keys %$feats_all) {
        $key1 = $key;
        my $feats = $feats_all->{$key};
        my $refcds_ids = [ keys %$feats ];
        $debug && print STDERR "$subname: \$refcds_ids=\n". Dumper(@$refcds_ids) . "End of \$refcds_ids\n";
        for my $k2 (@$refcds_ids) {
          $debug && print STDERR "$subname: \$key=$key \$k2=$k2 \$feats_all=\n". Dumper($feats_all->{$key}->{$k2}) . "End of \$feats_all\n";
          $n += $#{$feats_all->{$key}->{$k2}};
        }
    }
    $comment = "MSA returned $n mat_peptides";
    $debug && print STDERR "$subname: \$key1=$key1\n";

    # In case a genbank file w/ new annotation is needed, save it here. Example: NC_004162
    if ( 0 && $debug ) {
        for my $key (keys %$feats_all) {
          $debug && print STDERR "$subname: \$key=$key \$feats_all=\n". Dumper($feats_all->{$key}) . "End of \$feats_all\n";
          for my $k2 (keys %{$feats_all->{$key}}) {
            my $feats = $feats_all->{$key}->{$k2};

            # Look for a CDS in old annotation
            my $old_feats = [ $inseq->get_all_SeqFeatures ];
            $debug && print STDERR "$subname: \$key=$key \$k2=$k2 \$old_feats=\n". Dumper($old_feats) . "End of \$old_feats\n";
            for my $i (0 .. $#{$old_feats}) {
                $debug && print STDERR "$subname: \$key=$key \$k2=$k2 \$old_feats->[$i]=". $old_feats->[$i]->primary_tag . "\n";
                next if ($old_feats->[$i]->primary_tag ne 'CDS');
                $debug && print STDERR "$subname: \$key=$key \$k2=$k2 \$old_feats->[$i]=". $old_feats->[$i]->primary_tag . "\n";
                my $cds = $old_feats->[$i];
                my $old_cds_id = '';
                if ($cds->has_tag('db_xref')) {
                  my @id = $cds->get_tag_values('db_xref');
                  for my $id1 (@id) {
                    $old_cds_id = $id1 if ($id1 =~ /^GI:/i);
                  }
                }
                my $cds_id = '';
                if ($feats->[0]->has_tag('db_xref')) {
                  my @id = $feats->[0]->get_tag_values('db_xref');
                  for my $id1 (@id) {
                    $cds_id = $id1 if ($id1 =~ /^GI:/i);
                  }
                }
                next if ($cds_id ne $old_cds_id);
                for my $j (reverse 1 .. $#{$feats}) { # Skip the CDS
                  my $seen = 0;
                  for my $k ($i+1 .. $#{$old_feats}) {
                      next if ($old_feats->[$k]->primary_tag ne 'mat_peptide');
                      my $loc = $old_feats->[$k]->location;
                      $seen = 1 if ($loc->start == $feats->[$j]->location->start && $loc->end == $feats->[$j]->location->end);
                  }
                  if (!$seen) {
                      my $new_feats = [ @{$old_feats}[0 .. $i],
                                        $feats->[$j],
                                        @{$old_feats}[$i+1 .. $#{$old_feats}] ];
                      $old_feats = $new_feats;
                  }
                  $debug && print STDERR "$subname: \$key=$key \$k2=$k2 \$old_feats=\n". Dumper($old_feats) ."End of \$old_feats\n";
                }
                $inseq->remove_SeqFeatures;
                $inseq->add_SeqFeature(@$old_feats);
            }
            $debug && print STDERR "$subname: \$key=$key \$k2=$k2 \$inseq=\n". Dumper($inseq) . "End of \$inseq\n";
          }

          my $outgbfile = $inseq->accession_number."_msaa.gb";
          my $seq_out = Bio::SeqIO->new('-file' => ">$outgbfile", '-format' => 'genbank');
          $seq_out->write_seq($inseq);
        }
    }
    return ($feats_all->{$key1}, $comment);

} # sub annotate_1gbk


=head2 muscle_profile

Takes a hash of refseqs, and an array of target genomes, in the form of
 [ accession, [CDS, mat_peptides] ... ], []... ]
 run muscle on each CDS,
 Return an array of alignments CDS1, CDS2, CDS3 ...

=cut

sub muscle_profile {
    my ($refseqs, $inseqs, $aln_fn,$exe_dir, $dir_path) = @_;

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
        $debug && print STDERR "$subname: \$key=$key \$refseq=\n".Dumper($refseq)."End of \$refseq\n";
        $n_sets = $#{$refseq};
        $debug && print STDERR "$subname: \$n_sets=$n_sets\n";
    }
    $debug && print STDERR "$subname: \$n_sets=$n_sets\n";
    if ($aln_fn->[0] && $#{$aln_fn}!=$n_sets) {
        $debug && print STDERR "$subname: \$n_sets=$n_sets not equal to # of aln files $#{$aln_fn}\n";
        $debug && print STDERR "$subname: please check aln files, and re-try\n";
        return undef;
    }


    # for each unique CDS in refseq
    for my $n_set (0 .. $n_sets) {
      $debug && print STDERR "$subname: \$n_set=$n_set\n";
      $cds_all = [];
      my $aln;
#      if ($aln_fn->[$n_set] && -f $aln_fn->[$n_set] && $debug) {
#        print STDERR "$subname: \$n_set=$n_set alignment read from $aln_fn->[$n_set]\n";
#        my $str = Bio::AlignIO->new('-file' => $aln_fn->[$n_set]);
#        $debug && print STDERR "$subname: \$str=\n".Dumper($str)."End of \$str\n";
#        $aln = $str->next_aln();
#        $debug && print STDERR "$subname: \$aln=\n".Dumper($aln)."End of \$aln\n";
#
#      } else {
        my @param;

        # First, get the CDS from refseq
        foreach my $key (keys %$refseqs) {
            my $refcds = $refseqs->{$key}->[$n_set]->[1];
            $debug && print STDERR "$subname: \$n_set=$n_set \$key=$key \$refcds=\n". Dumper($refcds) . "End of \$refcds\n";

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

#            $cds_all = []; # Clear all cds, run MUSCLE separately for each input CDS
#            push @$cds_all, $f1;

            # Second, add CDS from target genomes
            for (my $i = 0; $i<=$#{$inseqs}; $i++) {
#                $cds_all = []; # Clear all cds, run MUSCLE separately for each input CDS
#                push @$cds_all, $f1;
#                my $cds_set = $inseqs->[$i]->[$n_set+1]->{'GI:22129793'};
                $debug && print STDERR "$subname: \$key=$key \$i=$i \$inseqs=\n". Dumper($inseqs) . "End of \$inseqs\n";
                $debug && print STDERR "$subname: \$key=$key \$i=$i \$n_set=$n_set \$key=$key $refseqs->{$key}->[$n_set]->[0]\n";
                if (!exists($inseqs->[$i]->[1]->{$refseqs->{$key}->[$n_set]->[0]})) {
                  $debug && print STDERR "$subname: $refseqs->{$key}->[$n_set]->[0] not found in \$inseqs\n";
                  $debug && print STDERR "$subname: \$key=$key \$i=$i \$inseqs=\n". Dumper($inseqs) . "End of \$inseqs\n";
                  next;
                }
#                my $cds_set = $inseqs->[$i]->[$n_set+1]->{$refseqs->{$key}->[$n_set]->[0]};
                my $cds_set = $inseqs->[$i]->[1]->{$refseqs->{$key}->[$n_set]->[0]};
                $debug && print STDERR "$subname: \$key=$key \$i=$i \$cds_set=\n". Dumper($cds_set) . "End of \$cds_set\n";
                $debug && print STDERR "$subname: \$key=$key \$i=$i \$inseqs=\n". Dumper($inseqs) . "End of \$inseqs\n";
                foreach my $j (0 .. $#{$cds_set}) {
                  $cds_all = []; # Clear all cds, run MUSCLE separately for each input CDS
                  push @$cds_all, $f1;
                  my $cds = $cds_set->[$j]->[0];
                  next if (!$cds);
                  $debug && print STDERR "$subname: \$j=$j \$cds=\n". Dumper($cds) . "End of \$cds\n";
                  my $acc = $cds->seq->accession_number;
#                  my $s3 = $cds->seq->translate->seq;
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
                  $s3 = s/[*]/./g if $s3 =~ /[*]/;
#                  if ($s2 ne $s3) {
                  if ($s2 !~ /$s3/) {
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
#                } # foreach my $j (0 .. $#{$cds_set})
#            } # for (my $i = 0; $i<=$#{$inseqs}; $i++) {

#        } # foreach my $key (keys %$refseqs) {
#        $debug && print STDERR "$subname: Finished reading CDS for \$n_set=$n_set\n";
                $debug && print STDERR "$subname: \$cds_all=\n". Dumper($cds_all) . "End of \$cds_all\n";

        # Run MUSCLE for each CDS in refseq
        # Returns a SimpleAlign object
                print STDERR "$subname: \$cds_all has $#{$cds_all}+1 genomes for \$n_set=$n_set\n";
                $debug && print STDOUT "\n$subname: \$cds_all has $#{$cds_all}+1 genomes for \$n_set=$n_set\n";
                for my $iii ( 0 .. $#{$cds_all}) {
                    print STDERR "$subname: \$cds_all->[$iii]=".$cds_all->[$iii]->display_id."\n";
                }
                if ($#{$cds_all}<1) {
                    print STDERR "$subname: Not enough sequences in \$cds_all $#{$cds_all}+1. Skip \$n_set=$n_set\n";
                    print STDERR "$subname: \$cds_all=\n". Dumper($cds_all) . "End of \$cds_all\n";
                    next;
                }
                #  Returns : Reference to a SimpleAlign object containing the sequence alignment
                my $factory;
                my $outfile_name = "$dir_path/test.afa";
                if ($debug) {
                    my $display_id = $cds_all->[1]->display_id;
                    my $display_id2 = '';
                    $display_id2 = $1 if $display_id =~ /ACC=(.+)[|]/;
                    if ($display_id =~ /\|CDS=\D*(\d+[.]).*([.][<>]*\d+)\D*$/) {
                        $display_id2 .= "_${1}${2}";
                    }
#                    $display_id =~ s/\|CDS=/_/i;
                    $display_id2 =~ s/[.]{2}/_/i;
                    $display_id2 =~ s/[<>]//i;
                    my $count = 1;
                    while (-f $outfile_name && $count <10) {
                        $debug && print STDERR "$subname: \$display_id2=$display_id2 \$count=$count\n";
                        $outfile_name = sprintf("$dir_path/test_%s_%02d.afa", $display_id2, $count);
                        $count++;
                    }
                }
                @param = (
                   '-stable' => '',
                   '-outfile_name' => "$outfile_name",
                         );
                $factory = Bio::Tools::Run::Alignment::Muscle->new(@param);
                $aln = $factory->align(
                         [@$cds_all]
                         );
                $debug && print STDERR "$subname: `cat $outfile_name`\n". `cat $outfile_name` . "End of $outfile_name\n";
#                $debug && print STDERR "$subname: \$aln=\n". Dumper($aln) . "End of \$aln\n\n";

#      } # if

      next if (!$aln);
      push @$alns, $aln;

              $feats = &MSA_annotate($refseqs, $inseqs, $n_set, $aln, $exe_dir);
#              print STDERR "$subname: \$n_set=$n_set MSA returned \$feats=$#{$feats}\n";
              $debug && print STDERR "$subname: \$n_set=$n_set \$feats=\n". Dumper($feats) . "End of \$feats\n";

#              push @$feats_all, $feats;
              # reverse the indexes
              for my $key (keys %$feats) {
                  for my $k2 (keys %{$feats->{$key}}) {
                    if (exists($feats_all->{$k2}->{$key})) {
                      for (my $kk = 1; $kk<=$#{$feats->{$key}->{$k2}}; $kk++) {
                        my $feat_kk = $feats->{$key}->{$k2}->[$kk];
                        my $product_kk = [ $feat_kk->get_tag_values('product') ];
                        $product_kk = $product_kk->[0];
                        $debug && print STDERR "$subname: \$key=$key \$k2=$k2 \$product_kk=$product_kk ".$feat_kk->location->to_FTstring." is being checked\n";
                        my $seen = 0;
                        my $feat_ff;
                        # Check if the product name has been seen
                        for my $ff (1 .. $#{$feats_all->{$k2}->{$key}}) {
                          $feat_ff = $feats_all->{$k2}->{$key}->[$ff];
                          my $product_f = [ $feat_ff->get_tag_values('product') ];
                          $product_f = $product_f->[0];
                          $debug && print STDERR "$subname: \$ff=$ff \$product_f =$product_f ". $feat_ff->location->to_FTstring . "\n";
#                          $debug && print STDERR "$subname: \$ff=$ff \$f=\n". Dumper($f) . "End of $f\n";
                          next if ($product_f ne $product_kk);
                          $seen = 1;
                          my $len_kk = $feat_kk->location->end   - $feat_kk->location->start +1;
                          my $len_ff = $feat_ff->location->end   - $feat_kk->location->start +1;
                          $debug && print STDERR "$subname: \$len_kk=$len_kk \$len_ff=$len_ff\n";
                          if ( $len_kk > $len_ff ) {
                            $feats_all->{$k2}->{$key}->[$ff] = $feat_kk;
                          }
                          last;
                        }
                        $debug && print STDERR "$subname: \$key=$key \$k2=$k2 ". $feats->{$key}->{$k2}->[$kk]->location->to_FTstring . " \$seen=$seen\n";
                        if ($seen) {
                          next;
                        }
                        $debug && print STDERR "$subname: \$key=$key \$k2=$k2 ". $feats->{$key}->{$k2}->[$kk]->location->to_FTstring . " \$seen=$seen, keep\n";
                        push @{$feats_all->{$k2}->{$key}}, $feats->{$key}->{$k2}->[$kk];
                      }
                    } else {
                      $feats_all->{$k2}->{$key} = $feats->{$key}->{$k2};
                    }
                  }
              }
              $debug && print STDERR "$subname: \$n_set=$n_set \$feats_all=\n". Dumper($feats_all) . "End of \$feats_all\n";
                } # foreach my $j (0 .. $#{$cds_set})
            } # for (my $i = 0; $i<=$#{$inseqs}; $i++) {

        } # foreach my $key (keys %$refseqs) {

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
    my $feats_all = {};
    my $inset = [];    # [CDS, mat_peptides] ... ]

    for (my $i = 0; $i<=$#{$inseqs}; $i++) {
      my $inseq = $inseqs->[$i];
      for my $key (sort keys %$refseqs) {
        my $refcds_id = $refseqs->{$key}->[$n_set]->[0];
        $debug && print STDERR "$subname: \$refcds_id=$refcds_id\n";
        next if (!exists($inseq->[1]->{$refcds_id}));

        $refset = $refseqs->{$key}->[$n_set];
        $debug && print STDERR "$subname: \$i=$i \$refset=\n". Dumper($refset) . "End of \$refset\n";

#        my $cds_sets = $inseq->[$n_set+1]->{$refcds_id};
        my $cds_sets = $inseq->[1]->{$refcds_id};
        $debug && print STDERR "$subname: \$i=$i \$cds_sets=\n". Dumper($cds_sets) . "End of \$cds_sets\n";
        foreach $inset (@$cds_sets) {
#          $inset  = $inseqs->[$i]->[$n_set+1]->[0];
          $debug && print STDERR "$subname: \$i=$i \$inset=\n". Dumper($inset) . "End of \$inset\n";

          # find the new annotations
          my $feats_new = Annotate_Muscle::MSA_annotate_1cds( $refset, $inset, $aln,$exe_dir);
          $debug && print STDERR "$subname: \$feats_new=\n". Dumper($feats_new) . "End of \$feats_new\n";
          $debug && print STDERR "$subname: \$feats_all=\n". Dumper($feats_all) . "End of \$feats_all\n";

#          push @$feats_all, $feats_new;
          my $acc = $inset->[0]->seq->accession_number;
          $debug && print STDERR "$subname: \$refcds_id=$refcds_id \$acc=$acc\n";
          if (!exists($feats_all->{$refcds_id}) && !exists($feats_all->{$refcds_id}->{$acc})) {
              $debug && print STDERR "$subname: \$refcds_id=$refcds_id found new mat_peptide=". $#{$feats_new} . "\n";
              $feats_all->{$refcds_id}->{$acc} = $feats_new;
          } else {
              for (my $kk = 0; $kk<=$#{$feats_new}; $kk++) {
                  # Skip the mat_peptide if it's already in the list
#                  my $seen = 0;
#                  $debug && print STDERR "$subname: \$refcds_id=$refcds_id ". $feats_new->[$kk]->location->to_FTstring . " is being checked\n";
#                  for my $ff (0 .. $#{$feats_all->{$refcds_id}->{$acc}}) {
#                    my $f = $feats_all->{$refcds_id}->{$acc}->[$ff];
#                    $debug && print STDERR "$subname: \$ff=$ff ". $f->location->to_FTstring . " is being checked\n";
#                    next if ($f->location->to_FTstring ne $feats_new->[$kk]->location->to_FTstring);
#                    $debug && print STDERR "$subname: \$kk=$kk \$f=". $f->location->to_FTstring . " already seen, skip\n";
#                    $seen = 1;
#                    last;
#                  }
#                  $debug && print STDERR "$subname: \$kk=$kk ". $feats_new->[$kk]->location->to_FTstring . " \$seen=$seen, keep\n";
#                  next if ($seen);
#                  $debug && print STDERR "$subname: \$kk=$kk ". $feats_new->[$kk]->location->to_FTstring . " is new, keep\n";
                  push @{$feats_all->{$refcds_id}->{$acc}}, $feats_new->[$kk];
              }
          }

        $debug && print STDERR "$subname: \$feats_all=\n". Dumper($feats_all) . "End of \$feats_all\n";
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
    $s3 = s/[*]/./g if $s3 =~ /[*]/;
#    if ($s2 ne $s3) {
    if ($s2 !~ /$s3/) {
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
    my $feats_all = [];
    my $aln_h = Annotate_Util::msa_get_aln( $aln, $cds_id);
    $debug && print STDERR "$subname: \$aln_h=\n". Dumper($aln_h) . "End of \$aln_h\n";
#    if ( !$aln_h || !$aln_h->isa('Bio::LocatableSeq')) {
    if ( !$aln_h ) {
        $debug && print STDERR "$subname: \$aln_h is empty\n";
        return $feats_all;
    }
    $feats_all = Annotate_Util::project_matpept(
                         $refset,
                         $inset,
                         $aln,
                         $refcds_id,
                         $cds_id,
                         $note,$exe_dir,
                         );

    $debug && print STDERR "$subname: \$feats_all=$#{$feats_all}\n". Dumper($feats_all) . "End of \$feats_all\n";

    return $feats_all;
} # sub MSA_annotate_1cds


1;
