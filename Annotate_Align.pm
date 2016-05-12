package Annotate_Align;

use strict;
use warnings;
use English;
use Carp;
use Data::Dumper;

use version; our $VERSION = qv('1.2.0'); # Apr 12 2013
use Bio::SeqIO;
use Bio::Seq;
use Bio::AlignIO;
use Bio::Tools::Run::Alignment::Muscle;
use Bio::Tools::Run::Alignment::Clustalw;
use IO::String;
use Annotate_Math;
use Annotate_Util;
use Annotate_Verify;

my $debug_all = 0;

####//README//####
#
# Annotate_Align contains the core functions to perform annotation based on MUSCLE alignment.
# The alignment is done with AA sequences of ref CDS and target CDS
#
#    Authors Chris Larsen, clarsen@vecna.com; Guangyu Sun, gsun@vecna.com
#    December 2010
#
##################


sub setDebugAll {
    my ($debug) = @_;
    $debug_all = $debug;
} # sub setDebugAll


=head2 saveNewGenbank
Takes the newly annotated mat_peptides, saves to a genbank file
 Return the number of saved files, or 0
=cut

sub saveNewGenbank {
    my ($inseq, $feats_all) = @_;

    my $debug = 0 || $debug_all;
    my $subn = 'saveNewGenbank';

    my $acc = $inseq->accession_number;

    # In case a genbank file w/ new annotation is needed, save it here. Example: NC_004162
    print STDERR "\n$subn: saving annotated mat_peptides to gbk files...\n";
    my $ct = 0;
    for my $key (keys %$feats_all) {
        $ct++;
        $debug && print STDERR "$subn: \$key=$key \$feats_all=\n". Dumper($feats_all->{$key}) . "End of \$feats_all\n";
        for my $k2 (keys %{$feats_all->{$key}}) {
            my $feats = $feats_all->{$key}->{$k2};

            # Look for a CDS in old annotation
            my $old_feats = [ $inseq->get_all_SeqFeatures ];
            $debug && print STDERR "$subn: \$key=$key \$k2=$k2 \$old_feats=\n". Dumper($old_feats) . "End of \$old_feats\n";
            for my $i (0 .. $#{$old_feats}) {
                $debug && print STDERR "$subn: \$key=$key \$k2=$k2 \$old_feats->[$i]=". $old_feats->[$i]->primary_tag . "\n";
                next if ($old_feats->[$i]->primary_tag ne 'CDS');
                $debug && print STDERR "$subn: \$key=$key \$k2=$k2 \$old_feats->[$i]=". $old_feats->[$i]->primary_tag . "\n";
                my $cds = $old_feats->[$i];
                my $old_cds_id = '';
                if ($cds->has_tag('db_xref')) {
                  my @id = $cds->get_tag_values('db_xref');
                  for my $id1 (@id) {
                    $old_cds_id = $id1 if ($id1 =~ /^GI:/i);
                  }
                }
                my $cds_id = '';
if ( 0 ) {
                if ($feats->[0]->has_tag('db_xref')) {
                  my @id = $feats->[0]->get_tag_values('db_xref');
                  for my $id1 (@id) {
                    $cds_id = $id1 if ($id1 =~ /^GI:/i);
                  }
                }
} else {
                  my @id = $feats->[0]->get_tag_values('note'); # 4/09/2013
                  for my $id1 (@id) {
                    $cds_id = $1 if ($id1 =~ /[|]CDS=([^|]+)[|]/i);
                  }
}
                $debug && print STDERR "$subn: \$key=$key \$k2=$k2 \$old_feats->[$i]=$old_cds_id \$cds_id=$cds_id\n";
                next if ($cds_id ne $old_cds_id);
                for my $j (reverse 0 .. $#{$feats}) { # All feats are mat_peptide, CDS is not included anymore
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
                  $debug && print STDERR "$subn: \$key=$key \$k2=$k2 \$old_feats=\n". Dumper($old_feats) ."End of \$old_feats\n";
                }
                $inseq->remove_SeqFeatures;
                $inseq->add_SeqFeature(@$old_feats);
            }
            $debug && print STDERR "$subn: \$key=$key \$k2=$k2 \$inseq=\n". Dumper($inseq) . "End of \$inseq\n";
        }

        my $outgbfile = $inseq->accession_number."_msaa.gb";
        my $seq_out = Bio::SeqIO->new('-file' => ">$outgbfile", '-format' => 'genbank');
        $seq_out->write_seq($inseq);
        print STDERR "$subn: File #$ct: saved to $outgbfile\n";
    }
    print STDERR "$subn: Saved $ct files\n";
    print STDERR "$subn: Finished, exit\n\n";
    
    return $ct;

} # sub saveNewGenbank


=head2 annotate_1gbk
Takes a genbank file, finds hash of refseqs, and an array of target genomes, in the form of
 [ accession, [CDS, mat_peptides] ... ], []... ]
 run muscle on each CDS,
 Return an array of alignments CDS1, CDS2, CDS3 ...
=cut

sub annotate_1gbk {
    my ($gbk, $exe_dir, $aln_fn, $dir_path, $progs) = @_;

    my $debug = 0 || $debug_all;
    my $subn = 'annotate_1gbk';

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
    $refpolyprots = Annotate_Def::get_refpolyprots( $refseqs, $inseq, $exe_dir);
    if ( !defined($refpolyprots) ) {
            my $species = Annotate_Def::getTaxonInfo( $inseq->species->ncbi_taxid);
            my $fam = ($species->[6]) ? $species->[6] : 'unknown';
            $species = ($species->[4]) ? $species->[4] : '';
            $comment = $acc." w/ taxid=".$inseq->species->ncbi_taxid."($fam:$species)";
            $comment .= " not covered in V$Annotate_Def::VERSION";
            return (undef, $comment);
    }
    my $n_reffeat = 0;
    for my $set (@$refpolyprots) {
        for my $i (1 .. $#{$set}) {
        $n_reffeat++ if ($set->[$i]->primary_tag eq 'mat_peptide');
        }
    }
    print STDERR "$subn: \$refpolyprots=$#{$refpolyprots} \$n_reffeat=$n_reffeat\n";
    $debug && print STDERR "$subn: \$refpolyprots=\n".Dumper($refpolyprots)."End of \$refpolyprots\n\n";

    # According to refseq, get the CDS/mat_peptides in inseq, use bl2seq to determine if the CDS matches
    my $num_cds;
    ($polyprots, $num_cds, $comment) = Annotate_Util::get_polyprots( $inseq, $refpolyprots);
    $debug && print STDERR "$subn:  ".$inseq->accession_number."  \$num_cds = $num_cds\n";
    $debug && print STDERR "$subn: \$polyprots=\n".Dumper($polyprots)."End of \$polyprots\n\n";
    # Skip the refseqs that has mat_peptide annotation from genbank
    if ($acc =~ /^NC_\d+$/i) {
        my $has_mat_peptide = 0;
        CHECK: for my $key (keys %$polyprots) {
            my $sets = $polyprots->{$key};
#            $debug && print STDERR "$subn: \$sets=\n".Dumper($sets)."End of \$sets\n\n";
            for my $j (0 .. $#{$sets}) {
                my $set = $sets->[$j];
#                $debug && print STDERR "$subn: \$set=\n".Dumper($set)."End of \$set\n\n";
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
        print STDERR "$subn: no polyprotein (suitable CDS) found in ".$inseq->accession_number." \$polyprots=$num_cds\n";
        $comment = "No polyprotein: ". $comment;
        return (undef, $comment);
    }
    push @$inseqs, [$inseq->accession_number, $polyprots];
    $debug && print STDERR "$subn: \$inseqs=\n". Dumper($inseqs) . "End of \$inseqs\n";

    my $feats_all;
    $feats_all = Annotate_Align::run_MSA( $refseqs, $inseqs, $aln_fn,$exe_dir, $dir_path, $progs);
    $debug && print STDERR "$subn: \$feats_all=\n". Dumper($feats_all) . "End of \$feats_all\n";

    my $n = 0;
    my $key1;
    for my $key (keys %$feats_all) {
        $key1 = $key;
        my $feats = $feats_all->{$key};
        my $refcds_ids = [ keys %$feats ];
        $debug && print STDERR "$subn: \$refcds_ids=\n". Dumper(@$refcds_ids) . "End of \$refcds_ids\n";
        for my $k2 (@$refcds_ids) {
          $debug && print STDERR "$subn: \$key=$key \$k2=$k2 \$feats_all=\n". Dumper($feats_all->{$key}->{$k2}) . "End of \$feats_all\n";
          $n += $#{$feats_all->{$key}->{$k2}} +1; # array $feats_all->{$key}->{$k2} contains CDS as 0th element
        }
    }
    $comment = "MSA returned $n/$n_reffeat mat_peptides vs refseq";
    $debug && print STDERR "$subn: \$key1=$key1\n";

    # In case a genbank file w/ new annotation is needed, save it here. Example: NC_004162
    my $SAVE_GBK_ANNOTATED = 0; # Saves a new version of genbank with annotated mat_peptide, then exit
    if ( $SAVE_GBK_ANNOTATED ) {
        Annotate_Align::saveNewGenbank( $inseq, $feats_all);
        exit(0);
    }

    return ($feats_all->{$key1}, $comment);

} # sub annotate_1gbk

sub msa_seqid {
    my ($str1, $cds) = @_;

    my $debug = 0 && $debug_all;
    my $subn = 'msa_seqid';

    $str1 .= $cds->seq->accession_number;
    if ( 0 ) {
        $str1 .= '|'. $cds->primary_tag.'='.$cds->location->to_FTstring;
    } else {
        $str1 .= '|'. $cds->primary_tag.'='.$cds->location->start.'..'.$cds->location->end;
    }
    $debug && print STDERR "$subn: \$str1='$str1'\n";
    return $str1;
} # msa_seqid


=head2 run_MSA
Takes a hash of refseqs, and an array of target genomes, in the form of
 [ accession, [CDS, mat_peptides] ... ], []... ]
 run clustalw (or muscle) on each CDS,
 Return an array of alignments CDS1, CDS2, CDS3 ...
=cut

sub run_MSA {
    my ($refseqs, $inseqs, $aln_fn,$exe_dir, $dir_path, $progs) = @_;

    my $debug = 0 || $debug_all;
    my $subn = 'run_MSA';

    # $runMUSCLE = 0; # 0: Muscle, 1: Clustalw
    my $runMUSCLE = 0;
    if ($runMUSCLE) {
        # Check if MUSCLE is accessible either from prompt or from user-defined location
        $progs->{muscle}  = (defined($progs) && $progs->{muscle})  ? $progs->{muscle}  : "muscle";
    } else {
        # Check if CLUSTALW is accessible either from prompt or from user-defined location
        # on gsun's workstation, clustalw is a softlink to clustalw2 for v2.0
        $progs->{clustalw}  = (defined($progs) && $progs->{clustalw})  ? $progs->{clustalw}  : "clustalw";
    }
    foreach my $k (keys %$progs) {
            my $err = Annotate_Util::check_program($progs->{$k});
            $debug && print STDERR "$subn: \$err='$err'\n";
            exit(1) if ($err && print STDERR "$subn: \$err='$err'\n");
    }

    my $feats_all = {};

    $debug && print STDERR "$subn: \$refseqs=\n". Dumper($refseqs) . "End of \$refseqs\n";
    $debug && print STDERR "$subn: \$inseqs=\n". Dumper($inseqs) . "End of \$inseqs\n";
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
        $debug && print STDERR "$subn: \$key=$key \$refseq=\n".Dumper($refseq)."End of \$refseq\n";
        $n_sets = $#{$refseq};
        $debug && print STDERR "$subn: \$n_sets=$n_sets\n";
    }
    $debug && print STDERR "$subn: \$n_sets=$n_sets\n";
    if ($aln_fn->[0] && $#{$aln_fn}!=$n_sets) {
        $debug && print STDERR "$subn: \$n_sets=$n_sets not equal to # of aln files $#{$aln_fn}\n";
        $debug && print STDERR "$subn: please check aln files, and re-try\n";
        return undef;
    }


    # for each unique CDS in refseq
    my $uniqId = 0;
    for my $n_set (0 .. $n_sets) {
        $debug && print STDERR "$subn: \$n_set=$n_set\n";
        $cds_all = [];
        my $aln;
        my @param;

        # First, get the CDS from refseq
        foreach my $key (keys %$refseqs) {
            my $refcds = $refseqs->{$key}->[$n_set]->[1];
            $debug && print STDERR "$subn: \$n_set=$n_set \$key=$key \$refcds=\n". Dumper($refcds) . "End of \$refcds\n";

            my @values = $refcds->get_tag_values('translation');
            my $s1 = $values[0];
            $s1 = $1 if ($s1 =~ /([^*]+)[*]$/); # to remove any trailing * in CDS
            my $str1 = Annotate_Align::msa_seqid('ref=', $refcds);
            $debug && print STDERR "$subn: \$str1=$str1\n";
            my $f1 = Bio::PrimarySeq->new(
                         -seq      => $s1,
                         -id       => $str1,	# id can't contain space
                         -alphabet => 'protein'
                                    );

            # Second, add CDS from target genomes
            for (my $i = 0; $i<=$#{$inseqs}; $i++) {
#                my $cds_set = $inseqs->[$i]->[$n_set+1]->{'GI:22129793'};
                $debug && print STDERR "$subn: \$key=$key \$i=$i \$inseqs=\n". Dumper($inseqs) . "End of \$inseqs\n";
                $debug && print STDERR "$subn: \$key=$key \$i=$i \$n_set=$n_set \$key=$key $refseqs->{$key}->[$n_set]->[0]\n";
                if (!exists($inseqs->[$i]->[1]->{$refseqs->{$key}->[$n_set]->[0]})) {
                  $debug && print STDERR "$subn: $refseqs->{$key}->[$n_set]->[0] not found in \$inseqs\n";
                  $debug && print STDERR "$subn: \$key=$key \$i=$i \$inseqs=\n".Dumper($inseqs)."End of \$inseqs\n";
                  next;
                }
#                my $cds_set = $inseqs->[$i]->[$n_set+1]->{$refseqs->{$key}->[$n_set]->[0]};
                my $cds_set = $inseqs->[$i]->[1]->{$refseqs->{$key}->[$n_set]->[0]};
                $debug && print STDERR "$subn: \$key=$key \$i=$i \$cds_set=\n".Dumper($cds_set)."End of \$cds_set\n";
                $debug && print STDERR "$subn: \$key=$key \$i=$i \$inseqs=\n".Dumper($inseqs)."End of \$inseqs\n";
                foreach my $j (0 .. $#{$cds_set}) {
                  my $cds = $cds_set->[$j]->[0];
                  next if (!$cds);
                  $debug && print STDERR "$subn: \$j=$j \$cds=\n". Dumper($cds) . "End of \$cds\n";
                  my $acc = $cds->seq->accession_number;
#                  my $s3 = $cds->seq->translate->seq;
                  my $s3 = Annotate_Util::get_new_translation( $cds, $cds);
                  next if (!defined $s3 && print STDERR "$subn: got undef result. Skip \$j=$j acc=$acc\n");
                  
                  $s3 = $1 if ($s3 =~ /([^*]+)[*]$/); # to remove any trailing * in CDS
                  my @values = $cds->get_tag_values('translation');
                  my $s2 = $values[0];
                  $s2 = $1 if ($s2 =~ /([^*]+)[*]$/); # to remove any trailing * in CDS
                  $s3 =~ s/[*]/./g if ($s3 =~ /[*]/);
#                  $s3 =~ s/X/./ig if ($s3 =~ /X/i); # What's this for? This causes problem for FR675275 where this is a dot in the CDS sequence -11/19/2013
                  $s3 =~ s/[.]/X/ig; # -11/19/2013
                  $s2 =~ s/[.]/X/ig; # -11/19/2013
#                  if ($s2 ne $s3) {
                  if ($s2 !~ /$s3/ && Annotate_Verify::diff_2str( $s3, $s2) !~ /^L[.]+$/i) {
                    print STDERR "$subn: \$s2='$s2'\n";
                    print STDERR "$subn: \$s3='$s3'\n";
                    my $sdiff = Annotate_Verify::diff_2str( $s2, $s3);
                    print STDERR "$subn: diff='$sdiff'\n";
                    $sdiff =~ s/[.]//gi;
                    if (length($sdiff)>1) {
                      print STDERR "$subn: translation tag and translate don't match in CDS of $acc. Skip.\n";
                      print STDOUT "$subn: translation tag and translate don't match in CDS of $acc. Skip.\n";
                      next;
                    } else {
                      print STDERR "$subn: diff='$sdiff'\n";
                      print STDERR "$subn: The mismatch in CDS of $acc is only of 1 AA, proceed.\n";
                    }
                  }
                  $s2 =~ s/[BJZ]/X/ig; # so that clustalw to align BJZ, instead of printing '.', eg DQ835769
                  $debug && print STDERR "$subn: \$s2=$s2\n";
                  my $str2 = Annotate_Align::msa_seqid('ACC=', $cds);
                  $debug && print STDERR "$subn: \$str2=$str2\n";
                  my $f2 = Bio::PrimarySeq->new(
                         -seq      => $s2,
                         -id       => $str2,	# id can't contain space
                         -alphabet => 'protein'
                                    );
                  $cds_all = []; # Clear all cds, run MUSCLE separately for each input CDS
                  push @$cds_all, $f1;
                  if ( 1 ) {
                      Annotate_Def::add_extra_refCDS( $cds_all, $exe_dir);
                  }
                  push @$cds_all, $f2;

if ( 0 ) {
                $debug && print STDERR "$subn: \$cds_all=\n". Dumper($cds_all) . "End of \$cds_all\n";
} else {
                for my $cdsi (@$cds_all) { $debug && print STDERR "$subn: \$cdsi=". $cdsi->display_id . "\n"; }
}

        # Run CLUSTALW (or MUSCLE) for each CDS in refseq
        # Returns a SimpleAlign object
                printf STDERR "$subn: \$cds_all has %d genomes:", $#{$cds_all}+1;
                for my $iii ( 0 .. $#{$cds_all}) {
                    print STDERR " $iii:".$cds_all->[$iii]->display_id;
                }
                print STDERR "\n";
                if ($#{$cds_all}<1) {
                    print STDERR "$subn: too few seqs in \$cds_all $#{$cds_all}+1. Skip \$n_set=$n_set\n";
                    print STDERR "$subn: \$cds_all=\n". Dumper($cds_all) . "End of \$cds_all\n";
                    next;
                }
                #  Returns : Reference to a SimpleAlign object containing the sequence alignment
                my $factory;
                my $outfile_ext = $runMUSCLE ? 'afa' : 'msf';
                my $outfile_name = sprintf("$dir_path/test_n%d_i%d.$outfile_ext", $n_set, $j);
                print STDERR "$subn: MSA \$n_set=$n_set \$i=$i \$outfile_name='$outfile_name'\n";
                if ($debug) {
                    my $display_id = $cds_all->[1]->display_id;
                    my $display_id2 = '';
                   if ($display_id =~ /ACC=(.+)[|]/) {
                        $display_id2 = $1;
                    }
                    if ($display_id =~ /\|CDS=\D*(\d+[.]).*([.][<>]*\d+)\D*$/) {
                        $display_id2 .= "_${1}${2}";
                    }
#                    $display_id =~ s/\|CDS=/_/i;
                    $display_id2 =~ s/[.]{2}/_/i;
                    $display_id2 =~ s/[<>]//i;
                    $outfile_name = sprintf("$dir_path/test_%s_%02d.$outfile_ext", $display_id2, ++$uniqId);
                    $debug && print STDERR "$subn: MSA \$outfile_name='$outfile_name'\n";
                    if (-e "./$outfile_name") {
                        my $copies = Annotate_Util::backupFiles( '.', $outfile_name, 8);
                        print STDERR "$subn: sub backupFiles made $copies backups\n";
                    } else {
                        print STDERR "$subn: file:'./$outfile_name' doesn't exist, sub backupFiles skipped\n";
                    }
#                    my $count = 1;
#                    while (-f $outfile_name && $count <10) {
#                        $debug && print STDERR "$subn: \$display_id2=$display_id2 \$count=$count file:$outfile_name exists\n";
#                        $outfile_name = sprintf("$dir_path/test_%s_%02d.$outfile_ext", $display_id2, $count);
#                        $count++;
#                    }
                }
                $debug && print STDERR "$subn: ready to run alignment\n";
    if ($runMUSCLE) {
                @param = (
                   '-stable' => '',
#                   '-quiet' => '',
                   '-outfile_name' => "$outfile_name",
                         );
                $factory = Bio::Tools::Run::Alignment::Muscle->new(@param);
                $aln = $factory->align(
                         [@$cds_all]
                         );
                $debug && print STDERR "$subn: `cat $outfile_name`\n". `cat $outfile_name` . "End of $outfile_name\n";
#                $debug && print STDERR "$subn: \$aln=\n". Dumper($aln) . "End of \$aln\n\n";
    } else {
                @param = (
#                   '-quiet' => '',
                   '-OUTFILE' => "$outfile_name",
#                   '-OUTput' => "fasta",
                         );
                $factory = Bio::Tools::Run::Alignment::Clustalw->new(@param);
                $aln = $factory->align(
                         [@$cds_all]
                         );
#                $debug && print STDERR "$subn: \$aln=\n". Dumper($aln) . "End of \$aln\n\n";
                $debug && print STDERR "$subn: \$aln->gap_char='". $aln->gap_char . "'\n";
                $aln->gap_char('.'); # gap_char='.' for clustalw, '-' for MUSCLE
                $debug && print STDERR "$subn: \$aln->gap_char='". $aln->gap_char . "'\n";
    }
                # Print out the resulting MSA file
                $debug && print STDERR "$subn: `cat $outfile_name`\n". `cat $outfile_name` . "End of $outfile_name\n";

      next if (!$aln);
      push @$alns, $aln;

              $feats = Annotate_Align::MSA_annotate( $refseqs, $inseqs, $n_set, $aln, $exe_dir);
#              print STDERR "$subn: \$n_set=$n_set MSA returned \$feats=$#{$feats}\n";
              $debug && print STDERR "$subn: \$n_set=$n_set \$feats=\n". Dumper($feats) . "End of \$feats\n";

#              push @$feats_all, $feats;
              # reverse the indexes
              # Incoorporate the new features to any existing list of features
                  combineFeatures( $feats, $feats_all);

              $debug && print STDERR "$subn: \$n_set=$n_set \$feats_all=\n". Dumper($feats_all) . "End of \$feats_all\n";
                } # foreach my $j (0 .. $#{$cds_set})
            } # for (my $i = 0; $i<=$#{$inseqs}; $i++) {

        } # foreach my $key (keys %$refseqs) {

    } # for my $n_set (0 .. $n_sets)

    return $feats_all;

} # sub run_MSA

=head2 combineFeatures
Takes the result from MSA_annotate, check if the existing result already has any new features
special cases are:
AJ299464, where duplicate mat_peptides resulted from alternative reference mat_peptides

=cut

sub combineFeatures {
    my ($feats, $feats_all) = @_;

    my $debug = 0 || $debug_all;
    my $subn = 'combineFeatures';

    if ( 0 ) {
        $debug && print STDERR "$subn: Entering $subn \$feats_all=\n". Dumper($feats_all) . "End of \$feats_all\n";
    } else {
                $debug && print STDERR "$subn: Entering $subn \$feats_all=\n";
                for my $k1 (keys %$feats_all) {
                    for my $k2 (keys %{$feats_all->{$k1}}) {
                      for my $f (@{$feats_all->{$k1}->{$k2}}) {
                        $debug && print STDERR "$subn: feats_all: k1=$k1 k2=$k2 ".$f->primary_tag.':'.$f->location->to_FTstring."\n";
                      }
                    }
                }
    }
    if ( 0 ) {
        $debug && print STDERR "$subn: Entering $subn \$feats=\n". Dumper($feats) . "End of \$feats\n";
    } else {
                $debug && print STDERR "$subn: Entering $subn \$feats=\n";
                for my $k1 (keys %$feats) {
                    for my $k2 (keys %{$feats->{$k1}}) {
                      for my $f (@{$feats->{$k1}->{$k2}}) {
                        $debug && print STDERR "$subn: feats: k1=$k1 k2=$k2 ".$f->primary_tag.':'.$f->location->to_FTstring."\n";
                      }
                    }
                }
    }

    # Incoorporate the new features to any existing list of features
    for my $nk1 (sort keys %$feats) { # refcds GI
        $debug && print STDERR "$subn: \$nk1=$nk1: refcds\n";
        for my $nk2 (keys %{$feats->{$nk1}}) { # accession
            $debug && print STDERR "$subn: \$nk1=$nk1 \$nk2=$nk2: input genome\n";
#              if (!exists($feats_all->{$nk2}->{$nk1})) {
#                $debug && print STDERR "$subn: \$nk1=$nk1 \$nk2=$nk2 doesn't exist in \$feats_all, simply add\n";
#                $feats_all->{$nk2}->{$nk1} = $feats->{$nk1}->{$nk2};
#                next;
#              }
            for (my $kk = 1; $kk<=$#{$feats->{$nk1}->{$nk2}}; $kk++) {
                my $feat_kk = $feats->{$nk1}->{$nk2}->[$kk];
                my $prod_kk = [ $feat_kk->get_tag_values('product') ];
                $prod_kk = $prod_kk->[0];
                $debug && print STDERR "$subn: \$nk1=$nk1 \$nk2=$nk2 \$kk=$kk:$prod_kk\n";
                my $prod_kk_loc = $feat_kk->location->to_FTstring;
                my $prod_kk_note = [ $feat_kk->get_tag_values('note') ];
                for (@$prod_kk_note) { next if ($_!~/^Desc:/i); $prod_kk_note = $_; }
                $debug && print STDERR "\n$subn: \$nk1=$nk1 \$nk2=$nk2 \$kk=$kk \$prod_kk='".$feat_kk->primary_tag.":$prod_kk_loc:$prod_kk' is being checked\n";
                my $seen = 0;
#                my $feat_ff;
                # Check if the product name has been seen
                my @feats_all_keys = keys(%$feats_all);
                $debug && print STDERR "$subn: \$nk1=$nk1 \$nk2=$nk2 \$feats_all has ".scalar @feats_all_keys." keys: '@feats_all_keys'\n";
                for my $ok1 (keys %{$feats_all->{$nk2}}) {
                $debug && print STDERR "$subn: \$nk1=$nk1 \$nk2=$nk2 \$kk=$kk \$ok1=$ok1: input genome\n";
                    $debug && print STDERR "$subn: \$nk1=$nk1 \$nk2=$nk2 \$ok1=$ok1\n";
                    $debug && print STDERR "$subn: \$kk=$kk \$prod_kk='$prod_kk' \$prod_kk_loc='$prod_kk_loc'\n";
                    for my $ff (0 .. $#{$feats_all->{$nk2}->{$ok1}}) {
                        my $feat_ff = $feats_all->{$nk2}->{$ok1}->[$ff];
                        my $prod_f = [ $feat_ff->get_tag_values('product') ];
                        $prod_f = $prod_f->[0];
                        my $prod_f_loc = $feat_ff->location->to_FTstring;
                        my $prod_f_note = [ $feat_ff->get_tag_values('note') ];
                        for (@$prod_f_note) { next if ($_!~/^Desc:/i); $prod_f_note = $_; }

#                        next if ($prod_f ne $prod_kk);
                        if (!($prod_f_loc eq $prod_kk_loc || $prod_f eq $prod_kk)) {
                            $debug && print STDERR "$subn: \$ff=$ff $prod_kk_loc:$prod_kk vs $prod_f_loc:$prod_f are different\n";
                            next;
                        }
                        $debug && print STDERR "$subn: \$ff=$ff \$prod_f='$prod_f' $prod_f_loc\n";

                        $debug && print STDERR "$subn: \$prod_kk_loc=$prod_kk_loc:$prod_kk \n";
                        $debug && print STDERR "$subn: \$prod_f_loc =$prod_f_loc:$prod_f have similar location/product\n";
                        $debug && print STDERR "$subn: \$kk=$kk \$feat_kk=\n". Dumper($feat_kk) . "End of \$feat_kk\n";
                        $debug && print STDERR "$subn: \$ff=$ff \$feat_ff=\n". Dumper($feat_ff) . "End of \$feat_ff\n";
                        my $len_kk = $feat_kk->location->end   - $feat_kk->location->start +1;
                        my $len_ff = 0;
                        # $len_ff = $feat_ff->location->end   - $feat_kk->location->start +1;
                        $len_ff = $feat_ff->location->end   - $feat_ff->location->start +1;
                        $debug && print STDERR "$subn: \$len_ff=$prod_f_loc:$len_ff \$len_kk=$prod_kk_loc:$len_kk\n";

                        if ($prod_f eq $prod_kk && $prod_f_loc eq $prod_kk_loc) {
                          $seen = 1;
                          print STDERR "$subn: found duplicate by product and location: $prod_f_loc vs. $prod_kk_loc \$seen=$seen\n";
                          $debug && print STDERR "$subn: duplicate: \$prod_kk_note='$prod_kk_note'\n";
                          $debug && print STDERR "$subn: duplicate: \$prod_f_note ='$prod_f_note'\n";
                          $debug && print STDERR "$subn: duplicate: Both product and location are same, skip\n";
                          last;
                        } elsif ($prod_f eq $prod_kk) {
#                          my $len_kk = $feat_kk->location->end   - $feat_kk->location->start +1;
#                          my $len_ff = 0;
#                          # $len_ff = $feat_ff->location->end   - $feat_kk->location->start +1;
#                          $len_ff = $feat_ff->location->end   - $feat_ff->location->start +1;
                          $seen = 1;
                          print STDERR "$subn: found duplicate by product: $prod_f:$prod_f_loc vs. $prod_kk:$prod_kk_loc \$seen=$seen\n";
                          print STDERR "$subn: duplicate: \$prod_kk_note='$prod_kk_note'\n";
                          print STDERR "$subn: duplicate: \$prod_f_note ='$prod_f_note'\n";
                          my $rm_kk = '';
                          my $rm_f = '';
                          $rm_kk = $1 if ($prod_kk_note =~ m/[|](RM=[^|]+)[|]/i);
                          $rm_f = $1  if ($prod_f_note =~ m/[|](RM=[^|]+)[|]/i);
                          $debug && print STDERR "$subn: $prod_f:$prod_f_loc:$rm_f vs. $prod_kk:$prod_kk_loc:$rm_kk\n";
                          if ( ($rm_kk || $rm_f) && $rm_kk ne $rm_f ) {
                            print STDERR "$subn: different reference mat_peptide: $rm_kk vs. $rm_f, keep both.\n";
                            $seen = 0;
                            last;
                          } elsif ( $len_kk > $len_ff ) {
                            print STDERR "$subn: replacing $prod_f:$prod_f_loc with $prod_kk:$prod_kk_loc \$seen=$seen\n";
                            $feats_all->{$nk2}->{$ok1}->[$ff] = $feat_kk;
                            next;
                          }
                        } elsif ($prod_f_loc eq $prod_kk_loc) {
                          my $pct = 0;
                          $pct = Annotate_Util::cmp_cds_bl2seq( $feat_kk, $feat_ff);
                          $debug && print STDERR "$subn: ".$prod_f_loc." with $prod_kk_loc \$pct=$pct\n";

                          if ( $pct > 0.99 ) {
                            $seen = 1;
                            print STDERR "$subn: found duplicate by location: $prod_f_loc vs. $prod_kk_loc \$seen=$seen\n";
                            print STDERR "$subn: duplicate: \$prod_kk_note='$prod_kk_note'\n";
                            print STDERR "$subn: duplicate:  \$prod_f_note='$prod_f_note'\n";
                            if ( $prod_f_note =~ /[|]Partial=Y[|]/i ) { # take new feature if there is no 'Partial=Y'
                              print STDERR "$subn: replacing $prod_f:$prod_f_loc with $prod_kk:$prod_kk_loc \$seen=$seen\n";
                              $feats_all->{$nk2}->{$ok1}->[$ff] = $feat_kk;
                            }
                          }
                          last;
#                        } else {
#                          next;
                        }
                    } # for my $ff (0 .. $#{$feats_all->{$nk2}->{$ok1}})
                    last if ($seen);
                } # for my $ok1 (keys %{$feats_all->{$nk2}})
                $debug && print STDERR "$subn: after checking \$feats_all \$nk1=$nk1 \$nk2=$nk2 ".$feat_kk->primary_tag.":$prod_kk_loc \$seen=$seen\n";
                if ($seen) {
                    $debug && print STDERR "$subn: after checking \$feats_all \$nk1=$nk1 \$nk2=$nk2 ".$feat_kk->primary_tag.":$prod_kk_loc skipped\n";
                    next;
                }
                $debug && print STDERR "$subn: Add to list of features: ".$feat_kk->primary_tag.":$prod_kk_loc \$seen=$seen, keep\n";
                push @{$feats_all->{$nk2}->{$nk1}}, $feat_kk;
        if ( 0 ) {
                $debug && print STDERR "$subn: Leaving $subn \$feats_all=\n". Dumper($feats_all) . "End of \$feats_all\n";
        } else {
                for my $k1 (keys %$feats_all) {
                    for my $k2 (keys %{$feats_all->{$k1}}) {
                      for my $f (@{$feats_all->{$k1}->{$k2}}) {
                        $debug && print STDERR "$subn: feats_all: k1=$k1 k2=$k2 ".$f->primary_tag.':'.$f->location->to_FTstring."\n";
                      }
                    }
                }
        }
            }
        } # for my $nk2 (keys %{$feats->{$nk1}}) {
    } # for my $nk1 (keys %$feats) {
    $debug && print STDERR "$subn: Done with $subn: \$feats_all=\n". Dumper($feats_all) . "End of \$feats_all\n";
    $debug && print STDERR "$subn: Done with $subn\n";

    return $feats_all;

} # combineFeatures


=head2 MSA_annotate
Takes a hash of refseqs, in the form of { accession => [ [CDS, mat_peptides] ... ], []... }
 and an array of target genomes, in the form of [ accession, [ [CDS, mat_peptides, ...], []... ] ]
 and alignment from MUSCLE
 Return an array of alignments CDS1, CDS2, CDS3 ...
=cut

sub MSA_annotate {
    my ($refseqs, $inseqs, $n_set, $aln,$exe_dir) = @_;

    my $debug = 0 || $debug_all;
    my $subn = 'MSA_annotate';

    my $refset = [];    # [CDS, mat_peptides] ... ]

    $debug && print STDERR "$subn: \$aln=\n". Dumper($aln) . "End of \$aln\n";
    # get $refset
    foreach my $key (sort keys %$refseqs) {
        $debug && print STDERR "$subn: \$key=$key\n";
        $refset = $refseqs->{$key}->[$n_set];
#        $debug && print STDERR "$subn: \$refset=\n". Dumper($refset) . "End of \$refset\n";
#        my $refcds = $refseqs->{$key}->[$n_set]->[0];
#        $debug && print STDERR "$subn: \$refcds=". $refcds->seq->accession_number ."\n";
#        $debug && print STDERR "$subn: \$refcds=\n". Dumper($refcds) ."End of \$refcds\n";
        last;
    }

    # get $inset
    my $feats_all = {};
    my $inset = [];    # [CDS, mat_peptides] ... ]

    for (my $i = 0; $i<=$#{$inseqs}; $i++) {
      my $inseq = $inseqs->[$i];
      for my $key (sort keys %$refseqs) {
        my $refcds_id = $refseqs->{$key}->[$n_set]->[0];
        $debug && print STDERR "$subn: \$refcds_id=$refcds_id\n";
        next if (!exists($inseq->[1]->{$refcds_id}));

        $refset = $refseqs->{$key}->[$n_set];
if (0) {
        $debug && print STDERR "$subn: \$key=$key \$refset=\n". Dumper($refset) . "End of \$refset\n";
} else {
        for my $k (0.. $#{$refset}) {
            my $fk = $refset->[$k];
            next if (!$fk->isa('Bio::SeqFeature::Generic'));
            $debug && print STDERR "$subn: \$k=$k \$refseq=".$fk->primary_tag.':'.$fk->location->to_FTstring."\n";
        }
}

#        my $cds_sets = $inseq->[$n_set+1]->{$refcds_id};
        my $cds_sets = $inseq->[1]->{$refcds_id};
        $debug && print STDERR "$subn: \$i=$i \$cds_sets=\n". Dumper($cds_sets) . "End of \$cds_sets\n";
        foreach $inset (@$cds_sets) {
#          $inset  = $inseqs->[$i]->[$n_set+1]->[0];
          $debug && print STDERR "$subn: \$i=$i \$inset=\n". Dumper($inset) . "End of \$inset\n";

          # find the new annotations
          my $feats_new = Annotate_Align::MSA_annotate_1cds( $refset, $inset, $aln,$exe_dir);
          $debug && print STDERR "$subn: MSA_annotate_1cds returns \$feats_new=\n".Dumper($feats_new)."End of \$feats_new\n";
          $debug && print STDERR "$subn: After MSA_annotate_1cds \$feats_all=\n".Dumper($feats_all)."End of \$feats_all\n";

#          push @$feats_all, $feats_new;
          my $acc = $inset->[0]->seq->accession_number;
          $debug && print STDERR "$subn: \$refcds_id=$refcds_id \$acc=$acc\n";
          if (!exists($feats_all->{$refcds_id}) && !exists($feats_all->{$refcds_id}->{$acc})) {
              $debug && print STDERR "$subn: \$refcds_id=$refcds_id found new mat_peptide=". $#{$feats_new} . "\n";
              $feats_all->{$refcds_id}->{$acc} = $feats_new;
          } else {
              for (my $kk = 0; $kk<=$#{$feats_new}; $kk++) {
                  # Skip the mat_peptide if it's already in the list
#                  my $seen = 0;
#                  $debug && print STDERR "$subn: \$refcds_id=$refcds_id ". $feats_new->[$kk]->location->to_FTstring . " is being checked\n";
#                  for my $ff (0 .. $#{$feats_all->{$refcds_id}->{$acc}}) {
#                    my $f = $feats_all->{$refcds_id}->{$acc}->[$ff];
#                    $debug && print STDERR "$subn: \$ff=$ff ". $f->location->to_FTstring . " is being checked\n";
#                    next if ($f->location->to_FTstring ne $feats_new->[$kk]->location->to_FTstring);
#                    $debug && print STDERR "$subn: \$kk=$kk \$f=". $f->location->to_FTstring . " already seen, skip\n";
#                    $seen = 1;
#                    last;
#                  }
#                  $debug && print STDERR "$subn: \$kk=$kk ". $feats_new->[$kk]->location->to_FTstring . " \$seen=$seen, keep\n";
#                  next if ($seen);
#                  $debug && print STDERR "$subn: \$kk=$kk ". $feats_new->[$kk]->location->to_FTstring . " is new, keep\n";
                  push @{$feats_all->{$refcds_id}->{$acc}}, $feats_new->[$kk];
              }
          }

          $debug && print STDERR "$subn: After \$inset \$feats_all=\n". Dumper($feats_all) . "End of \$feats_all\n";
        }
      }
    }
    $debug && print STDERR "$subn: Finished with MSA_annotate \$feats_all=\n". Dumper($feats_all) . "End of \$feats_all\n";

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
    my $subn = 'MSA_annotate_1cds';
    my $feats_all = [];

    my $refcds = $refset->[1];
    $debug && print STDERR "$subn: \$refcds='".$refcds->seq->accession_number."'\n";
#    $debug && print STDERR "MSA_annotate_1cds: \$refcds=\n". Dumper($refcds) . "End of \$refcds\n";

    my @values = $refcds->get_tag_values('translation');
    my $s1 = $values[0];
    $s1 = $1 if ($s1 =~ /([^*]+)[*]$/); # to remove any trailing * in CDS
    my $refcds_id = Annotate_Align::msa_seqid('ref=', $refcds);
    my $aln_h = Annotate_Util::msa_get_aln( $aln, $refcds_id);
    if (!$aln_h) {
        $refcds_id =~ s/[(]/_/g;
        $refcds_id =~ s/[)]/_/g;
        $refcds_id =~ s/[,]/_/g;
        $aln_h = Annotate_Util::msa_get_aln( $aln, $refcds_id);
    }
    if (!$aln_h) {
        print STDERR "$subn: null alignment for \$refcds_id='$refcds_id', quit";
        exit;
    }

    my $cds = $inset->[0];
    next if (!$cds);
    my $acc = $cds->seq->accession_number;
#    $debug && print STDERR "$subn: \$cds=\n". Dumper($cds) . "End of \$cds\n";
#    $debug && print STDERR "$subn: \$inset=\n". Dumper($inset) . "End of \$inset\n";

    my $s3 = Annotate_Util::get_new_translation( $cds, $cds);
    if (!defined $s3) {
       print STDERR "$subn: sub get_new_translation returned undef result. Skip accession=$acc\n";
       print STDOUT "$subn: sub get_new_translation returned undef result. Skip accession=$acc\n";
       return;
    }
    $s3 = $1 if ($s3 =~ /([^*]+)[*]$/); # to remove any trailing * in CDS
    @values = $cds->get_tag_values('translation');
    my $s2 = $values[0];
    $s2 = $1 if ($s2 =~ /([^*]+)[*]$/); # to remove any trailing * in CDS
    $s3 =~ s/[*]/./g if ($s3 =~ /[*]/);
    $s3 =~ s/[X]/./ig if ($s3 =~ /[X]/i);
#    if ($s2 ne $s3) {
    if ($s2 !~ /$s3/ && Annotate_Verify::diff_2str( $s3, $s2) !~ /^L[.]+$/i) {
        print STDERR "$subn: \$s2='$s2'\n";
        print STDERR "$subn: \$s3='$s3'\n";
        print STDERR "$subn: diff='".Annotate_Verify::diff_2str( $s2, $s3)."'\n";
        print STDERR "$subn: translation tag and translate don't match in CDS of '$acc'.\n";
        print STDERR "$subn: Skip.\n\n";
        return $feats_all;
#        next;
    }
    my $cds_id = Annotate_Align::msa_seqid('ACC=', $cds);
    $debug && print STDERR "$subn: \$refcds_id='$refcds_id' \$cds_id='$cds_id'\n";

    my $note = "Annotated by VIPRBRC, MSA, refseq=".$refcds->seq->accession_number();
    $aln_h = Annotate_Util::msa_get_aln( $aln, $cds_id);
    if (!$aln_h) {
        $cds_id =~ s/[(]/_/g;
        $cds_id =~ s/[)]/_/g;
        $cds_id =~ s/[,]/_/g;
        $aln_h = Annotate_Util::msa_get_aln( $aln, $cds_id);
    }
    $debug && print STDERR "$subn: \$aln_h=\n". Dumper($aln_h) . "End of \$aln_h\n";
#    if ( !$aln_h || !$aln_h->isa('Bio::LocatableSeq')) {
    if (!$aln_h) {
        $debug && print STDERR "$subn: \$aln_h is empty\n";
        return $feats_all; # This is required --10/16/2012
    }
    $feats_all = Annotate_Util::project_matpept(
                         $refset,
                         $inset,
                         $aln,
                         $refcds_id,
                         $cds_id,
                         $note,$exe_dir,
                         );

    $debug && print STDERR "$subn: \$feats_all=$#{$feats_all}\n". Dumper($feats_all) . "End of \$feats_all\n";

    return $feats_all;
} # sub MSA_annotate_1cds


1;
