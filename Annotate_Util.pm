package Annotate_Util;

use strict;
use warnings;
use English;
use Carp;
use Data::Dumper;

use File::Temp qw/ tempfile tempdir /;
use Bio::SeqIO;
use Bio::Seq;
use Bio::AlignIO;
use Bio::Tools::Run::StandAloneBlast;
use IO::String;

use Annotate_Math;
use Annotate_Verify;
use Annotate_Def;
use version; our $VERSION = qv('1.1.7'); # Feb 05 2013

my $debug_all = 0;

my $MIN_CDS_LENGTH = 30; # any CDS or mat_peptide shorter than this will be discarded
#$MIN_CDS_LENGTH = 24; # set to 24, for genomes such as AB363975, whose CDS=48

####//README//####
#
# Annotate_Util contains the core functions to perform annotation based on both MUSCLE and bl2seq alignment
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


=head2
  sub getPerlVersions finds out the version of the important BioPerl modules used
/perl -MBio::Root::Version -e 'print $Bio::Root::Version::VERSION,"\n"'/
=cut
sub getPerlVersions{
    my $debug = 0 || $debug_all;
    my $subname = "Annotate_Util::getPerlVersions";
    my $msgs = [];
    my $cmds = [ ['perl', ''],
                 ['Bio::Root::Version', '$Bio::Root::Version::VERSION'],
                 ['Bio::Seq',   '$Bio::Seq::VERSION'],
                 ['Bio::SeqIO',   '$Bio::SeqIO::VERSION'],
                 ['Bio::AlignIO', '$Bio::AlignIO::VERSION'],
#                 ['Bio::TreeIO',  '$Bio::TreeIO::VERSION'],
#                 ['Bio::Tools::Run',  'Bio::Tools::Run::VERSION'],
                 ['Bio::Tools::Run::StandAloneBlast',  'Bio::Tools::Run::StandAloneBlast::VERSION'],
               ];
    for my $i (0 .. $#{$cmds}) {
        my $cmd;
        if ($cmds->[$i]->[0] eq 'perl') {
            $cmd = 'perl -v | head -n3';
        } else {
            $cmd = "perl -M$cmds->[$i]->[0] -e 'print $cmds->[$i]->[1]'"
        }
        $cmd = ` $cmd `;
        $cmd =~ s/(^\n|\n$)//ig; # Remove the leading return
        push @$msgs, sprintf("%-14s'", "$cmds->[$i]->[0]="). "$cmd'";
    }

    for (@$msgs) {
        $debug && print STDERR "$subname: $_\n";
    }
    return $msgs;
} # sub getPerlVersion


=head2
my $TAXON = {
               taxon_loaded => 0,
               taxon_fn     => "Annotate_taxon_records.txt",
            };
my $REFSEQS = {
            };
my $gene_symbol2 = {
               symbol_loaded => 0,
               symbol_fn     => "Annotate_symbol_records.txt",
            };

=cut

=head2
sub backupFiles takes a directory path, a filename, and a number, backs up the file to <filename>.bakx with x=1..number
=cut

sub backupFiles {
    my ($exe_dir, $filename, $numBak) = @_;

    my $debug = 0 || $debug_all;
    my $subn = 'Annotate_Util::backupFiles';

    my $numCopies = 0;
    $numBak = 1 if (!$numBak || $numBak<1);
    $exe_dir = './' if (!$exe_dir);
    $debug && print STDERR "$subn: \$exe_dir='$exe_dir' \$numBak=$numBak\n";
    return $numCopies if (!$filename);
    if (!-e "$exe_dir/$filename") {
        print STDERR "$subn: file='$exe_dir/$filename' doesn't exist, nothing to backup\n";
        return $numCopies;
    }

    # First, see if file $filename exists, then see if 1st backup exists, then further backups
    my $result = '';
    my $bak = "$filename.bak1";
    my $cmds = [ ];
    if (!-e "$exe_dir/$bak") {
        $debug && print STDERR "$subn: file='$exe_dir/$bak' doesn't exist, backup one copy\n";
        $result = ` ls -l $exe_dir/$bak 2>&1`;
    } else {
        $result = "diff $exe_dir/$filename $exe_dir/$bak 2>&1";
        $debug && print STDERR "$subn: command='$result'\n";
        $result = ` $result `;
        $debug && print STDERR "$subn: \$result='\n$result'\n";
        chomp $result;
    }
    if (!$result) {
        print STDERR "$subn: no difference between file $filename vs. $bak\n";
    } else {
        $cmds->[$#{$cmds}+1] = "mv $exe_dir/$filename $exe_dir/$bak 2>&1";
    }
    $debug && print STDERR "$subn: \$cmds=".Dumper($cmds)."\n";
    for my $i (1 .. $numBak) {
            $bak = sprintf("${filename}.bak%d", $i);
            $debug && print STDERR "$subn: \$i=$i \$bak=$bak\n";
            if (!-e"$exe_dir/$bak") {
                $debug && print STDERR "$subn: \$i=$i \$bak=$bak doesn't exist.\n";
                last;
            }
            my $bak2 = sprintf("${filename}.bak%d", $i+1);
#            if (!-e"$exe_dir/$bak2") {
#                $debug && print STDERR "$subn: \$i=$i \$bak=$bak2 doesn't exist.\n";
#                last;
#            }
            my $cmd = "diff $exe_dir/$bak $exe_dir/$bak2 2>&1";
            $result = ` $cmd `;
            chomp $result;
            $debug && print STDERR "$subn: \$cmd='$cmd' \$result='\n$result'\n";
            if (!$result) {
                $debug && print STDERR "$subn: \$i=$i $exe_dir/$bak $exe_dir/$bak2 have no difference.\n";
                last;
            }
            $cmds->[$#{$cmds}+1] = "mv $exe_dir/$bak $exe_dir/$bak2 2>&1";
    }
    $debug && print STDERR "$subn: \$cmds=".Dumper($cmds)."\n";

    $result = '';
    for my $i (reverse 0 .. $#{$cmds}) {
        my $cmd = $cmds->[$i];
        if ($cmd) {
            $numCopies++;
            $debug && print STDERR "$subn: \$i=$i \$numCopies=$numCopies \$cmd='$cmd'\n";
            $result = ` $cmd `;
            $debug && print STDERR "$subn: \$i=$i \$numCopies=$numCopies \$result='\n$result'\n";
            chomp $result;
            $debug && print STDERR "$subn: \$cmd=$cmd\n";
            print STDERR "$subn: \$result='\n$result'\n" if ($result);
        }
    }
    $debug && print STDERR "$subn: \$numCopies=$numCopies \$result='$result'\n";

    $result = "ls -l $exe_dir/$filename* 2>&1";
    $result .= "\n" . `ls -l $exe_dir/$filename* 2>&1`;
    $debug && print STDERR "$subn: \$result='\n$result'\n";
    $debug && print STDERR "$subn: backed up $numCopies copies for file:$filename\n";
    return $numCopies;
} # sub backupFiles


=head2 msa_get_aln
Takes an alignment, and an id, return the gaps within the alignment with such id
=cut

sub msa_get_aln {
    my ($aln, $id) = @_;

    my $debug = 0 ;#|| $debug_all;
    my $subn = 'msa_get_aln';
    my $refaln = [ $aln->each_seq_with_id($id) ]; # 
#    $debug && print STDERR "$subn: \$refaln=\n".Dumper($refaln)."End of \$refaln\n\n";
    if ($#{$refaln} < 0) {
        $debug && print STDERR "$subn: Couldn't find id=$id in alignment file.\n";
        return undef;
    }
    my $seq = $refaln->[0];
    if (!$seq) {
        $id =~ s/[(]/_/g;
        $id =~ s/[)]/_/g;
        $id =~ s/[,]/_/g;
        $refaln = [ $aln->each_seq_with_id($id) ]; # 
        $seq = $refaln->[0];
    }
    $debug && print STDERR "$subn: \$refaln = ".$seq->display_id()."\n";
    $debug && print STDERR "$subn: \$refaln = ".$seq->start() ."..". $seq->end()."\n";
    $debug && print STDERR "$subn: \$refaln = ".$seq->alphabet() ."..". $seq->length."\n";
#    $debug && print STDERR "$subn: \$refaln = ".$seq->seq."\n";

    my @gaps = (''); # we don't need element 0, but need it to suppress error from print
    $seq = $refaln->[0];
    $debug && print STDERR "$subn: \$seq='$seq'\n";

    return ($seq);
} # sub msa_get_aln


=head2 run_1MSA
 Takes an array of CDS, runs either muscle or clustalw
 Returns a SimpleAlign object
=cut

sub run_1MSA {
    my ($cds_all, $runMUSCLE, $outfile_name, $dir_path, $progs) = @_;

    my $debug = 1 || $debug_all;
    my $subn = 'Annotate_Util::run_1MSA';

    # $runMUSCLE = 0; # 0: Muscle, 1: Clustalw
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

    my $aln;
    my @param;
    if ( 0 ) {
                $debug && print STDERR "$subn: \$cds_all=\n". Dumper($cds_all) . "End of \$cds_all\n";
    } else {
                for my $cdsi (@$cds_all) { $debug && print STDERR "$subn: \$cdsi=". $cdsi->display_name . "\n"; }
    }

    # Run CLUSTALW (or MUSCLE) for each CDS in refseq
    # Returns a SimpleAlign object
    printf STDERR "$subn: \$cds_all has %d genomes:", $#{$cds_all}+1;
    for my $iii ( 0 .. $#{$cds_all}) {
        print STDERR " $iii:".$cds_all->[$iii]->display_name;
    }
    print STDERR "\n";
    if ($#{$cds_all}<1) {
        print STDERR "$subn: too few seqs in \$cds_all $#{$cds_all}+1. Skip\n";
        print STDERR "$subn: \$cds_all=\n". Dumper($cds_all) . "End of \$cds_all\n";
        next;
    }
    #  Returns : Reference to a SimpleAlign object containing the sequence alignment
    print STDERR "$subn: MSA \$outfile_name='$outfile_name'\n";
    $debug && print STDERR "$subn: ready to run alignment\n";
    my $factory;
    if ($runMUSCLE) {
        @param = (
           '-stable' => '',
#           '-quiet' => '',
                 );
        ($outfile_name) && push @param, ('-outfile_name' => "$outfile_name",);
        $factory = Bio::Tools::Run::Alignment::Muscle->new(@param);
        $aln = $factory->align(
                 [@$cds_all]
                 );
        $debug && print STDERR "$subn: `cat $outfile_name`\n". `cat $outfile_name` . "End of $outfile_name\n";
#        $debug && print STDERR "$subn: \$aln=\n". Dumper($aln) . "End of \$aln\n\n";
    } else {
        @param = (
#           '-quiet' => '',
#           '-OUTput' => "fasta",
                 );
        ($outfile_name) && push @param, ('-OUTFILE' => "$outfile_name",);
        $factory = Bio::Tools::Run::Alignment::Clustalw->new(@param);
        $aln = $factory->align(
                 [@$cds_all]
                 );
#        $debug && print STDERR "$subn: \$aln=\n". Dumper($aln) . "End of \$aln\n\n";
#        $debug && print STDERR "$subn: \$aln->gap_char='". $aln->gap_char . "'\n";
        $aln->gap_char('.'); # gap_char='.' for clustalw, '-' for MUSCLE
        $debug && print STDERR "$subn: \$aln->gap_char='". $aln->gap_char . "'\n";
    }
    # Print out the resulting MSA file
    $debug && print STDERR "$subn: `pwd`=". `pwd `;
    $debug && print STDERR "$subn: `cat $outfile_name`\n". `cat $outfile_name` . "End of $outfile_name\n";

    return $aln;

} # sub run_1MSA


sub msa_seqid {
    my ($str1, $cds) = @_;

    my $debug = 0 && $debug_all;
    my $subn = 'msa_seqid';

    $str1 .= $cds->seq->accession_number if ($cds && $cds->seq);
    if ( 0 ) {
        $str1 .= '|'. $cds->primary_tag.'='.$cds->location->to_FTstring;
    } else {
        $str1 .= '|'. $cds->primary_tag.'='.$cds->location->start.'..'.$cds->location->end;
    }
    $debug && print STDERR "$subn: \$str1='$str1'\n";
    return $str1;
} # msa_seqid


sub cds2PrimarySeq {
    my ($cds) = @_;

    my $debug = 0 && $debug_all;
    my $subn = 'cds2PrimarySeq';

    my $f2 = undef;
    my $ecode = '';

    $debug && print STDERR "$subn: \$cds=\n". Dumper($cds) . "End of \$cds\n";
    
    my $acc = '';
    if ($cds->seq) {
        $acc = $cds->seq->accession_number
    } else {
        for my $note ($cds->get_tag_values('note')) {
            $acc = $1 if ($note =~ m/^Desc:.+[|]ACC=([^|]+)[|]/i);
            $debug && print STDERR "$subn: \$acc='$acc' \$note='$note'\n";
        }
    }

    my @values = $cds->get_tag_values('translation');
    my $s2 = $values[0];
    $s2 = $1 if ($s2 =~ /([^*]+)[*]$/); # to remove any trailing * in CDS
    $s2 =~ s/[.]/X/ig; # -11/19/2013
#    my $s3 = $cds->seq->translate->seq;
    my $s3 = Annotate_Util::get_new_translation( $cds, $cds);
    if (!defined $s3) {
        $ecode = "sub get_new_translation returned undef for acc=$acc";
        return ($f2, $ecode);
    } elsif($s3) {
      $s3 = $1 if ($s3 =~ /([^*]+)[*]$/); # to remove any trailing * in CDS
      $s3 =~ s/[*]/./g if ($s3 =~ /[*]/);
#      $s3 =~ s/X/./ig if ($s3 =~ /X/i); # What's this for? This causes problem for FR675275 where this is a dot in the CDS sequence -11/19/2013
      $s3 =~ s/[.]/X/ig; # -11/19/2013
#      if ($s2 ne $s3) {
      if ($s2 !~ /$s3/ && Annotate_Verify::diff_2str( $s3, $s2) !~ /^L[.]+$/i) {
        print STDERR "$subn: \$s2='$s2'\n";
        print STDERR "$subn: \$s3='$s3'\n";
        my $sdiff = Annotate_Verify::diff_2str( $s2, $s3);
        print STDERR "$subn: diff='$sdiff'\n";
        $sdiff =~ s/[.]//gi;
        if (length($sdiff)>1) {
            $ecode = "translation tag and translate don't match in CDS of $acc. Skip.";
            return ($f2, $ecode);
        } else {
            print STDERR "$subn: diff='$sdiff'\n";
            print STDERR "$subn: The mismatch in CDS of $acc is only of 1 AA, proceed.\n";
        }
      }
    }
    $s2 =~ s/[BJZ]/X/ig; # so that clustalw to align BJZ, instead of printing '.', eg DQ835769
    $debug && print STDERR "$subn: \$s2=$s2\n";
    my $str2 = Annotate_Util::msa_seqid('ACC=', $cds);
    $debug && print STDERR "$subn: \$str2=$str2\n";
    $f2 = Bio::PrimarySeq->new(
                         -seq      => $s2,
                         -id       => $str2,    # id can't contain space
                         -alphabet => 'protein'
                              );

    $debug && print STDERR "$subn: \$f2=\n".Dumper($f2)."End of \f2$\n";
    return ($f2, $ecode);
} # cds2PrimarySeq


sub refcds2PrimarySeq {
    my ($refcds) = @_;

    my $debug = 0 && $debug_all;
    my $subn = 'refcds2PrimarySeq';

    my @values = $refcds->get_tag_values('translation');
    my $s1 = $values[0];
    $s1 = $1 if ($s1 =~ /([^*]+)[*]$/); # to remove any trailing * in CDS
    my $str1 = Annotate_Util::msa_seqid('ref=', $refcds);
    $debug && print STDERR "$subn: \$str1=$str1\n";
    my $f1 = Bio::PrimarySeq->new(
                         -seq      => $s1,
                         -id       => $str1,    # id can't contain space
                         -alphabet => 'protein'
                                    );

    $debug && print STDERR "$subn: \$f1=\n".Dumper($f1)."End of \f1$\n";
    return $f1;
} # refcds2PrimarySeq


=head2 get_polyprots
 Takes an inseq object, and array of [CDS, mat_peptide 1, mat_peptide 2, ...]
 return the {CDS_id => [[CDS, mat_peptide 1, mat_peptide 2, ...] [] ...] in the sequence and refseq. Empty if there is any problem
=cut

sub get_polyprots {
    my ($inseq, $refpolyprots, $exe_dir, $dir_path, $progs) = @_;

    my $debug = 1 || $debug_all;
    my $subn = 'get_polyprots';
    my $polyprots = {};
    my $num_cds = -1;
    my $comment = 'No refpolyprots provided';
#    $debug && print STDERR "$subn: \$refpolyprots=\n".Dumper($refpolyprots)."end of \$refpolyprots\n\n";
    if ($debug) {
        for my $gi (0 .. $#{$refpolyprots}) {
            my $g = $refpolyprots->[$gi];
            for my $fi (0 .. $#{$g}) {
                my $f = $g->[$fi];
                print STDERR "$subn: gi=$gi fi=$fi:\t";
                if ($fi==0) {
                    print STDERR "$f\n";
                    next;
                }
                print STDERR $f->seq_id.'|'.$f->primary_tag.'='.$f->location->to_FTstring."\n";
            }
        }
#        print STDERR "$subn: \$refpolyprots=\n".Dumper($refpolyprots)."End of \$refpolyprots\n\n";
    }


    if ($#{$refpolyprots} < 0) {
            return ($polyprots, $num_cds, $comment);
    }
    $comment = 'No CDS defined in genome file';

    my $cds_total = 0;
    my $acc = $inseq->accession_number;
    my $found_polyprot = 0;
    my $matchHigh = 0;
    my $feats = [ $inseq->get_SeqFeatures ];

    for (my $ct = 0; $ct<=$#{$feats}; $ct++) {
        my $feat = $feats->[$ct];
        $debug && printf STDERR ("$subn: \$ct=$ct \$feat=%-8s", $feat->primary_tag);
        $debug && printf STDERR (" \$loc=%s", $feat->location->to_FTstring);
        my $len = $feat->location->end - $feat->location->start +1;
        if ($feat->primary_tag ne 'CDS') {
            $debug && print STDERR " is not CDS. Skip.\n";
            next; # Skip those features such as source, gene, 5'URT
        } elsif (!$feat->has_tag('translation')) {
          if ( 0 ) {
            $debug && print STDERR " has no translation. Skip.\n";
            next; # Skip
          } else {
             $feat->add_tag_value('translation', Annotate_Util::get_new_translation( $feat, $feat)); 
          }
        } elsif ($len < $MIN_CDS_LENGTH * 1.6) { # why do we have a factor of 2?
            $debug && print STDERR " lengh=$len, shorter than required $MIN_CDS_LENGTH*1.6. Skip.\n";
            next; # Skip those CDS that are too short
        } else {
            $debug && print STDERR " (potential polyprotein)\n";
        }
        $cds_total++;
        $comment = "Found $cds_total CDS, but doesnot align with refseq";

        my $cds = $feat;
        my $s_cds_loc = $cds->location->to_FTstring;
        my $s_cds_id = $acc .'|'. $s_cds_loc;
        # polyprotein is determined by bl2seq similarity search, not by label of polyprotein

        # double check to ensure the translation of CDS is right. For instance, DQ430819
        my $tr2 = &get_new_translation($cds, $cds);
        if (!$tr2) {
            print STDERR "$subn: get_new_translation returned translation='undef' for CDS=".$cds->location->to_FTstring.". Skip.\n";
            $comment = 'Problem with translation of CDS='.$cds->location->to_FTstring;
            next; # Skip if no translation can be obtained
        }
        $debug && print STDERR "$subn: translation='$tr2'\n";

        my $PCT_CONSERVED_MIN = 0.6667;
        $PCT_CONSERVED_MIN = 0.5667; # the earlier value of 2/3 would cause problem in FJ872155:CDS<1..>121
        $PCT_CONSERVED_MIN = 0.5 if ($debug);
        my $pct_conserved = 0;
        my $refcds_set;
        $matchHigh = 0; # reset for each feat in $inseq
        my $bestMatch = undef; # hold the best match so far
        $debug && print STDERR ("$subn: \$ct=$ct Now look through ".scalar($refpolyprots)." refcds for best cookie cutter\n");
        for (my $i = 0; $i<=$#{$refpolyprots}; $i++) {
            my $refset = $refpolyprots->[$i];
            my $reffeat = $refset->[1]; # [0] is CDS_id
            my $s_ref_loc = $reffeat->location->to_FTstring;
            my $s_ref_id = $reffeat->seq->accession_number.'|'.$s_ref_loc;
            $s_ref_id = $refset->[0] .'|'. $s_ref_id;
            $debug && print STDERR "$subn: \$i=$i \$refseq=$s_ref_id \$refcds=$s_ref_loc\n";
            next if ($reffeat->primary_tag ne 'CDS');

            my $match = 0.0;
if ( 0 ) {
            $match = Annotate_Util::cmp_cds_bl2seq($cds, $reffeat);
            $debug && printf STDERR ("$subn: CDS=$s_cds_id $s_ref_id \$match=%6.4f \$matchHigh=$matchHigh from sub cmp_cds_bl2seq\n", $match);
} else {
            my $runMUSCLE = 0;
            my $outfile_name = sprintf("$dir_path/test_ref%d_cds%d", $i, $ct);
            $outfile_name .= $runMUSCLE ? '.afa' : '.msf';
            $match = Annotate_Util::cmp_cds_clustalw($cds, $reffeat, $outfile_name, $runMUSCLE, $exe_dir, $dir_path, $progs);
#            $debug && printf STDERR ("$subn: CDS=$s_cds_id $s_ref_id \$match=%6.4f \$matchHigh=$matchHigh from sub cmp_cds_clustalw\n", $match);
}
#            $match = 100 if ($debug);
            $debug && printf STDERR ("$subn: CDS=$s_cds_id \$pct_conserved=%6.4f \$match=%6.4f \$refcds=$s_ref_id\n", $pct_conserved, $match);
            if ($match < 0.5) {
                $debug && printf STDERR ("$subn: \$match=%6.4f is too low for $s_ref_id to be considered a good refcds\n", $match);
                next;
            } else {
                $debug && printf STDERR ("$subn: \$match=%6.4f (>0.5) for $s_ref_id is a good refcds\n", $match);
            }
#            $matchHigh = $match if ($matchHigh < $match);
            if ($matchHigh < $match) {
                $matchHigh = $match;
            } elsif ($matchHigh == $match) { # NC_001489 has alternative start at 735..7418 and 741..7418. This line keeps both
                $matchHigh = $match;
            } else {
                $debug && print STDERR (" Old refcds is better. Skip\n");
                next;
            }

#            if ($match > $pct_conserved) {
#                $debug && print STDERR (" Take new refcds\n");
                $pct_conserved = $match;
                $refcds_set = $refset;
#            } else {
#                $debug && print STDERR (" Old refcds is better. Skip\n");
#            }

#        }
            if ($pct_conserved > $PCT_CONSERVED_MIN) {
        # Interesting case, in FJ588686: the CDS=246..12815 match with refcds very well, however, there is a large
        # gap at the middle, causing bl2seq to report 2 HSP, and neither is >0.67. This in turn causes the CDS to be
        # eliminated as candidate for annotation
        # 0.75 seems too high, eg AB005702, AF052446. 0.667 seems a good value --10/10/2012
                $found_polyprot = 1;
                $debug && print STDERR "$subn: \$ct=$ct \$cds=$s_cds_id \$refcds=$s_ref_id \$found_polyprot=$found_polyprot\n";
            } else {
                $comment = "CDS sequence doesn't match refseq with conserved=$pct_conserved";
                $debug && print STDERR "$subn: $comment  \$pct_conserved=$pct_conserved. Move to next feature\n\n";
                next;
            }

            # Found a polyprotein CDS, now collect all mat_peptide and sig_peptide following it
#            $debug && print STDERR ("$subn: Now look for any existing mat_peptide in genbank\n");
            my $ct_cds = $ct;
            my $matps = [];
            my %allowed_tags = ( 'mat_peptide' => 1, 'sig_peptide' => 1, 'misc_feature' => 1,);
            while ($ct<$#{$feats} && $allowed_tags{$feats->[$ct+1]->primary_tag}) {
                ++$ct;
                next if ($feats->[$ct]->primary_tag eq 'misc_feature');
                my $feat = $feats->[$ct];
                $debug && print STDERR "$subn: \$ct=$ct \$feat=".$feat->primary_tag." ".$feat->location->to_FTstring."\n";
                push @$matps, $feats->[$ct];
            }
#            $debug && print STDERR ("$subn: Done look for existing mat_peptide in genbank, found ".scalar(@$matps)."\n");

        # check if this is a new CDS, e.g. M55506 has a 2nd CDS that's part of 1st one with identical sequence
        # However, NC_004718 has a 2nd CDS that mostly same to 1st, but with small segment different at the end
            $debug && print STDERR "$subn: Check if new cds is just part of any existing cds, by (\$s1 =~ /\$s2/i)\n";
            my $seen=0;

            $debug && print STDERR "$subn: CDS=$s_cds_id";
            $seen = 0;
            if (!$seen) {
              $debug && print STDERR " is NOT part of any existing cds\n";
              push @{$polyprots->{$refcds_set->[0]}}, [$cds, @$matps]; # $refset->[0] is CDS_id of refseq
#              $polyprots->{$refcds_set->[0]}->[0] = [$cds, @$matps]; # $refset->[0] is CDS_id of refseq #2015Apr30
#              $bestMatch = [$cds, @$matps]; # $refset->[0] is CDS_id of refseq #2015Apr30
              $num_cds++;
            } else {
              $debug && print STDERR " IS part of an existing cds, see above for action taken\n";
            }

#        push @{$polyprots->{$refcds_set->[0]}}, [$cds, @$matps]; # $refset->[0] is CDS_id of refseq
#        $num_cds++;
            $debug && print STDERR "$subn: \$#polyprots=". keys(%$polyprots)."\n";
#            $debug && print STDERR "$subn: \$polyprots=\n".Dumper($polyprots)."end of \$polyprots\n\n";

            print STDERR "$subn: \$ct=$ct polyprotein: $acc:$ct_cds=$s_cds_id #$num_cds for refseq=$s_ref_id|$refset->[0]\n";

        } # for (my $i = 0; $i<=$#{$refpolyprots}; $i++)
#        $polyprots->{$refcds_set->[0]}->[0] = $bestMatch if ($bestMatch); # $refset->[0] is CDS_id of refseq #2015Apr30
#        push @{$polyprots->{$refcds_set->[0]}}, $bestMatch if ($bestMatch); # $refset->[0] is CDS_id of refseq #2015Apr30
        $debug && print STDERR ("$subn: CDS=$s_cds_id Finished looking through ".scalar(@$refpolyprots)." refcds\n\n");
      for my $refcds (sort keys %$polyprots) {
        for my $cdss (@{$polyprots->{$refcds}}) {
            my $cds = $cdss->[0];
            $debug && print STDERR "$subn: summary: refcds=$refcds cds=$acc|".$cds->primary_tag.'|'.$cds->location->to_FTstring."\n";
            
        }
      }

    } #     for (my $ct = 0; $ct<=$#{$feats}; $ct++)

    my $n = scalar (keys %$polyprots);
    $debug && print STDERR "$subn: polyprotein CDS for acc=$acc \$n=$n\n";
    $debug && print STDERR "$subn: \$polyprots=\n".Dumper($polyprots)."end of \$polyprots\n\n";
    print STDERR "$subn: Summary:\n";
    for my $refcds (sort keys %$polyprots) {
        for my $cdss (@{$polyprots->{$refcds}}) {
            my $cds = $cdss->[0];
            print STDERR "$subn: summary: refcds=$refcds cds=$acc|".$cds->primary_tag.'|'.$cds->location->to_FTstring."\n";
            
        }
    }
#    $comment = "Found $#{$n} CDS for acc=$acc suitable for annotation";
    $comment = "In $acc, highest match is $matchHigh" if ($n<0);

    return ($polyprots, $num_cds, $comment);
} # sub get_polyprots


=head2 cmp_cds_clustalw
 Takes 2 CDS feature
 Compares the translation of 2 CDS by clustalw, see if the conserved length is at least 80% of the length of the shorter seq
 Return 1 if the 2 seqs are similar
=cut

sub cmp_cds_clustalw {
    my ($cds, $refcds, $outfile_name, $runMUSCLE, $exe_dir, $dir_path, $progs) = @_;

    my $debug = 0 || $debug_all;
    my $subn = 'cmp_cds_clustalw';

    my $match_cds = 0;
    my $pct_conserved = 0;
    my $emsgs = '';

    my @values = $refcds->get_tag_values('translation');
    my $s1 = $values[0];

    @values = $cds->get_tag_values('translation');
    my $s2 = $values[0];

    my $cds_all = [ ];
#    $debug && print STDERR "$subn: \$refcds=\n".Dumper($refcds)."end of \$refcds\n\n";
    if ($refcds->isa('Bio::PrimarySeq')) {
        $cds_all = [ $refcds ];
    } else {
        $cds_all = [ Annotate_Util::refcds2PrimarySeq($refcds) ];
        if ( 1 ) {
            Annotate_Def::add_extra_refCDS( $cds_all, $exe_dir); # to increase the stability of the alignment, eg. MERS
        }
    }

#    $debug && print STDERR "$subn: \$cds=\n".Dumper($cds)."end of \$cds\n\n";
    my ($f2, $ecode);
    if (# $cds->isa('Bio::SeqFeature::Generic') || 
        $cds->isa('Bio::PrimarySeq')) {
        $f2 = $cds;
    } else {
        ($f2, $ecode) = Annotate_Util::cds2PrimarySeq($cds);
    }
    if ($ecode) {
        print STDERR "$subn: \$ecode='$ecode'\n";
        return $pct_conserved;
    }
    push @$cds_all, $f2;
    for my $i (0 .. $#{$cds_all}) { $debug && print STDERR "$subn: #$i \$cds='". $cds_all->[$i]->display_name . "'\n"; }

    my $clustalw_result = Annotate_Util::run_1MSA($cds_all, $runMUSCLE, $outfile_name, $dir_path, $progs);
#    $debug && print STDERR "$subn: \$clustalw_result=\n".Dumper($clustalw_result)."end of \$clustalw_result\n\n";
if (1) {
    my $aln = $clustalw_result;
    my $refcds_id = Annotate_Util::msa_seqid('ref=', $refcds);
    my $aln_q = Annotate_Util::msa_get_aln( $aln, $refcds_id);
    my $cds_id = Annotate_Util::msa_seqid('ACC=', $cds);
    $debug && print STDERR "$subn: \$refcds_id='$refcds_id' \$cds_id='$cds_id'\n";

    my $aln_h = Annotate_Util::msa_get_aln( $aln, $cds_id);
#    $debug && print STDERR "$subn: \$aln_h=\n". Dumper($aln_h) . "End of \$aln_h\n";

    my $qseq = $aln_q->seq;
    my $hseq = $aln_h->seq;
#    $debug && print STDERR "$subn: \$refcds_id='$refcds_id' \nqseq='$qseq'\n";
#    $debug && print STDERR "$subn:    \$cds_id='$cds_id' \nhseq='$hseq'\n";

    # Trim leading gaps
    my $gap = 0;
    for (my $i=0; $i<length($hseq); $i++) {
        my $hc = substr($hseq, $i, 1);
        if ($i==0) {
            if ($hc eq '.') {
                $gap++;
            } else {
                last;
            }
        } elsif ($hc eq '.') {
            $gap++;
        } else {
            $qseq = substr($qseq, $gap, length($qseq)-$gap-1);
            $hseq = substr($hseq, $gap, length($hseq)-$gap-1);
            last;
        }
    }
    # Trim trailing gaps
    $gap = 0;
    for (my $i=length($hseq)-1; $i>=0; $i--) {
        my $hc = substr($hseq, $i, 1);
        if ($i==length($hseq)-1) {
            if ($hc eq '.') {
                $gap++;
            } else {
                last;
            }
        } elsif ($hc eq '.') {
            $gap++;
        } else {
            $qseq = substr($qseq, 0, length($qseq)-$gap);
            $hseq = substr($hseq, 0, length($hseq)-$gap);
            last;
        }
    }
#    $debug && print STDERR "$subn: \$refcds_id='$refcds_id' \nqseq='$qseq'\n";
#    $debug && print STDERR "$subn:    \$cds_id='$cds_id' \nhseq='$hseq'\n";

    my $cseq = ''; # build conserved sequence since clustalw doesn't have one
    my ($gap_mid, $gap_mid_total) = (0, 0); # For any significant gap in middle, need to take care of it
    my ($matched, $conserved, $identical) = (0, 0, 0);
    for (my $i=0; $i<length($qseq); $i++) {
        my $qc = substr($qseq, $i, 1);
        my $hc = substr($hseq, $i, 1);

        if ($qc eq '.') {
            next;
        } elsif ($hc eq '.') {
            $gap_mid++;
            next;
        } elsif ($gap_mid) {
            if ($gap_mid > 20) {
                $gap_mid_total += $gap_mid;
            }
            $gap_mid = 0;
        }
        $matched++;

        if ($qc eq $hc) {
            $identical++;
            $conserved++;
            $cseq .= $qc;
        } elsif (Annotate_Def::conservedAA($hc, $qc)) {
            $conserved++;
            $cseq .= '+';
        } else {
            $cseq .= ' ';
        
        }
    }
    $debug && print STDERR "$subn: \$matched=$matched\t\$gap_mid_total=$gap_mid_total\n";
#    $matched -= $gap_mid_total; # Subtract gap in the middle from total length of alignment
    if ( 0 ) {
        $matched = length($f2->seq); # Take the length of target CDS as base of convserved/identical
    }
    $debug && print STDERR "$subn: \$conserved=$conserved / $matched\t\$identical=$identical / $matched\n";
#    $debug && print STDERR "$subn: qseq='".length($qseq)."'\n";
#    $debug && print STDERR "$subn: \$refcds_id='$refcds_id' \nqseq='$qseq'\n";
#    $debug && print STDERR "$subn: cseq='".length($cseq)."'\n";
#    $debug && print STDERR "$subn: conserved='conserved' \ncseq='$cseq'\n";
#    $debug && print STDERR "$subn: hseq='".length($hseq)."'\n";
#    $debug && print STDERR "$subn:    \$cds_id='$cds_id' \nhseq='$hseq'\n";

    my $pct_conserved = $conserved / $matched;
    $debug && printf STDERR "$subn: \$pct_conserved=%6.4f\n", $pct_conserved;

    if ( 1 ) {
        return $pct_conserved;
    } else {
        exit(0);
    }
}

#    my $clustalw_result = Annotate_Util::run_bl2seq_search($s1, $s2, $debug);
    $s1 = length($s1);
    $s2 = length($s2);
    $debug && print STDERR "$subn: refcds: \$s1=$s1 cds: \$s2=$s2\n\n";

    # Look at the hit
    my $hit_cds = $clustalw_result->next_hit;
    $debug && print STDERR "$subn: \$hit_cds=\n".Dumper($hit_cds)."end of \$hit_cds\n\n";
    return 0.88;
    return $match_cds if (!defined($hit_cds));

    my $hsp_cds_conserved = 0;
    my $hsp_cds_hit_length = 0;
    my $CONSERVED_RESIDUES_REQUIRED = 0.66667; # The target must have 66.667% conserved, needed for cases like KC175410 -7/8/2013
#    my $CONSERVED_RESIDUES_REQUIRED = 0.75; # The target must have 75% residues conserved wrt refseq
    $CONSERVED_RESIDUES_REQUIRED = 0.400 if ($debug); # 0.75 seems too high, eg AB005702, AF052446. --10/10/2012
    while (my $hsp_cds = $hit_cds->next_hsp) {
        $debug && print STDERR "$subn: \$hsp_cds=\n".Dumper($hsp_cds)."end of \$hsp_cds\n\n";
        # ignore any HSP shorter than 10% of hit_length
        $debug && print STDERR "$subn: \$hsp_cds=".$hsp_cds->{'HIT_START'}."..".$hsp_cds->{'HIT_END'}." has     conserved=".$hsp_cds->{'CONSERVED'}."\n";
        if (($hsp_cds->{'HIT_END'} - $hsp_cds->{'HIT_START'} +1) < 0.2*$hsp_cds->{'HIT_LENGTH'}) {
            $debug && print STDERR "$subn: \$hsp_cds=".$hsp_cds->{'HIT_START'}."..".$hsp_cds->{'HIT_END'}." too short. Skip\n";
            next;
        }
        # ignore any HSP shorter than 66.7% of hit_length
#        if ($hsp_cds->{'CONSERVED'}  < 0.667*($hsp_cds->{'HIT_END'} - $hsp_cds->{'HIT_START'} +1)) {
        if ($hsp_cds->{'CONSERVED'}  < $CONSERVED_RESIDUES_REQUIRED*($hsp_cds->{'HIT_END'} - $hsp_cds->{'HIT_START'} +1)) {
            $debug && print STDERR "$subn: \$hsp_cds=".$hsp_cds->{'HIT_START'}."..".$hsp_cds->{'HIT_END'}." too few conserved=".$hsp_cds->{'CONSERVED'}.". Skip\n";
            next;
        }
        # get the accumulative conserved length, in case there is significant gap in the middle, such as 2nd CDS in FJ588686
        $hsp_cds_conserved += (exists($hsp_cds->{'CONSERVED'})) ? $hsp_cds->{'CONSERVED'} : $hsp_cds->{'num_conserved'};
        $hsp_cds_hit_length += $hsp_cds->{'HIT_END'} - $hsp_cds->{'HIT_START'} +1;
    }
    $debug && print STDERR "$subn: refcds: \$s1=$s1 \$s2=$s2 length of conserved=$hsp_cds_conserved\n";
    # This check assumes that the target CDS should never be longer than refcds, I think this is reasonable -gsun
    my $s3 = ($s1<$s2) ? $s1 : $s2;
    if ($hsp_cds_conserved >= $s3 * $CONSERVED_RESIDUES_REQUIRED) {
            $pct_conserved = $hsp_cds_conserved/$s2 if ($hsp_cds_conserved/$s2 > $pct_conserved);
            $pct_conserved = $hsp_cds_conserved/$s3 if ($hsp_cds_hit_length && $hsp_cds_conserved/$hsp_cds_hit_length > $pct_conserved);

            if ($match_cds) {
              print STDERR "$subn: WARNING: found more matched refcds when \$match_cds=$match_cds\n";
            } else {
              $match_cds = 1;
            }
            $debug && print STDERR "$subn: CDS length=$s2 refcds length=$s1 conserved=$hsp_cds_conserved.";
            $debug && print STDERR " Conserved meets required $CONSERVED_RESIDUES_REQUIRED\n";
            $debug && print STDERR "$subn: \$match_cds=$match_cds\n";
#            last;
    } else {
            my $cs = $cds->location->to_FTstring;
            my $rs = $refcds->location->to_FTstring;
            $emsgs .= "$subn: CDS=$cs refcds=$rs length=$s2 conserved=$hsp_cds_conserved.";
            $emsgs .= " Conserved doesn't meet required $CONSERVED_RESIDUES_REQUIRED\n";
            $emsgs .= "$subn: \$match_cds=$match_cds. See following \$hit_cds\n";
            $emsgs .= "$subn: \$hit_cds=\n".Dumper($hit_cds)."end of \$hit_cds\n\n";
    }
    $debug && (!$match_cds) && print STDERR "$emsgs";
    $debug && print STDERR "$subn: \$pct_conserved=$pct_conserved, leaving subroutine\n";

#    return $match_cds;
    return $pct_conserved;
} # sub cmp_cds_clustalw


=head2 cmp_cds_bl2seq
 Takes 2 CDS feature
 Compares the translation of 2 CDS by bl2seq, see if the conserved length is at least 80% of the length of the shorter seq
 Return 1 if the 2 seqs are similar
=cut

sub cmp_cds_bl2seq {
    my ($cds, $refcds) = @_;

    my $debug = 0 || $debug_all;
    my $subn = 'cmp_cds_bl2seq';

    my $match_cds = 0;
    my $pct_conserved = 0;
    my $emsgs = '';

    my @values = $refcds->get_tag_values('translation');
    my $s1 = $values[0];

    @values = $cds->get_tag_values('translation');
    my $s2 = $values[0];
    my $bl2seq_result = &run_bl2seq_search($s1, $s2, $debug);

if ( 1 ) { # count the conserved pairs
    
    
    
    
    
    
    
    
    
    
    
    
    
   exit(0); 
}

    $s1 = length($s1);
    $s2 = length($s2);
    $debug && print STDERR "$subn: \$bl2seq_result=\n".Dumper($bl2seq_result)."end of \$bl2seq_result\n\n";
    $debug && print STDERR "$subn: refcds: \$s1=$s1 cds: \$s2=$s2\n\n";

    # Look at the hit
    my $hit_cds = $bl2seq_result->next_hit;
    $debug && print STDERR "$subn: \$hit_cds=\n".Dumper($hit_cds)."end of \$hit_cds\n\n";
    return 0.88;
    return $match_cds if (!defined($hit_cds));

    my $hsp_cds_conserved = 0;
    my $hsp_cds_hit_length = 0;
    my $CONSERVED_RESIDUES_REQUIRED = 0.66667; # The target must have 66.667% conserved, needed for cases like KC175410 -7/8/2013
#    my $CONSERVED_RESIDUES_REQUIRED = 0.75; # The target must have 75% residues conserved wrt refseq
    $CONSERVED_RESIDUES_REQUIRED = 0.400 if ($debug); # 0.75 seems too high, eg AB005702, AF052446. --10/10/2012
    while (my $hsp_cds = $hit_cds->next_hsp) {
        $debug && print STDERR "$subn: \$hsp_cds=\n".Dumper($hsp_cds)."end of \$hsp_cds\n\n";
        # ignore any HSP shorter than 10% of hit_length
        $debug && print STDERR "$subn: \$hsp_cds=".$hsp_cds->{'HIT_START'}."..".$hsp_cds->{'HIT_END'}." has     conserved=".$hsp_cds->{'CONSERVED'}."\n";
        if (($hsp_cds->{'HIT_END'} - $hsp_cds->{'HIT_START'} +1) < 0.2*$hsp_cds->{'HIT_LENGTH'}) {
            $debug && print STDERR "$subn: \$hsp_cds=".$hsp_cds->{'HIT_START'}."..".$hsp_cds->{'HIT_END'}." too short. Skip\n";
            next;
        }
        # ignore any HSP shorter than 66.7% of hit_length
#        if ($hsp_cds->{'CONSERVED'}  < 0.667*($hsp_cds->{'HIT_END'} - $hsp_cds->{'HIT_START'} +1)) {
        if ($hsp_cds->{'CONSERVED'}  < $CONSERVED_RESIDUES_REQUIRED*($hsp_cds->{'HIT_END'} - $hsp_cds->{'HIT_START'} +1)) {
            $debug && print STDERR "$subn: \$hsp_cds=".$hsp_cds->{'HIT_START'}."..".$hsp_cds->{'HIT_END'}." too few conserved=".$hsp_cds->{'CONSERVED'}.". Skip\n";
            next;
        }
        # get the accumulative conserved length, in case there is significant gap in the middle, such as 2nd CDS in FJ588686
        $hsp_cds_conserved += (exists($hsp_cds->{'CONSERVED'})) ? $hsp_cds->{'CONSERVED'} : $hsp_cds->{'num_conserved'};
        $hsp_cds_hit_length += $hsp_cds->{'HIT_END'} - $hsp_cds->{'HIT_START'} +1;
    }
    $debug && print STDERR "$subn: refcds: \$s1=$s1 \$s2=$s2 length of conserved=$hsp_cds_conserved\n";
    # This check assumes that the target CDS should never be longer than refcds, I think this is reasonable -gsun
    my $s3 = ($s1<$s2) ? $s1 : $s2;
    if ($hsp_cds_conserved >= $s3 * $CONSERVED_RESIDUES_REQUIRED) {
            $pct_conserved = $hsp_cds_conserved/$s2 if ($hsp_cds_conserved/$s2 > $pct_conserved);
            $pct_conserved = $hsp_cds_conserved/$s3 if ($hsp_cds_hit_length && $hsp_cds_conserved/$hsp_cds_hit_length > $pct_conserved);

            if ($match_cds) {
              print STDERR "$subn: WARNING: found more matched refcds when \$match_cds=$match_cds\n";
            } else {
              $match_cds = 1;
            }
            $debug && print STDERR "$subn: CDS length=$s2 refcds length=$s1 conserved=$hsp_cds_conserved.";
            $debug && print STDERR " Conserved meets required $CONSERVED_RESIDUES_REQUIRED\n";
            $debug && print STDERR "$subn: \$match_cds=$match_cds\n";
#            last;
    } else {
            my $cs = $cds->location->to_FTstring;
            my $rs = $refcds->location->to_FTstring;
            $emsgs .= "$subn: CDS=$cs refcds=$rs length=$s2 conserved=$hsp_cds_conserved.";
            $emsgs .= " Conserved doesn't meet required $CONSERVED_RESIDUES_REQUIRED\n";
            $emsgs .= "$subn: \$match_cds=$match_cds. See following \$hit_cds\n";
            $emsgs .= "$subn: \$hit_cds=\n".Dumper($hit_cds)."end of \$hit_cds\n\n";
    }
    $debug && (!$match_cds) && print STDERR "$emsgs";
    $debug && print STDERR "$subn: \$pct_conserved=$pct_conserved, leaving subroutine\n";

#    return $match_cds;
    return $pct_conserved;
} # sub cmp_cds_bl2seq


=head2 get_new_translation
Takes a new feature object, and a CDS
 return translation
=cut

sub get_new_translation {
    my ($feat, $cds) = @_;

    my $debug = 0;
#    $debug = 0 || $debug_all;
    my $subn = 'get_new_translation';
#    my $translation = '';
    my $translation = undef; # emptry string '' suppresses error msg, but cause other problems
#    my $acc = $cds->seq->accession_number;
    my $acc = '';
    if ($cds->seq) {
        $acc = $cds->seq->accession_number
    } else {
        for my $note ($cds->get_tag_values('note')) {
            $acc = $1 if ($note =~ m/^Desc:.+[|]ACC=([^|]+)[|]/i);
            $debug && print STDERR "$subn: \$acc='$acc' \$note='$note'\n";
        }
    }

    $debug && print STDERR "$subn: \$feat=\n".Dumper($feat)."End of \$feat\n\n";
    if ( 0 ) {
        my $overhang1 = ($feat->location->start - $cds->location->start) % 3; # any overhang at start?
        my $overhang2 = ($feat->location->end   - $cds->location->start + 1) % 3; # any overhang at end?
        if ($overhang1 || $overhang2 ) {
            $debug && print STDERR "$subn: ERROR: \$acc=$acc \$feat='".$feat->location->to_FTstring."' \$cds='".$cds->location->to_FTstring."'\n";
            $debug && print STDERR "$subn: ERROR: \$acc=$acc \$overhang1=$overhang1 \$overhang2=$overhang2\n";
#            croak "$subn: ERROR: \$overhang1=$overhang1 \$overhang2=$overhang2";
        }
    }

    my $s = '';
#    if ($feat->strand() != -1) {
#        $s = Annotate_Math::get_dna_byloc_complement( $feat, $cds->{_gsf_seq}->seq);
#    } else {
    if ($cds->{_gsf_seq}) {
        $s = Annotate_Math::get_dna_byloc( $feat, $cds->{_gsf_seq}->seq);
    } else {
        for my $v ($cds->get_tag_values('translation')) { next if (!$v); $acc = $v; last; }
    }
#    }
    # Add generic N to end of DNA, if the last codon has only 2 nucleotides but is able to define AA - 03/19/2013
    if (length($s) % 3 ==2 && substr($s, length($s)-2, 2) =~/([TU]C|C[TU]|CC|CG|AC|G[TU]|GC|GG)/i ) {
        $s .= 'N';
        $debug && print STDERR "$subn: Appended generic nucleotide 'N' to \$s=".length($s)." \$s='$s'\n";
    }
    $debug && print STDERR "$subn: \$s  ='$s'\n";
    my $f = Bio::PrimarySeq->new(
                         -seq      => $s,
                         -id       => '',	# id can't contain space
                         -alphabet => 'dna'
                                );
    $debug && print STDERR "$subn: \$feat->strand() =".$feat->strand()."\n";
#    $f->revcom() if ($feat->strand() == -1);
    $debug && print STDERR "$subn: \$f=\n".Dumper($f)."End of \$f\n\n";

# The flag '-complete => 1' is added to accommodate the change in BioPerl version. -3/08/2013
# The '-complete' flag failed, since it translates CTG at the start of CDS to M, instead of L (eg. AB127995). -3/19/2013
if (0) {
    $s = $f->translate(-complete => 1); 
    $debug && print STDERR "$subn: \$f=\n".Dumper($f)."End of \$f\n\n";
    $debug && print STDERR "$subn: \$s=\n".Dumper($s)."End of \$s\n\n";
    $s = $s->seq; 
} else {
    $s = $f->translate(); 
    $debug && print STDERR "$subn: \$f=\n".Dumper($f)."End of \$f\n\n";
    $debug && print STDERR "$subn: \$s=\n".Dumper($s)."End of \$s\n\n";
    $s = $s->seq; 
}
    $s =~ s/[*]$// if ($s =~ /[*]$/);
    $s =~ s/[*]/./g if ($s =~ /[*]/);
#    $s =~ s/[X]/./ig if ($s =~ /[X]/i);
    my $s1 = $s;
    $s1 =~ s/[X]/./ig if ($s1 =~ /[X]/i);
    $debug && print STDERR "$subn: \$s  ='$s'\n";
    $debug && print STDERR "$subn: \$s  =".length($s)."\n";

    # Only keep the annotated mat_peptide that confirms to CDS. It should, just to double check
    my @parent_cds_seq = ();
    if ($cds->has_tag('translation')) {
      @parent_cds_seq = $cds->get_tag_values('translation');
#      $parent_cds_seq[0] =~ s/[X]/./;
    } else {
        my $s1 = Annotate_Math::get_dna_byloc( $cds, $cds->{_gsf_seq}->seq);
        if (length($s1) % 3 ==2 && substr($s1, length($s1)-2, 2) =~/([TU]C|C[TU]|CC|CG|AC|G[TU]|GC|GG)/i ) {
            $s1 .= 'N';
            print STDERR "$subn: Appended generic nucleotide 'N' to \$s1=".length($s1)." \$s1='$s1'\n";
        }
        $debug && print STDERR "$subn: \$s1  ='$s1'\n";
        my $f = Bio::PrimarySeq->new(
                         -seq      => $s1,
                         -id       => '',	# id can't contain space
                         -alphabet => 'dna'
                                );
        $debug && print STDERR "$subn: \$cds->strand() =".$cds->strand()."\n";
#        $f->revcom() if ($feat->strand() == -1);
        $debug && print STDERR "$subn: \$f=\n".Dumper($f)."End of \$f\n\n";
        $s1 = $f->translate(); 
        $debug && print STDERR "$subn: \$f=\n".Dumper($f)."End of \$f\n\n";
        $debug && print STDERR "$subn: \$s=\n".Dumper($s)."End of \$s\n\n";
        $s1 = $s1->seq; 
        $s1 =~ s/[*]$// if ($s1 =~ /[*]$/);
        $s1 =~ s/[*]/./g if ($s1 =~ /[*]/);
        $s1 =~ s/[X]/./ig if ($s1 =~ /[X]/i);
        @parent_cds_seq = ($s1);
    }
    $debug && print STDERR "$subn: \$cds  =$parent_cds_seq[0]\n";
    $debug && print STDERR "$subn: \$cds  =".length($parent_cds_seq[0])."\n";
    if ($parent_cds_seq[0] =~ /($s1)/i
           || $s1 =~ /$parent_cds_seq[0]/i
           || Annotate_Verify::diff_2str( $s1, substr($parent_cds_seq[0], 0, length($s1))) =~ /^L[.]+$/i) {
        $translation = $s;
        $debug && print STDERR "$subn: \$s  ='$s'\n";
        $debug && print STDERR "$subn: \$cds='".$parent_cds_seq[0]."'\n";
        $debug && print STDERR "$subn: diff='".Annotate_Verify::diff_2str( $s, $parent_cds_seq[0])."'\n";
    } else {
        print STDERR "$subn: ERROR: \$acc=$acc translation for feature=".$feat->location->to_FTstring." doesn't match CDS.\n";
        print STDERR "$subn: \$s  ='$s'\n";
        print STDERR "$subn: \$cds='".$parent_cds_seq[0]."'\n";
        print STDERR "$subn: diff='".Annotate_Verify::diff_2str( $s, $parent_cds_seq[0])."'\n";
#        return undef;
    }

    $debug && print STDERR "$subn: \$acc=$acc \$translation  ='$translation'\n";
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
 note="Desc:CDS=GI:11559447|Loc=311..883|symbol=unk|mat_peptide=core protein"

=cut

sub assemble_new_feature {
    my ($refcds, $reffeat, $cds, $aln, $refcds_id, $cds_id, $note,$exe_dir) = @_;

    my $debug = 0 || $debug_all;
    my $subn = 'assemble_new_feature';

    my $acc = $cds->seq->accession_number;
    my $aln_q = Annotate_Util::msa_get_aln( $aln, $refcds_id);
    my $aln_h = Annotate_Util::msa_get_aln( $aln, $cds_id);
    $debug && print STDERR "$subn: \$aln_q=\n".Dumper($aln_q)."End of \$aln_q\n\n";
    $debug && print STDERR "$subn: \$aln_h=\n".Dumper($aln_h)."End of \$aln_h\n\n";
    if (!$aln_q ) {
        $debug && print STDERR "$subn: \$aln_q is empty\n";
        return undef;
    }
    if ( !$aln_h) {
        $debug && print STDERR "$subn: \$aln_h is empty\n";
        return undef;
    }
#     my ($reffeat, $refcds, $cds, $aln, $aln_q, $aln_h) = @_;
    my ($loc2, $errcode) = Annotate_Math::msa_get_feature_loc(
                                 $reffeat,
                                 $refcds,
                                 $cds,
                                 $aln,
                                 $aln_q,
                                 $aln_h,
                                 );
    $debug && print STDERR "$subn: after Annotate_Math::msa_get_feature_loc \$loc2=\n".Dumper($loc2)."End of \$loc2\n\n";

    Annotate_Util::checkHeadTailSimilarity($loc2, $errcode, $refcds, $reffeat, $cds, $aln, $refcds_id, $cds_id, $note,$exe_dir);
    $debug && print STDERR "$subn: after checkHeadTailSimilarity \$loc2=\n".Dumper($loc2)."End of \$loc2\n\n";

    if (!exists($errcode->{alnstart}) || !$errcode->{alnstart} || !exists($errcode->{alnend}) || !$errcode->{alnend}) {
        $debug && print STDERR "$subn: ERROR \$loc2 has problem with alnstart=$errcode->{alnstart} or alnend=$errcode->{alnend}\n";
    }
    if (!$loc2) {
        $debug && print STDERR "$subn: \$loc2 is empty\n";
        return undef;
    }
    if ($errcode->{OutsideCDS}==1) {
        print STDERR "$subn: ACC=$acc \$loc2=".$loc2->to_FTstring." doesn't match \$cds=".$cds->location->to_FTstring."\n";
        return undef;
    }
    if ($errcode->{long_internal_gap}==1) {
        print STDERR "$subn: ACC=$acc \$loc2=".$loc2->to_FTstring." has too long internal gap in alignment, discard.\n";
        return undef;
    }
=head2
#    if (($loc2->end-$loc2->start+1) < 0.15 * ($reffeat->location->end-$reffeat->location->start+1)) {
    # kick out those seqs that are too short, sometimes caused by artifacts in MUSCLE alignment
    # previously set at 50% of reffeat, which turns out too high. 10 AA might be a good one --10/22/2012
    # Also, switching to CLUSTALW should reduce MSA artifacts
    my $len_min = ($reffeat->location->end-$reffeat->location->start+1) * 0.33;
    $len_min = $MIN_CDS_LENGTH if ($len_min > $MIN_CDS_LENGTH);
    if (($loc2->end-$loc2->start+1) < $len_min) {
        if ($loc2->end-$loc2->start+1>2) {
          print STDERR "$subn: ACC=$acc \$loc2=".$loc2->to_FTstring."=".($loc2->end-$loc2->start+1)." is too short for ref=".$reffeat->location->to_FTstring."=".($reffeat->location->end-$reffeat->location->start+1).", discard.\n";
        }
        return undef;
    }
=cut
    my $feat;
    $feat = Bio::SeqFeature::Generic->new();
    $feat->primary_tag($reffeat->primary_tag); # Take the primary_tag from reffeat
    $debug && print STDERR "$subn: \$feat is ".ref($feat)." \n";

    $feat->location($loc2);
    $debug && print STDERR "$subn: \$loc2=\n".Dumper($loc2)."end of \$loc2\n";
    $debug && print STDERR "$subn: \$feat=\n".Dumper($feat)."end of \$feat\n\n";

#    $feat->attach_seq($cds->seq);

    $feat->add_tag_value('note', $note);
#    $feat->strand($reffeat->strand()); # set the strand same as reffeat
#    $feat->strand($cds->strand()) if ($cds->strand==-1); # in case if target CDS differs from refseq, such as HQ647170 and HQ647168
    $feat->strand($cds->strand()); # in case if target CDS differs from refseq, such as HQ647170 and HQ647168
    $debug && print STDERR "$subn: \$cds=".$cds->strand." \$feat=".$feat->strand."\n";

    if ($reffeat->has_tag('product')) {
        my $tag = [$reffeat->get_tag_values('product')];
        # take care of the wrong product 2C (4195..5157) in NC_012800
        if ($reffeat->seq->accession_number eq 'NC_012800' && $reffeat->location->to_FTstring eq '4195..5157' && $tag->[0] eq '3C') {
            print STDERR "$subn: refseq accession=".$reffeat->seq->accession_number.
               " mat_peptide=".$reffeat->location->to_FTstring.
               " found erranous 3C, changed to '2C'\n";
            $tag->[0] = '2C';
        }
        # take care of the wrong product 2C (3801..4763) in NC_012801
        if ($reffeat->seq->accession_number eq 'NC_012801' && $reffeat->location->to_FTstring eq '3801..4763' && $tag->[0] eq '3C') {
            print STDERR "$subn: refseq accession=".$reffeat->seq->accession_number.
               " mat_peptide=".$reffeat->location->to_FTstring.
               " found erranous 3C, changed to '2C'\n";
            $tag->[0] = '2C';
        }
        # take care of the wrong product 2C (3787..4749) in NC_012802
        if ($reffeat->seq->accession_number eq 'NC_012802' && $reffeat->location->to_FTstring eq '3787..4749' && $tag->[0] eq '3C') {
            print STDERR "$subn: refseq accession=".$reffeat->seq->accession_number.
               " mat_peptide=".$reffeat->location->to_FTstring.
               " found erranous 3C, changed to '2C'\n";
            $tag->[0] = '2C';
        }
        $feat->add_tag_value('product', $tag->[0]);
    }

    if ($cds->has_tag('codon_start')) {
        my $codon_start = [$cds->get_tag_values('codon_start')];
        $codon_start = $codon_start->[0];
        $debug && print STDERR "$subn: \$cds=".$cds->location->to_FTstring." has codon_start=$codon_start\n";
        if ($cds->strand!=-1) {
          if ($codon_start>=1 && $cds->location->start == $feat->location->start - $codon_start +1) {
        # codon_start is handled separately, adjust the start of the first mat_peptide, and set codon_start, V1.1.5 --10/9/2012
            if ($cds->location->isa('Bio::Location::Split')) {
                my $locs = [ $cds->location->sub_Location ];
                $locs->[0]->start($locs->[0]->start - $codon_start +1);
                $debug && print STDERR "$subn: Added codon_start=$codon_start to start=".$locs->[0]->to_FTstring."\n";
            } else {
                $feat->location->start($feat->location->start - $codon_start +1);
                $debug && print STDERR "$subn: Added codon_start=$codon_start to ".$feat->location->to_FTstring."\n";
            }
            $feat->add_tag_value('codon_start', $codon_start) if ($codon_start>1);
            $debug && print STDERR "$subn: After adding codon_start \$feat=".$feat->location->to_FTstring."\n";
          }
        } else {
          if ($codon_start>=1 && $cds->location->end == $feat->location->end + $codon_start -1) {
        # codon_start is handled separately, adjust the start of the first mat_peptide, and set codon_start, V1.1.5 --10/9/2012
            if ($cds->location->isa('Bio::Location::Split')) {
                my $locs = [ $cds->location->sub_Location ];
                $locs->[$#{$locs}]->end($locs->[$#{$locs}]->end + $codon_start -1);
                $debug && print STDERR "$subn: Added codon_start=$codon_start to start=".$locs->[0]->to_FTstring."\n";
            } else {
                $feat->location->end($feat->location->end + $codon_start -1);
                $debug && print STDERR "$subn: Added codon_start=$codon_start to ".$feat->location->to_FTstring."\n";
            }
            $feat->add_tag_value('codon_start', $codon_start) if ($codon_start>1);
            $debug && print STDERR "$subn: After adding codon_start \$feat=".$feat->location->to_FTstring."\n";
          }
        }
    }
    $debug && print STDERR "$subn: \$feat=\n".Dumper($feat)."end of \$feat\n\n";

    # now get translation
    my $translation = '';
    if ($cds) {
        my $s = Annotate_Util::get_new_translation( $feat, $cds);
        if ($s) {
            $translation = $s;
        } else {
            $debug && print STDERR "$subn: \$cds =\n".Dumper($cds )."End of \$cds \n\n";
            $debug && print STDERR "$subn: \$feat=\n".Dumper($feat)."End of \$feat\n\n";
            print STDERR "$subn: WARNING: translation for ".$feat->primary_tag.":".$feat->location->to_FTstring." returned null.\n";
            $debug && print STDERR "$subn: \$cds='".$cds->seq->translate->seq."'\n";
            return undef;
        }

    }
    $translation =~ s/[.]/X/ig;
    $feat->add_tag_value('translation', $translation);
#    print STDERR "$subn: \$feat=\n".Dumper($feat)."end of \$feat\n\n";

#    if (($loc2->end-$loc2->start+1) < 0.15 * ($reffeat->location->end-$reffeat->location->start+1)) {
    # kick out those seqs that are too short, sometimes caused by artifacts in MUSCLE alignment
    # previously set at 50% of reffeat, which turns out too high. 10 AA might be a good one --10/22/2012
    # Also, switching to CLUSTALW should reduce MSA artifacts
    my $len_min = ($reffeat->location->end-$reffeat->location->start+1) * 0.33;
    $len_min = $MIN_CDS_LENGTH if ($len_min > $MIN_CDS_LENGTH);
    if (($loc2->end-$loc2->start+1) < $len_min) {
        if ($loc2->end-$loc2->start+1>2) {
          print STDERR "$subn: ACC=$acc \$loc2=".$loc2->to_FTstring."=".($loc2->end-$loc2->start+1).":'$translation' is too short for ref=".$reffeat->location->to_FTstring."=".($reffeat->location->end-$reffeat->location->start+1).", discard.\n";
        }
        return undef;
    }

    if ($cds->has_tag('locus_tag')) {
        my $tag = [$cds->get_tag_values('locus_tag')];
        $feat->add_tag_value('locus_tag', $tag->[0]);
    }

    my $gene_symbol = Annotate_Def::get_gene_symbol( $reffeat,$exe_dir);
    $debug && print STDERR "$subn: ready to enter sub get_feature_id_desc\n";
    my ($id, $desc) = Annotate_Util::get_feature_id_desc( $feat, $errcode, $gene_symbol, $cds, $refcds, $reffeat);
    my @notes = ();
    $notes[0] = 'Desc:'.$desc;
    $notes[1] = "MSA=$errcode->{alnstart}..$errcode->{alnend}";
    $feat->add_tag_value('note', @notes);

    $debug && print STDERR "$subn: \$loc=\n".Dumper($reffeat->location)."End of \$loc\n\n";
    $debug && print STDERR "$subn: \$feat=\n".Dumper($feat)."End of \$feat\n\n";

    return ($feat);
} # sub assemble_new_feature


=head2 checkHeadTailSimilarity
Takes a new location, checks if there is enough similarity at either end between the RefSeq and target
 Requires a minimum similarity of 2 out of 5+5 AA at head+tail of the new mat_peptide
 change $loc to undef if the requirement is not met
=cut

sub checkHeadTailSimilarity {
    my ($loc2, $errcode, $refcds, $reffeat, $cds, $aln, $refcds_id, $cds_id, $note,$exe_dir) = @_;

    my $debug = 0 || $debug_all;
    my $subn = 'checkHeadTailSimilarity';

    if (!exists($errcode->{low_similarity})) {
        $errcode->{low_similarity} = 0;
    }
    my $min_similar_within5 = 2;
    my $MAX_LENGTH = 05;
    my $changed = 0;
#    $debug && print STDERR "$subn: \$cds=\n".Dumper($cds)."end of \$cds\n";
    $debug && print STDERR "$subn: \$errcode=\n".Dumper($errcode)."end of \$errcode\n";

  # This for loop is used to take care of cases where 2 mat_peptides end at same
  # location and there is a gap, eg. preM and M
    my $acc = $cds->seq->accession_number;

    # For loc, find the N-term and C-term
    my ($str, $l, $n, $f);
    my $ts = ['', ''];
    my $rs = ['', ''];
    my $aln_q = Annotate_Util::msa_get_aln( $aln, $refcds_id);
    my $aln_h = Annotate_Util::msa_get_aln( $aln, $cds_id);

    $debug && print STDERR "$subn: \$loc2=\n".Dumper($loc2)."end of \$loc2\n";
    $rs->[0] = substr($aln_q->seq, $errcode->{alnstart}-1, $MAX_LENGTH);
    $ts->[0] = substr($aln_h->seq, $errcode->{alnstart}-1, $MAX_LENGTH);
    $rs->[1] = substr($aln_q->seq, $errcode->{alnend}-$MAX_LENGTH, $MAX_LENGTH);
    $ts->[1] = substr($aln_h->seq, $errcode->{alnend}-$MAX_LENGTH, $MAX_LENGTH);

    # Determine the reffeat for each feature
    $debug && print STDERR "$subn: ACC=$acc \$rs='@$rs'\n";
    $debug && print STDERR "$subn: ACC=$acc \$ts='@$ts'\n";

    # Update the dscription to include GapAtCleavage

    # The mat_peptide have to have some minimum similarity w.r.t. the ref mat_peptide. -May 17, 2013

    my $similar_first5 = 0;
    my $similar_last5 = 0;

    $n = (length($rs->[0])<=length($ts->[0])) ? length($rs->[0]) : length($ts->[0]);
    my ($rc, $tc);
    foreach my $j (reverse 1 .. $n) {
        $rc = uc substr($rs->[0], length($rs->[0])-$j, 1);
        $tc = uc substr($ts->[0], length($ts->[0])-$j, 1);
        next if ($rc eq '.' || $tc eq '.');
        $similar_first5++ if ($j<=5 && Annotate_Def::similarAA($rc, $tc));
        $debug && print STDERR "$subn: \$j=$j \$rc-\$tc:$rc-$tc \$similar_first5=$similar_first5 \$similar_last5=$similar_last5\n";
    }
    $rc = uc substr($rs->[0], length($rs->[0])-1, 1);
    $tc = uc substr($ts->[0], length($ts->[0])-1, 1);

    $n = (length($rs->[1])<=length($ts->[1])) ? length($rs->[1]) : length($ts->[1]);
    foreach my $j (0 .. $n-1) {
        $rc = uc substr($rs->[1], $j, 1);
        $tc = uc substr($ts->[1], $j, 1);
        next if ($rc eq '.' || $tc eq '.');
        $similar_last5++ if ($j<=5 && Annotate_Def::similarAA($rc, $tc));
        $debug && print STDERR "$subn: \$j=$j \$rc-\$tc:$rc-$tc \$similar_first5=$similar_first5 \$similar_last5=$similar_last5\n";
    }
    my $sim = $similar_first5 + $similar_last5;
    $debug && print STDERR "$subn: ACC=$acc \$rs='@$rs'\n";
    $debug && print STDERR "$subn: ACC=$acc \$ts='@$ts' there are $sim/". $MAX_LENGTH*2 ." similar AA";
    $debug && print STDERR " among the first $MAX_LENGTH and last $MAX_LENGTH residues\n";

    if ($sim < $min_similar_within5) {
        print STDERR "$subn: ACC=$acc \$rs='@$rs'\n";
        print STDERR "$subn: ACC=$acc \$ts='@$ts' there are $sim/". $MAX_LENGTH*2 ." similar AA";
        print STDERR " among the first $MAX_LENGTH and last $MAX_LENGTH residues\n";
        print STDERR "$subn: ERROR: less than required $min_similar_within5 similar residues, changing \$loc2=undef\n";
        $loc2 = undef;
        $errcode->{low_similarity} = 1;
    }
    $debug && print STDERR "$subn: \$errcode=\n".Dumper($errcode)."end of \$errcode\n";
    return $changed;
} # sub checkHeadTailSimilarity


=head2 get_feature_id_desc

Takes $feat, $errcode, $refcds, $reffeat, $cds, $note
 Return a suitable id and desc for the new feature
 note="Annotated by VIPRBRC, MUSCLE, refseq=$reffeat->accession_number"
 note="Desc:CDS=GI:11559447|Loc=311..883|symbol=unk|mat_peptide=core protein"

=cut

sub get_feature_id_desc {
    my ($feat, $errcode, $gene_symbol, $cds, $refcds, $reffeat) = @_;

    my $debug = 0 || $debug_all;
    my $subn = 'get_feature_id_desc';

    $debug && print STDERR "$subn: \$feat=\n".Dumper($feat)."end of \$feat\n";
    $debug && print STDERR "$subn: \$gene_symbol='$gene_symbol'\n";
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
#    $id = "CDS=unknown|" if (($id !~ /^CDS=/) ); # indicates a potential problem
    $id = "CDS=".$cds->location->to_FTstring."|" if (($id !~ /^CDS=/) ); # indicates a potential problem
    $debug && print STDERR "$subn: \$id='$id'\n";

    # add the accession of target genome
    if ( 1 ) {
        my $desc = '';
        $desc .= 'src=MSA|'; # Add method to description
        $desc .= 'ACC='. $cds->seq->accession_number; # Add accession to description
        $desc .= '|Ver='. $cds->seq->accession_number; # Add accession to description
        $desc = $desc .'.'.$cds->entire_seq->{'_version'}; # Add version of accession to description
        $debug && print STDERR "$subn: \$cds=\n".Dumper($cds)."End of \$cds\n\n";
        $debug && print STDERR "$subn: \$seq=\n".Dumper($cds->entire_seq)."End of \$seq\n\n";
        $id = $desc .'|'. $id;

        # add accession.version of refseq
        $desc = 'ref='. $refcds->seq->accession_number; # Add accession to description
        $desc .= '.'.$refcds->entire_seq->{'_version'}; # Add version of accession to description
#        $debug && print STDERR "$subn: \$refcds=\n".Dumper($refcds)."End of \$refcds\n\n";
#        $debug && print STDERR "$subn: \$refseq=\n".Dumper($refcds->entire_seq)."End of \$refseq\n\n";
        $id .= $desc .'|';

        # add reffeat location, this is used in part to detect any gaps between annotated mat_peptides
#        $desc = 'reffeat='. $reffeat->location->start.'..'.$reffeat->location->end;
        $desc = 'RM='. $reffeat->location->start.'..'.$reffeat->location->end;
        $debug && print STDERR "$subn: \$desc=$desc \$reffeat=\n".Dumper($reffeat)."End of \$reffeat\n\n";
#        $debug && print STDERR "$subn: \$refseq=\n".Dumper($refcds->entire_seq)."End of \$refseq\n\n";
        $id .= $desc .'|';
    }

    # add the range of mat_peptide
    $id .= 'Loc='. &get_DNA_loc($feat) .'|';
    $debug && print STDERR "$subn: \$id='$id'\n";

    if ($feat->has_tag('codon_start')) {
        my $codon_start = [$feat->get_tag_values('codon_start')];
        $codon_start = $codon_start->[0];
        $id .= 'cstart='. $codon_start;
    }
    $id .= '|';
    $debug && print STDERR "$subn: \$id='$id'\n";

    # Add range in polyprotein. This is simply calculated from existing DNA coordinates.
    # One caveat is, in some cases, the last AA of polyprotein may have 2 coden, need special handling.
    my $s = '';
    $s = &get_AA_loc($feat, $cds);
    $id .= 'AA='. $s .'|';
    $debug && print STDERR "$subn: \$id='$id'\n";

    # add gene symbol
    my $desc = $id;
    my $productName = '';
  if (0) {
    my @tags = ('product');
    for my $tag (@tags) {
        if ($feat->has_tag($tag)) {
           @id = $feat->get_tag_values($tag);
           $productName = $id[0];
           last;
        }
    }
#    $gene_symbol = $productName if (!$gene_symbol && length($productName)<=6);
  }
    if ( 0 ) {
        $desc .= 'gene_symbol='. $gene_symbol .'|';
    } else {
        $desc .= 'symbol='. $gene_symbol .'|';
    }
    $debug && print STDERR "$subn: \$desc='$desc'\n";

    # error code: 0: none; 1: gaps around cleavage site; 2: partial mat_peptide; 3: outside of CDS
    if (exists($errcode->{GapAtCleavage}) && $errcode->{GapAtCleavage}==1) {
        $desc = $desc. 'GapAtCleavage=Y|'; # This could be removed in fixCleaveGap
    }
    if (exists($errcode->{OutsideCDS}) && $errcode->{OutsideCDS}==1) {
       $desc = $desc. 'OutsideCDS=Y|';
    }
    if (exists($errcode->{partial_mat_peptide}) && $errcode->{partial_mat_peptide}==1) {
       $desc = $desc. 'Partial=Y|';
    } else {
       $desc = $desc. 'Partial=N|'; # fixing the field --10/12/2012
    }
    $debug && print STDERR "$subn: \$desc='$desc'\n";

    # add the product of mat_peptide
  if (1) {
    my $productName = '';
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
  } else {
    $desc .= "product=$productName";
    $id .= "product=$productName";
  }
    $debug && print STDERR "$subn: \$desc='$desc'\n";

    return ($id, $desc);
} # sub get_feature_id_desc


=head2 get_DNA_loc
Takes $feat
 Returns the location range in a string
 This basically is $feat->location->to_FTstring, but we are turning off the <> notation
=cut

sub get_DNA_loc {
    my ($feat) = @_;

    my $debug = 0 || $debug_all;
    my $subn = 'Annotate_Util::get_DNA_loc';

    # get the range of mat_peptide
    my $location_allow_split = 1;
    my $location_allow_fuzzy = 0;
    my $s;
    if ($location_allow_split) {
        $s = $feat->location->to_FTstring;
        $debug && print STDERR "$subn: may contain '()' in \$s=$s.\n";
    } else {
        $s = $feat->location->start .'..'. $feat->location->end;
        $debug && print STDERR "$subn: no '()' in \$s=$s.\n";
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

    my $debug = 0 || $debug_all;
    my $subn = 'Annotate_Util::get_AA_loc';

    # Add range in polyprotein. This is simply calculated by comparing the sequence of mat_peptide with
    # that of the CDS, and finding the index.
    $debug && print STDERR "$subn: Entering sub: \$feat=".$feat->location->to_FTstring." \$cds=".$cds->location->to_FTstring."\n";
    my $s = '0..0';
    if ( 1 ) {
        my $s_temp = -1;
        my $s1 = Annotate_Util::get_new_translation( $feat, $cds);
        $s1 = '' if (!$s1);
        $s1 =~ s/[.BJZ]/X/g; # Change all .BJZ to X
        my $s0 = '';
        if ( 0 ) {
            $s0 = [ $cds->get_tag_values('translation') ]->[0];
        } else {
            $s0 = Annotate_Util::get_new_translation( $cds, $cds);
        }
        $s0 =~ s/[.BJZ]/X/g; # Change all .BJZ to X
        if ( $cds->strand!=-1 ) {
            $s_temp = $feat->location->start - $cds->location->start;
        } else {
            $s_temp = $cds->location->end - $feat->location->end;
        }
        $debug && print STDERR "$subn: \$s_temp=$s_temp\n";
        $s_temp = $s_temp / 3 - 10;
        $debug && print STDERR "$subn: \$s_temp=$s_temp\n";

        $s_temp = index($s0, $s1, $s_temp);
        $debug && print STDERR "$subn: \$s_temp=$s_temp\n";
        $s_temp += 1;
        $debug && print STDERR "$subn: \$s_temp=$s_temp \$s1='$s1'\n";
        $debug && print STDERR "$subn: \$s_temp=$s_temp \$s0='$s0'\n";
        $s = sprintf("%d..%d", $s_temp, $s_temp+length($s1)-1) if ($s_temp>0);
        $debug && print STDERR "$subn: \$s=$s \$s1=$s1\n";
        $debug && print STDERR "$subn: \$s=$s \$s0=$s0\n";
    } else {
        my $codon_start_feat = 0;
        $codon_start_feat = [$feat->get_tag_values('codon_start')]->[0] -1 if ($feat->has_tag('codon_start'));
        my $codon_start_cds = 0;
        $codon_start_cds = [$cds->get_tag_values('codon_start')]->[0] -1 if ($cds->has_tag('codon_start'));
        my $gap1 = 0;
        my $gap2 = 0;
        if ( $cds->location->isa('Bio::Location::Split') ) {
#            print "location=". $feat_obj->location->to_FTstring . "\n";
            my @locs = ($cds->location->sub_Location);
            for my $iloc (0 .. $#locs-1) {
                my $loc = $locs[$iloc];
#                $debug && print STDERR "$subn: \$loc=\n".Dumper($loc)."End of \$loc\n\n";
#                $debug && print STDERR "$subn: \$loc->strand()=".$loc->strand()."\n";
#                $debug && print STDERR "$subn: \$loc->strand()=".$feat_obj->location->strand()."\n";
                if ($locs[$iloc+1]->start < $feat->location->start ) {
                    $gap1 += $locs[$iloc+1]->start-1 - $loc->end;
                    $gap2 += $locs[$iloc+1]->start-1 - $loc->end;
                } elsif ($locs[$iloc+1]->end < $feat->location->end) {
                    $gap2 += $locs[$iloc+1]->start-1 - $loc->end;
                }
            }
        }
        $debug && print STDERR "$subn: \$gap1=$gap1 \$gap2=$gap2 after checking CDS\n";

        my $b0 = 0;
        if ($feat->location->start==$cds->location->start) {
            $b0 = $feat->location->start +$codon_start_feat - $cds->location->start -$codon_start_cds -$gap1;
        } else {
            $b0 = $feat->location->start                    - $cds->location->start -$codon_start_cds -$gap1;
        }
        $debug && print STDERR "$subn: \$b0=$b0 \$b0/3=". $b0 % 3 ."\n";
        $b0 = (($b0 % 3) ==0) ? $b0 /3 +1 : 0;
        $debug && print STDERR "$subn: \$b0=$b0 after converting from DNA to AA\n";

        my $e0 = 0;
        if ($feat->location->start==$cds->location->start) {
            #$e0 = $feat->location->end +$codon_start_feat + 1 - $cds->location->start -$codon_start_cds;
            $e0 = $feat->location->end                    + 1 - $cds->location->start -$codon_start_cds -$gap2;
        } else {
            $e0 = $feat->location->end                    + 1 - $cds->location->start -$codon_start_cds -$gap2;
        }
        $debug && print STDERR "$subn: \$e0=$e0\n";
        $e0 +=1 if (($e0 +1) % 3 ==0); # some of the polyproteins have 2 nucleotides for the last codon
        $debug && print STDERR "$subn: \$e0=$e0\n";
        #$e0 = ($e0 % 3 ==0) ? $e0 /3 : 0;
        $e0 = (($e0-1) % 3 ==0) ? $e0-1 : $e0;
        $e0 = $e0 /3;
        $debug && print STDERR "$subn: \$e0=$e0 after converting from DNA to AA\n";
        $debug && print STDERR "$subn: \$b0=$b0 \$e0=$e0  after converting from DNA to AA\n";
        $s = $e0 - $b0 +1; # This is the length calculated from the DNA location

        # For split locations, any gap between sublocations needs to be accounted for
        if ($feat->location->isa('Bio::Location::Split')) {
            my $locs = [ $feat->location->sub_Location ];
            my $gap = 0;
            for my $i (1 .. $#{$locs}) {
               $gap = $locs->[$i]->start - $locs->[$i-1]->end -1; 
               $debug && print STDERR "$subn: \$gap=$gap\n";
            }
            $debug && print STDERR "$subn: \$s=$s \$gap=$gap\n";
            $s -= $gap/3 if (($gap%3)==0);
        }

        my $b = $b0;
        my $e = $e0;
        $debug && print STDERR "$subn: \$feat=\n".Dumper($feat)."End of \$feat\n\n";
        # check the length of existing translation against length from coordinates, just in case
        if ($feat->has_tag('translation')) {
            my $t = [ $feat->get_tag_values('translation') ];
            $t = $t->[0];
            $debug && print STDERR "$subn: \$seq=$t\n";
            if ($s!=length($t)) {
                print STDERR  "$subn: ".$cds->seq->accession_number ." has problem between existing translation and coordinate \$s=$s length=".length($t)."\n";
                print STDERR "$subn: \$s=$s \$t=".length($t)."\n";
                print STDERR "$subn: \$feat=\n".Dumper($feat)."End of \$feat\n\n";
                print STDERR "$subn: \$cds=\n".Dumper($cds)."End of \$cds\n\n";
                croak "$subn: ".$cds->seq->accession_number ." has problem between existing translation and coordinate \$s=$s length=".length($t)."\n";
            }
        }
        $s = $b .'..'. $e;
    }
    $debug && print STDERR "$subn: \$s=$s\n";

    return $s;
} # sub get_AA_loc


=head2 project_matpept
Takes references to refset, $inset, alignment, ids of refcds and cds, and a note
 returns the list of features of CDS, mat_peptide, sig_peptide.
=cut

sub project_matpept {
    my ($refset, $inset, $aln, $refcds_id, $cds_id, $note,$exe_dir) = @_;

    my $debug = 0 || $debug_all;
    my $subn = 'project_matpept';

    my $feats_all = []; # Holds all new features for the target genome
#    $debug && print STDERR "$subn: \$refset=\n".Dumper($refset)."End of \$refset\n";
    # Go through all features in refseq, map the corresponding features in targe genome,
    # first, check CDS that is labeled as polyprotein, then all mat_pepride (and sig_peptide) after such CDS
    
    $debug && print STDERR "$subn: \$refset=\n".Dumper($refset)."End of \$refset\n\n";
    for (my $ct = 1; $ct <= $#{$refset}; $ct++) {
        # First get CDS from refseq & target
        my $reffeat = $refset->[$ct];
        $debug && print STDERR "$subn: #$ct is ".$reffeat->primary_tag." \t$reffeat\n";
        if ($reffeat->primary_tag ne "CDS") {
            $debug && print STDERR "$subn: ERROR: refset($ct) is not CDS\n";
            next;
        }
        my @prod = $reffeat->get_tag_values('product');
#        $debug && print STDERR "project_matpept: ct=$ct \$reffeat is ".$reffeat->primary_tag." prod=$prod[0]\n";
        # Skip any CDS not labeled as "polyprotein". Potential problem as some polyproteins are not so labeled
#        next if ($prod[0] !~ /polyprotein/i);
        next if (!Annotate_Def::is_polyprotein( $reffeat, ['product', 'note'], $refset->[$ct+1]));
        $debug && print STDERR "$subn: \$ct=$ct \$reffeat=\n".Dumper($reffeat)."End of \$reffeat\n\n";

        my $refcds = $reffeat;
        my $cds = $inset->[0];
        $debug && print STDERR "$subn: \$cds=".$cds->seq->accession_number."\n";
        my %allowed_feats = ('mat_peptide' => 1, 'sig_peptide' => 1);
        while (($refset->[$ct+1]) && ($allowed_feats{$refset->[$ct+1]->primary_tag})) {
            $reffeat = $refset->[++$ct];
#            $debug && print STDERR "$subn: \$reffeat #$ct is ".$reffeat->primary_tag."=\n".Dumper($reffeat)."End of \$reffeat\n\n";

            $debug && print STDERR "$subn: \$refcds_id=$refcds_id \$cds_id=$cds_id reffeat=".$reffeat->location->to_FTstring."\n";

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
                $debug && print STDERR "$subn: #$ct assemble_new_feature returned \$feat undef, skip\n";
                next;
            }
            $debug && print STDERR "$subn: #$ct \$feat=\n".Dumper($feat)."End of \$feat\n\n";

            # See if this is entirely new annotation wrt genbank
            my $new = 1;
            my $feat_gbk;
            ($new, $feat_gbk) = Annotate_Util::is_new_annotation( $feat, $inset);
            if ($new) {
                # add *new* to the description to indicate that GBK doesn't have this annotation
                my @tags = $feat->get_tag_values('note');
                $feat->remove_tag('note');
                for (my $i = 0; $i<=$#tags; $i++) {
                    $debug && print STDERR "$subn: \$tags[$i] = $tags[$i]\n";
                    if ($tags[$i] =~ /^Desc:/i) {
                        $tags[$i] = $tags[$i] .'|*new*';
                        $feat->add_tag_value('note', $tags[$i]);
                        $debug && print STDERR "$subn: \$tags[$i] = $tags[$i]\n";
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

            # Check if new $feat has same location with existing $feats. If so, omit the one with '=Y'
            my $seen = 0;
            for my $feat_oldi (0..$#{$feats_all}) {
                next if ($feat_oldi>$#{$feats_all}); # In case $feats_all has been changed inside this block
                my $feat_old = $feats_all->[$feat_oldi];

                # Check if the gene_symbol/product is same, skip if not
                # If gene/symbod/product is same and the length are different, take the long feature
                my $symbol = '';
                my $symbol2 = '';
                my @tags = $feat->get_tag_values('note');
                for (my $ii = 0; $ii<=$#tags; $ii++) {
                    $debug && print STDERR "$subn: \$tags[$ii] = $tags[$ii]\n";
                    $symbol = $1 if ($tags[$ii] =~ /symbol=([^|]+[|]product=[^|]+)([|]|$)/i);
                    $symbol = $1 if ($tags[$ii] =~ /symbol=([^|]+[|])/i);
                }
                @tags = $feat_old->get_tag_values('note');
                for (my $ii = 0; $ii<=$#tags; $ii++) {
                    $debug && print STDERR "$subn: \$tags[$ii] = $tags[$ii]\n";
                    $symbol2 = $1 if ($tags[$ii] =~ /symbol=([^|]+[|]product=[^|]+)([|]|$)/i);
                    $symbol2 = $1 if ($tags[$ii] =~ /symbol=([^|]+[|])/i);
                }
                $debug && print STDERR "$subn: \$symbol ='$symbol' \$symbol2='$symbol2'\n";
                my $feat_loc = $feat->location->to_FTstring;
                my $feat_old_loc = $feat_old->location->to_FTstring;
                if ($symbol && $symbol eq $symbol2) {
                    $seen = 1;
                    print STDERR "$subn: Found duplicate based on symbol: $symbol:$feat_old_loc vs $symbol2:$feat_old_loc\n";
                } elsif ($feat_loc eq $feat_old_loc) {
                    $seen = 1;
                    print STDERR "$subn: Found duplicate based on location: $symbol:$feat_old_loc vs $symbol2:$feat_old_loc\n";
                }
                next if (!$seen);

                my $length = $feat->location->end - $feat->location->start +1;
                my $length2 = $feat_old->location->end - $feat_old->location->start +1;
                next if ($length < $length2); # $seen=1 if new feature is shorter

                my $tags = [ $feat_old->get_tag_values('note') ];
                $debug && print STDERR "$subn: \$tags=\n".Dumper($tags)."End of \$tags\n\n";
                for my $tag (@$tags) {
                    next if ($tag !~ /^Desc:/i);
                    $tags = $tag;
                }
                $debug && print STDERR "$subn: \$tags=$tags\n";
                if ($tags =~ /=Y/i) {
                    my @ff = @$feats_all;
                    my @s = ();
                    push @s, @ff[0 .. $feat_oldi-1] if ($feat_oldi>0);
#                    push @s, $feat;  # Include the new feature
                    push @s, @ff[$feat_oldi+1 .. $#{$feats_all}] if ($feat_oldi<$#{$feats_all});
                    $debug && print STDERR "$subn: \@s=\n".Dumper(@s)."End of \@s\n\n";
                    print STDERR "$subn: Found duplicate $symbol with '=Y' in mat_peptide $feat_oldi:$feat_old_loc, removed\n";
                    $feats_all = [ @s ];
                    $seen = 0;
                    last;
                } else {
                    $debug && print STDERR "$subn: Don't see '=Y[|]\n";
                }
                last;
            } # for my $feat_oldi (0..$#{$feats_all})
            push @$feats_all, $feat if (!$seen);
            $debug && print STDERR "$subn: \$feat=\n".Dumper($feat)."End of \$feat\n\n";
        } # while (($refset->[$ct+1]) && ($allowed_feats{$refset->[$ct+1]->primary_tag}))
        unshift @$feats_all, $cds; # prepend CDS to the array
        $debug && print STDERR "$subn: result \$feats_all=\n".Dumper($feats_all)."End of \$feats_all\n\n";

        # check if there is any gaps in new annotation, any gaps is shown as >=?=<
        my ($has_gap, $str) = Annotate_Verify::check_ranges( $feats_all);
        $debug && print STDERR "$subn: \$has_gap=$has_gap \$str=$str\n";
        if ( 1 && $has_gap ) {
#            my $feats = Annotate_Util::fixCleaveGap( $feats_all, $refset, $inset, $aln, $refcds_id, $cds_id);
            my $changed = Annotate_Util::fixCleaveGap( $feats_all, $refset, $inset, $aln, $refcds_id, $cds_id);
            ($has_gap, $str) = Annotate_Verify::check_ranges( $feats_all) if ($changed);
            $debug && print STDERR "$subn: \$has_gap=$has_gap \$str=$str\n";
        }

        # check if there is any partial mat_peptide in the middle of polyprotein
        ($has_gap, $str) = Annotate_Verify::check_partial( $feats_all);
        $debug && print STDERR "$subn: \$has_gap = $has_gap \$str=$str\n";
        if ($has_gap) {
#            my $feats = Annotate_Util::fix_partial( $feats_all, $refset, $inset,);
#            ($has_gap, $str) = Annotate_Verify::check_partial( $feats_all);
#            $debug && print STDERR "$subn: \$has_gap = $has_gap \$str=$str\n";
        }

    } # while ($ct <= $#{$reffeats})

    $debug && print STDERR "$subn: finishing \$feats_all=\n".Dumper($feats_all)."End of \$feats_all\n\n";

    return $feats_all;
} # sub project_matpept


=head2 addNoteGap
  Takes a feature, adds 'GapAtCleavage' field to its note
=cut

sub addNoteGap {
#    my ($feats_all, $refset, $gaps_q, $gaps_h, $note) = @_;
    my ($feat1) = @_;

    my $debug = 0 || $debug_all;
    my $subn = 'addNoteGap';

    my @values = $feat1->get_tag_values('note');
    $debug && print STDERR "$subn: \@values=\n".Dumper(@values)."End of \@values\n";
    for (my $k=0; $k<=$#values; $k++) {
        my $value = $values[$k];
        next if ($value !~ /^Desc:(.+)$/i);
        my $pattn = 'symbol=[a-zA-Z0-9]+\|[^G]';
        next if ($value !~ /$pattn/i);
        $pattn = 'symbol=[a-zA-Z0-9]+';
        if ($value =~ s/($pattn)/$1\|GapAtCleavage=Y/i) {
            $values[$k] = $value;
        }
    }
    $debug && print STDERR "$subn: \@values=\n".Dumper(@values)."End of \@values\n";
    $feat1->remove_tag('note');
    $feat1->add_tag_value('note', @values);

} # sub addNoteGap


=head2 removeNoteGap
  Takes a feature, removes the 'GapAtCleavage' field from its note
=cut

sub removeNoteGap {
#    my ($feats_all, $refset, $gaps_q, $gaps_h, $note) = @_;
    my ($feat2) = @_;

    my $debug = 0 || $debug_all;
    my $subn = 'removeNoteGap';


    # After the modification, update the dscription
    my @values = $feat2->get_tag_values('note');
    $debug && print STDERR "$subn: \@values=\n".Dumper(@values)."End of \@values\n";
    for (my $k=0; $k<=$#values; $k++) {
        my $value = $values[$k];
        next if ($value !~ /^Desc:(.+)$/i);
        my $pattn = 'GapAtCleavage=Y\|';
        if ($value =~ s/($pattn)//i) {
            $values[$k] = $value;
            $debug && print STDERR "$subn: \$value = '$value'\n";
        }
    }
    $debug && print STDERR "$subn: \@values=\n".Dumper(@values)."End of \@values\n";
    $feat2->remove_tag('note');
    $feat2->add_tag_value('note', @values);

} # sub removeNoteGap


=head2 removeNoteNew
  Takes a feature, removes the '*new*' field from its note
=cut

sub removeNoteNew {
#    my ($feats_all, $refset, $gaps_q, $gaps_h, $note) = @_;
    my ($feat2, $inset) = @_;

    my $debug = 0 || $debug_all;
    my $subn = 'removeNoteNew';

    # After the modification, update the dscription
    my @values = $feat2->get_tag_values('note');
    $debug && print STDERR "$subn: \@values=\n".Dumper(@values)."End of \@values\n";
    for (my $k=0; $k<=$#values; $k++) {
        my $value = $values[$k];
        next if ($value !~ /^Desc:(.+)$/i);
        my ($new, $feat_gbk) = &is_new_annotation( $feat2, $inset);
        $debug && print STDERR "$subn: \$new=$new\n";
        my $pattn = '\|[*]new[*]';
        if (!$new && $value =~ s/($pattn)//i) {
            $values[$k] = $value;
            $debug && print STDERR "$subn: \$value = '$value'\n";
        }
    }
    $debug && print STDERR "$subn: \@values=\n".Dumper(@values)."End of \@values\n";
    $feat2->remove_tag('note');
    $feat2->add_tag_value('note', @values);

} # sub removeNoteNew


=head2 getGapSeqs
  Takes 2 mat_peptides and their CDS, returns the last $MAX_LENGTH of the 1st mat_peptide, 
 the first $MAX_LENGTH of the 2nd mat_peptide, and the gap in between
=cut

sub getGapSeqs {
    my ($feat1, $feat2, $MAX_LENGTH, $aln, $refcds_id, $cds_id) = @_;

    my $debug = 0 || $debug_all;
    my $subn = 'getGapSeqs';

    # For CDS, find the C-term of 1st mat_peptide, gap, and N-term of 2nd mat_peptide
    my $aln_q = Annotate_Util::msa_get_aln( $aln, $refcds_id);
    my $aln_h = Annotate_Util::msa_get_aln( $aln, $cds_id);
    my ($str, $l, $n, $f, $gapstart, $gapend);
    my $ts = ['', '', ''];
    my $rs = ['', '', ''];
    # Get the last 10AA of feat1
    my @values = $feat1->get_tag_values('note');
    for (my $k=0; $k<=$#values; $k++) {
        next if ($values[$k] !~ /^MSA=(\d+)\.\.(\d+)/i);
        $l = $2;
        $n = ($l>=$MAX_LENGTH) ? $MAX_LENGTH : $l;
        $rs->[0] = substr($aln_q->seq, $l-$n, $n);
        $ts->[0] = substr($aln_h->seq, $l-$n, $n);
        $gapstart = $l+1;
        last;
    }
    # Get the first 10AA of feat2
    @values = $feat2->get_tag_values('note');
    for (my $k=0; $k<=$#values; $k++) {
        next if ($values[$k] !~ /^MSA=(\d+)\.\.(\d+)/i);
        $l = $1;
        $n = ($l>=$MAX_LENGTH) ? $MAX_LENGTH : $l;
        $rs->[2] = substr($aln_q->seq, $l-1, $n);
        $ts->[2] = substr($aln_h->seq, $l-1, $n);
        $gapend = $l-1-1;
        last;
    }
    # Get the gap sequence in feat, and corresponding seq in reffeat
    $debug && print STDERR "$subn: \$gapstart=$gapstart \$gapend=$gapend\n";
    $debug && print STDERR "$subn: \$feat1=\n".Dumper($feat1)."End of \$feat1\n";
    $debug && print STDERR "$subn: \$feat2=\n".Dumper($feat2)."End of \$feat2\n";
    my $q = substr($aln_q->seq, $gapstart-1, $gapend-$gapstart+2);
    my $h = substr($aln_h->seq, $gapstart-1, $gapend-$gapstart+2);
    my $c;
    my $gap_char = $aln->gap_char;
    for my $i (1..length($q)) {
        $c = substr($h, $i-1, 1);
        next if ($c eq $gap_char);
        $rs->[1] .= substr($q, $i-1, 1);
        $ts->[1] .= $c;
    }

    $debug && print STDERR "$subn: \@\$rs='@$rs'\n";
    $debug && print STDERR "$subn: \@\$ts='@$ts'\n";
    return ($ts, $rs);
} # sub getGapSeqs


=head2 fixCleaveGap

Takes the new annotation in an array, refcds, cds, and the gaps.
 Check if all new mat_peptides are connected
 if not, decide if the missing sequence at the cleavage site should go to the tail of previous mat_peptide,
 or head of the next mat_peptide

If the alignment is not definite, as in EF407458 EF407463 EF407467 when run together and separately, the annotation will change according to the MSA. This is a consequence of the instability in the alignment.

=cut

sub fixCleaveGap {
    my ($feats_all, $refset, $inset, $aln, $refcds_id, $cds_id) = @_;

    my $debug = 0 || $debug_all;
    my $subn = 'fixCleaveGap';

    my $cds = $feats_all->[0];
#    my $refcds = $refset->[0];
    my $refcds = $refset->[1]; # #0 is id
    my $MAX_LENGTH = 10;
    my $has_gap = 0;
    my $feats;
    my $changed = 0;
    $debug && print STDERR "$subn: \$feats_all=\n".Dumper($feats_all)."end of \$feats_all\n";
    $debug && print STDERR "$subn: \$refset=\n".Dumper($refset)."end of \$refset\n";

  # This for loop is used to take care of cases where 2 mat_peptides end at same
  # location and there is a gap, eg. preM and M
  for (my $nfeat=1; $nfeat<$#{$feats_all}; $nfeat++) {
    my $seen_gap = 0; # Used to indicate if we need to update $nfeat
    my $feat1_n = $nfeat; # #0 is CDS
    my $acc = $feats_all->[0]->seq->accession_number;
#    for (my $i=2; $i<=$#{$feats_all}; $i++) {
    for (my $i=$nfeat+1; $i<=$#{$feats_all}; $i++) {
        my $feat1 = $feats_all->[$feat1_n];
        my $feat2 = $feats_all->[$i];
        $debug && print STDERR "$subn: \$feat1_n=$feat1_n \$i=$i\n";

        my $loc1 = $feat1->location;
        my $loc2 = $feat2->location;
        # Check the 1st end vs 2nd start
        if ($loc1->end+1 == $loc2->start) {
            # Skip if there is no gap at cleavage site
            $feat1_n = $i;
            $nfeat = $i if (!$seen_gap);
            $debug && print STDERR "$subn: ACC=$acc \$nfeat=$nfeat \$i=$i \$seen_gap=$seen_gap\n";
            next;
        } elsif ($loc1->end > $loc2->start) {
            # Skip if $feat2 starts before the end of $feat1, mostly indicating $feat2 is a part of $feat1
            $seen_gap = 1;
            $debug && print STDERR "$subn: ACC=$acc \$nfeat=$nfeat \$i=$i \$seen_gap=$seen_gap\n";
            next;
        }
        print STDERR "$subn: ACC=$acc Found a gap before #$i: between ".$loc1->end ."..". $loc2->start."\n";

        $seen_gap = 1;
        $has_gap = 1;
        print STDERR "$subn: ACC=$acc \$feat1_n=$feat1_n \$i=$i, found a gap\n";
        $debug && print STDERR "$subn: \$feat1=\n".Dumper($feat1)."end of \$feat1\n";
        $debug && print STDERR "$subn: \$feat2=\n".Dumper($feat2)."end of \$feat2\n";

        # For CDS, find the C-term of 1st mat_peptide, gap, and N-term of 2nd mat_peptide
        my ($str, $l, $n, $f);
        my $ts = ['', '', ''];
        my $rs = ['', '', ''];
        ($ts, $rs) = Annotate_Util::getGapSeqs( $feat1, $feat2, $MAX_LENGTH, $aln, $refcds_id, $cds_id);

        # Determine the reffeat for each feature
        my ($reffeat1, $reffeat2);
        if ( 1 ) {
            my ($str);
            # Determine the reffeat for feat1
            my @values;
            @values = $feat1->get_tag_values('note');
            for (my $k=0; $k<=$#values; $k++) {
#                $str = $1 if ($values[$k] =~ /reffeat=(\d+\.\.\d+)/i);
                $str = $1 if ($values[$k] =~ /RM=(\d+\.\.\d+)/i);
            }
            for my $ref (1 .. $#{$refset}) {
                $ref = $refset->[$ref];
                my $s = $ref->location->start. '..' .$ref->location->end;
                $debug && print STDERR "$subn: \$str=$str \$s=$s\n";
                next if ($s ne $str);
                $reffeat1 = $ref;
                last;
            }

            # Determine the reffeat for feat2
            @values = $feat2->get_tag_values('note');
            for (my $k=0; $k<=$#values; $k++) {
#                $str = $1 if ($values[$k] =~ /reffeat=(\d+\.\.\d+)/i);
                $str = $1 if ($values[$k] =~ /RM=(\d+\.\.\d+)/i);
            }
            for my $ref (1 .. $#{$refset}) {
                $ref = $refset->[$ref];
                my $s = $ref->location->start. '..' .$ref->location->end;
                $debug && print STDERR "$subn: \$str=$str \$s=$s\n";
                next if ($s ne $str);
                $reffeat2 = $ref;
                last;
            }
        }
        $debug && print STDERR "$subn: \$reffeat1=\n".Dumper($reffeat1)."end of \$reffeat1\n";
        $debug && print STDERR "$subn: \$reffeat2=\n".Dumper($reffeat2)."end of \$reffeat2\n";

        my $rlen = 0;
        if ($reffeat2->strand==-1) {
            $rlen = ($reffeat1->location->start - $reffeat2->location->end -1)/3;
        } else {
            $rlen = ($reffeat2->location->start - $reffeat1->location->end -1)/3;
        }
        $debug && print STDERR "$subn: \$rs=".length($rs->[1])." \$rlen=$rlen\n";
        if ( $rlen && abs(length($rs->[1]) - $rlen) /$rlen < 0.2 ) {
            # Skip if there is gap in the RefSeq, such as ACC=AB127995|ref=NC_005219.1
            print STDERR "$subn: ACC=$acc \$rs='@$rs'\n";
            print STDERR "$subn: ACC=$acc \$ts='@$ts'\n";
            print STDERR "$subn: ACC=$acc gap b/w ".$loc1->to_FTstring ." and ". $loc2->to_FTstring." is caused by gap in RefSeq, skip\n";
            $feat1_n = $i;
            next;
        } elsif (!$rs->[1]) {
            $rs->[1] = $ts->[1];
            $rs->[1] =~ s/./-/g; # pad with '-'
        }
        print STDERR "$subn: ACC=$acc \$rs='@$rs'\n";
        print STDERR "$subn: ACC=$acc \$ts='@$ts'\n";

        # Update the dscription to include GapAtCleavage
        Annotate_Util::addNoteGap( $feat1);
        Annotate_Util::addNoteGap( $feat2);

        # The gap would be lumped into the beginning of next mat_peptide provided that
        # 1) there are min. of 5 identical residues in the last 10 residues of the previous mat_peptide;
        # 2) the last residue in refseq and target are identical or from same class: (AG) (ST) (RK) (C);
        # 3) the gap is 8 or less residues long. -Ver1.1.5
        # New decision tree: --Ver1.1.6, Oct 24, 2012
        # 1. If there is still a gap after atempt to fix it, give null result for the genome
        # 2. A gap has to be 8AA or shorter. (If RefSeq has a gap between mat_peptide, no problem)
        # 3. Within the nearest 5AA before and after the gap. Whichever has fewer similar AA would get
        # the gap sequence. Similar residues are defined as AG, KR, NQ, DE, ILV, FY, ST
        # 4. If both sides have same similar residues, the gap goes to 2nd mat_peptide.
        # 5. If neither side has 2 or more similar residues, the gap goes to 2nd mat_peptide.
        my $max_gap_length = 8;
        my $min_identical_within10 = 2;
        my $min_similar_within5 = 2;

        my $length_gap = length($ts->[1]);

        my $identical_last10 = 0;
        my $similar_last5 = 0;
        my $id_last_residue = 0;

        my $identical_first10 = 0;
        my $similar_first5 = 0;
        my $id_first_residue = 0;

        $n = (length($rs->[0])<=length($ts->[0])) ? length($rs->[0]) : length($ts->[0]);
        my ($rc, $tc);
        foreach my $j (reverse 1 .. $n) {
            $rc = uc substr($rs->[0], length($rs->[0])-$j, 1);
            $tc = uc substr($ts->[0], length($ts->[0])-$j, 1);
            $identical_last10++ if ($rc eq $tc);
            $similar_last5++ if ($j<=5 && Annotate_Def::similarAA($rc, $tc));
            $debug && print STDERR "$subn: \$j=$j \$rc=$rc \$tc=$tc \$identical_last10=$identical_last10 \$similar_last5=$similar_last5\n";
        }
        $rc = uc substr($rs->[0], length($rs->[0])-1, 1);
        $tc = uc substr($ts->[0], length($ts->[0])-1, 1);
        $id_last_residue = 1 if (Annotate_Def::similarAA($rc, $tc));

        $n = (length($rs->[2])<=length($ts->[2])) ? length($rs->[2]) : length($ts->[2]);
        foreach my $j (0 .. $n-1) {
            $rc = uc substr($rs->[2], $j, 1);
            $tc = uc substr($ts->[2], $j, 1);
            $identical_first10++ if ($rc eq $tc);
            $similar_first5++ if ($j<5 && Annotate_Def::similarAA($rc, $tc));
            $debug && print STDERR "$subn: \$j=$j \$rc=$rc \$tc=$tc \$identical_first10=$identical_first10 \$similar_first5=$similar_first5\n";
        }
        $rc = uc substr($rs->[2], 0, 1);
        $tc = uc substr($ts->[2], 0, 1);
        $id_first_residue = 1 if (Annotate_Def::similarAA($rc, $tc));

        my $msg;
    if ( 1 ) {
        $msg =  "$subn: ERROR: ACC=$acc There is a gap between mat_peptides #$feat1_n and #$i: '@$ts'\n";
#        $msg .= "$subn: ACC=$acc mat_peptide #$feat1_n has $identical_last10 identical AA in last 10, #$i has $identical_first10 in first 10\n";
        $msg .= "$subn: ACC=$acc mat_peptide #$feat1_n has $similar_last5 similar AA in last 5, #$i has $similar_first5 in first 5\n";
        $msg .= ($length_gap>$max_gap_length) ? "$subn: 'X': " : "$subn: 'V': ";
        $msg .= "Length of gap between #$feat1_n-$i is $length_gap, ";
        $msg .= ($length_gap>$max_gap_length) ? ">$max_gap_length(required).\n" : "<=$max_gap_length(required).\n";

        if ($length_gap>$max_gap_length) {
            $msg .= "$subn: ACC=$acc mat_peptide #$feat1_n-#$i gap is too long. Skip the genome ACC=$acc\n";
            print STDERR $msg;
            return undef;
        } elsif ($similar_last5<$similar_first5 && $min_similar_within5<=$similar_first5) {
            $msg .= ($similar_first5<$min_similar_within5) ? "$subn: 'X': " : "$subn: 'V': ";
            $msg .= "$similar_first5 of first 5 AAs are similar in mat_peptide #$i, ";
            $msg .= "required is $min_similar_within5 or more\n";

            my $msg1 = Annotate_Util::lump_gap_to_feat1( $feat1, $feat2, $cds, $length_gap, $ts, $feat1_n, $i, $inset);
            $changed = 1;
            if (!$msg) {
                $msg .= "$subn: lump_gap_to_feat1 returned undef. Skip the genome ACC=$acc\n";
                print STDERR $msg;
                return undef;
            } else {
                $msg .= "$subn: Gap has been attached to the end of 1st mat_peptide\n";
            }
            $msg .= $msg1;
        } else {
            $msg .= ($similar_last5<$min_similar_within5) ? "$subn: 'X': " : "$subn: 'V': ";
            $msg .= "$similar_last5 of last 5 AAs are similar in mat_peptide #$feat1_n, ";
            $msg .= "required is $min_similar_within5 or more\n";

            my $msg1 = Annotate_Util::lump_gap_to_feat2( $feat1, $feat2, $cds, $length_gap, $ts, $feat1_n, $i, $inset);
            $changed = 1;
            if (!$msg) {
                $msg .= "$subn: lump_gap_to_feat1 returned undef. Skip the genome ACC=$acc\n";
                print STDERR $msg;
                return undef;
            } else {
                $msg .= "$subn: Gap has been attached to the start of 2nd mat_peptide\n";
            }
            $msg .= $msg1;
        }
        print STDERR $msg;
#        print STDOUT $msg;

        # After the modification, update the dscription
        if ($feat1->location->end+1 == $feat2->location->start) {
            Annotate_Util::removeNoteGap($feat1);
            Annotate_Util::removeNoteNew($feat1, $inset);
            Annotate_Util::removeNoteGap($feat2);
            Annotate_Util::removeNoteNew($feat2, $inset);

            $msg =  "$subn: WARNING: the gap ($ts->[1]) between mat_peptides #$feat1_n and #$i:$ts->[1].";
            $msg .= " The gap has been lumped into the next mat_peptide\n";

        }
    }

        $feat1_n = $i; # Keep the current mat_peptide as starting point.
        $debug && print STDERR "$subn: \$feat1_n = $feat1_n\n";
    } # for (my $i=2; $i<=$#{$feats_all}; $i++)
  } # for (my $nfeat=1; $nfeat<$#{$feats_all}; $nfeat++)

    $has_gap && $debug && print STDERR "$subn: \$feats_all=\n".Dumper($feats_all)."End of \$feats_all\n\n";

#    return ($feats);
    return $changed;
} # sub fixCleaveGap


=head2 insolatePttn
sub insolatePttn adds a \ in front of ()
=cut

sub insolatePttn {
    my ($pttn) = @_;

    my $debug = 0 || $debug_all;
    my $subn = 'insolatePttn';

    my $pttn1 = $pttn;
    $pttn1 =~ s/\(/\\\(/g;
    $pttn1 =~ s/\)/\\\)/g;
    $debug && print STDERR "$subn: \$pttn='$pttn' \$pttn1='$pttn1'\n";

    return $pttn1;
} # sub insolatePttn


=head2 lump_gap_to_feat1
sub lump_gap_to_feat1 takes the gap at the cleavage site in target sequence and add it to the mat_peptide before cleavage site.
=cut

sub lump_gap_to_feat1 {
    my ($feat1, $feat2, $cds, $length_gap, $ts, $feat1_n, $i, $inset) = @_;

    my $debug = 0 || $debug_all;
    my $subn = 'lump_gap_to_feat1';

    # Move the start of next mat_peptide by the length of gap
    my $oldDNA = &insolatePttn( 'Loc='. &get_DNA_loc($feat1) );
    my $oldAA = &insolatePttn( 'AA='. &get_AA_loc($feat1, $cds) );
    $debug && print STDERR "$subn: \$oldDNA='$oldDNA' \$oldAA='$oldAA'\n";
    my $loc2 = $feat1->location;
    if ($loc2->isa('Bio::Location::Simple') ) {
#        my $start = $loc2->start - $length_gap*3;
        $loc2->end($loc2->end + $length_gap*3);
    } elsif ($loc2->isa('Bio::Location::Fuzzy')) {
#        my $min_start;
#        $min_start = $loc2->min_start - $length_gap*3;
#        $loc2->min_start($min_start);
        $loc2->min_end($loc2->min_end + $length_gap*3);
#        my $max_start;
#        $max_start = $loc2->max_start - $length_gap*3;
#        $loc2->max_start($max_start);
        $loc2->max_end($loc2->max_end + $length_gap*3);
    }
    $debug && print STDERR "$subn: \$loc2=\n".Dumper($loc2)."end of \$loc2\n";

    # Now, update the translation
    my $translation = '';
    if ($cds) {
        my $s = Annotate_Util::get_new_translation( $feat1, $cds);
        if ($s) {
            $translation = $s;
            $feat1->remove_tag('translation');
            $feat1->add_tag_value('translation', $translation);
        } else {
            # undo the change in start if this fails
            if ($loc2->isa('Bio::Location::Simple') ) {
                $loc2->end($loc2->end - $length_gap*3);
            } elsif ($loc2->isa('Bio::Location::Fuzzy')) {
                $loc2->min_end($loc2->min_end - $length_gap*3);
                $loc2->max_end($loc2->max_end - $length_gap*3);
            }
            print STDERR "$subn: \$feat1=\n".Dumper($feat1)."End of \$feat1\n\n";
            print STDERR "$subn: ERROR: translation for mat_peptide doesn't match CDS.\n";
            print STDERR "$subn: \$s='$s'\n";
            print STDERR "$subn: \$cds='".$cds->seq->translate->seq."'\n";
            return undef;
        }

    }

    # Update the dscription
    my $newDNA = 'Loc='. &get_DNA_loc($feat1);
    my $newAA = 'AA='. &get_AA_loc($feat1, $cds);
    $debug && print STDERR "$subn: \$newDNA='$newDNA' \$newAA='$newAA'\n";
    my @values = $feat1->get_tag_values('note');
    for (my $k=0; $k<=$#values; $k++) {
        my $value = $values[$k];
        next if ($value !~ /^Desc:(.+)$/i);
        # Change the location
 #       if ($value =~ /(Loc=[0-9]+[.]{2,2}[0-9]+)/i) {
        $debug && print STDERR "$subn: \$value = '$value'\n";
        if ($value =~ /$oldDNA/i) {
            $value =~ s/$oldDNA/$newDNA/;
            $values[$k] = $value;
            $debug && print STDERR "$subn: \$value = '$value'\n";
        }
        if ($value =~ /$oldAA/i) {
            $value =~ s/$oldAA/$newAA/;
            $values[$k] = $value;
            $debug && print STDERR "$subn: \$value = '$value'\n";
        }
        my ($new, undef) = &is_new_annotation( $feat1, $inset);
        $debug && print STDERR "$subn: \$new='$new'\n";
        my $pattn = '\|[*]new[*]';
        if (!$new && $value =~ s/($pattn)//i) {
            $values[$k] = $value;
            $debug && print STDERR "$subn: \$value = '$value'\n";
            $debug && print STDERR "$subn: \@values=\n".Dumper(@values)."End of \@values\n\n";
        }
    }
    $feat1->remove_tag('note');
    $feat1->add_tag_value('note', @values);

    my $msg =  "$subn: WARNING: the gap ($ts->[1]) between mat_peptides #$feat1_n and #$i:$ts->[1]. ";
    $msg .= " The gap has been lumped into the first mat_peptide (before cleavage site)\n";

    return $msg;
} # sub lump_gap_to_feat1


=head2 lump_gap_to_feat2
sub lump_gap_to_feat2 takes the gap at the cleavage site in target sequence and add it to the next mat_peptide.
=cut

sub lump_gap_to_feat2 {
    my ($feat1, $feat2, $cds, $length_gap, $ts, $feat1_n, $i, $inset) = @_;

    my $debug = 0 || $debug_all;
    my $subn = 'lump_gap_to_feat2';

    # Move the start of next mat_peptide by the length of gap
    my $oldDNA = &insolatePttn( 'Loc='. &get_DNA_loc($feat2) );
    my $oldAA = &insolatePttn( 'AA='. &get_AA_loc($feat2, $cds) );
    $debug && print STDERR "$subn: \$oldDNA='$oldDNA' \$oldAA='$oldAA'\n";
    my $loc2 = $feat2->location;
    if ($loc2->isa('Bio::Location::Simple') ) {
#        my $start = $loc2->start - $length_gap*3;
        $loc2->start($loc2->start - $length_gap*3);
    } elsif ($loc2->isa('Bio::Location::Fuzzy')) {
#        my $min_start;
#        $min_start = $loc2->min_start - $length_gap*3;
#        $loc2->min_start($min_start);
        $loc2->min_start($loc2->min_start - $length_gap*3);
#        my $max_start;
#        $max_start = $loc2->max_start - $length_gap*3;
#        $loc2->max_start($max_start);
        $loc2->max_start($loc2->max_start - $length_gap*3);
    }
    $debug && print STDERR "$subn: \$loc2=\n".Dumper($loc2)."end of \$loc2\n";

    # Now, update the translation
    my $translation = '';
    if ($cds) {
        my $s = Annotate_Util::get_new_translation( $feat2, $cds);
        if ($s) {
            $translation = $s;
            $feat2->remove_tag('translation');
            $feat2->add_tag_value('translation', $translation);
            $debug && print STDERR "$subn: \$translation='$translation'\n";
        } else {
            # undo the change in start if this fails
            if ($loc2->isa('Bio::Location::Simple') ) {
                $loc2->start($loc2->start + $length_gap*3);
            } elsif ($loc2->isa('Bio::Location::Fuzzy')) {
                $loc2->min_start($loc2->min_start + $length_gap*3);
                $loc2->max_start($loc2->max_start + $length_gap*3);
            }
            print STDERR "$subn: \$feat2=\n".Dumper($feat2)."End of \$feat2\n\n";
            print STDERR "$subn: ERROR: translation for mat_peptide doesn't match CDS.\n";
            print STDERR "$subn: \$s='$s'\n";
            print STDERR "$subn: \$cds='".$cds->seq->translate->seq."'\n";
            return undef;
        }
    }

    # Update the dscription
    my $newDNA = 'Loc='. &get_DNA_loc($feat2);
    my $newAA = 'AA='. &get_AA_loc($feat2, $cds);
    $debug && print STDERR "$subn: \$newDNA='$newDNA' \$newAA='$newAA'\n";
    my @values = $feat2->get_tag_values('note');
    for (my $k=0; $k<=$#values; $k++) {
        my $value = $values[$k];
        next if ($value !~ /^Desc:(.+)$/i);
        # Change the location
#        if ($value =~ /(Loc=[0-9]+[.]{2,2}[0-9]+)/i) {
        $debug && print STDERR "$subn: \$value = '$value'\n";
        if ($value =~ /$oldDNA/i) {
            $value =~ s/$oldDNA/$newDNA/;
            $values[$k] = $value;
            $debug && print STDERR "$subn: \$value = '$value'\n";
        }
        if ($value =~ /$oldAA/i) {
            $value =~ s/$oldAA/$newAA/;
            $values[$k] = $value;
            $debug && print STDERR "$subn: \$value = '$value'\n";
        }
        my ($new, undef) = &is_new_annotation( $feat2, $inset);
        $debug && print STDERR "$subn: \$new='$new'\n";
        my $pattn = '\|[*]new[*]';
        if (!$new && $value =~ s/($pattn)//i) {
            $values[$k] = $value;
            $debug && print STDERR "$subn: \$value = '$value'\n";
            $debug && print STDERR "$subn: \@values=\n".Dumper(@values)."End of \@values\n\n";
        }
    }
    $feat2->remove_tag('note');
    $feat2->add_tag_value('note', @values);

    my $msg =  "$subn: WARNING: the gap ($ts->[1]) between mat_peptides #$feat1_n and #$i:$ts->[1]. ";
    $msg .= " The gap has been lumped into the next mat_peptide\n";

    return $msg;
} # sub lump_gap_to_feat2


=head2 do_bl2seq_search
Takes 2 sequence strings to align against each other.
Runs a pairwise bl2seq query looking for a segment of the second sequence
which matches the first string.
Returns a object of Bio::Search::Result::BlastResult
=cut

sub run_bl2seq_search {
    # Get the shortstring and the query string
    my ($s1string, $s2string, $debug) = @_;

    $debug = $debug || $debug_all;
    my $subn = 'run_bl2seq_search';

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
    my $cmd_bl2seq = "~/prog/blast/blast-2.2.20/bin/bl2seq -p blastp -i $s1temp -j $s2temp -o $blast_outfile -F F ";
#    `bl2seq -p blastp -i $s1temp -j $s2temp -o $blast_outfile -F F`;
#    `$blast_path/bl2seq -p tblastn -i $s1temp -j $s2temp -o $blast_outfile -F F`;
    $cmd_bl2seq = "bl2seq -p blastp -i $s1temp -j $s2temp -o $blast_outfile -F F ";
    ` $cmd_bl2seq `;
    print STDERR "$subn: \$cmd_bl2seq='$cmd_bl2seq'\n";

    if ( 1 && $debug ) {
        open my $debugfile, '<', $blast_outfile;
        print STDERR "\n";
        print STDERR "$subn: >>>>>>>>>\$s1=\n$s1string\n";
        print STDERR "$subn: >>>>>>>>>\$s2=\n$s2string\n";
        print STDERR "$subn: >>>>>>>>>do_bl2seq_search: $blast_outfile=\n";
        print STDERR <$debugfile>;
        print STDERR "$subn: >>>>>>>>>do_bl2seq_search: end of $blast_outfile\n\n";
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

    my $debug = 0 || $debug_all;
    my $subn = 'is_new_annotation';

    $debug && print STDERR "\n$subn: \$feat=\n".Dumper($feat)."End of \$feat\n";
    $debug && print STDERR "\n$subn: \$inset=\n".Dumper($inset)."End of \$inset\n";
    my $new = 1;
    my $feat_gbk;
    $debug && print STDERR "$subn: \$new = $new\n";
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
        $debug && print STDERR "$subn: \$str2='$str2' \$str1='$str1'\n";
        if ($str1 eq $str2) {
            $new = 0;
            $feat_gbk = $gbmatp;
            last;
        }
    }

    return ($new, $feat_gbk);
} # sub is_new_annotation

=head1
sub check_program, takes a program, tests if the program exists or can be executed.
  Could add more tests.
=cut

sub check_program {
    my ($program, $debug) = @_;
    my $subn = "Annotate_Util::check_program";
        
       $debug = 0;
    my $errcode = '';
    my $result = '';
    
    my $command = "which $program";
    $result = `$command`;
    chomp $result;
    $debug && print STDERR "$subn: Checking program with command '$command'\n";
    $debug && print STDERR "$subn: \$result='$result'\n";
    if ($result =~ /which: no /) {
        $errcode = "$subn: Program '$program' not accessible as defined\n";
    } else {
        $debug && print STDERR "$subn: Program '$program' is accessible\n";
    }
    
    return $errcode;
} # sub check_program


1;
