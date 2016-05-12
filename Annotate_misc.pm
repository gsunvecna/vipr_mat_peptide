package Annotate_misc;

use strict;
use warnings;
use English;
use Carp;
use Data::Dumper;

use version; our $VERSION = qv('1.1.2'); # February 09 2011
use File::Temp qw/ tempfile tempdir /;
use Bio::SeqIO;
use Bio::Seq;
use IO::String;

use Annotate_Math;
use Annotate_Verify;
use Annotate_Muscle;
use Annotate_gbk;		# for annotation from genbank

my $debug_all = 1;

####//README//####
#
# Annotate_misc contains the misc subroutine
#
#    Authors Guangyu Sun, gsun@vecna.com; Chris Larsen, clarsen@vecna.com
#    February 2010
#
##################

## //EXECUTE// ##


=head2 generate_fasta

Takes an array of SeqFeature.
 Write the translations to a virtual file
 Returns the fasta file in a string

=cut

sub generate_fasta {
    my ($feats_all) = @_;

    my $debug = 0 && $debug_all;
    my $subname = 'generate_fasta';

    my $seq_out;
    my $faa1 = '';
    my $vfile = IO::String->new($faa1);
    $seq_out = Bio::SeqIO->new(
                  '-fh'     => $vfile,
                  '-format' => 'fasta'
                              );
    # Following sort causes wrong ordering in AF126284
#    foreach my $feat (sort {$a->location->start <=> $b->location->start} @$feats_all) {
    foreach my $feat (@$feats_all) {

        $debug && print STDERR "$subname: \$feat=\n". Dumper($feat) . "End of \$feat\n";
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
            $desc .= '.'.$feat->seq->version; # Add version of accession to description
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


sub process_list1 {
    my ($accs, $aln_fn, $dbh_ref, $exe_dir, $exe_name, $dir_path, $test1, $refseq_required) = @_;

    my $debug = 0 && $debug_all;
    my $subname = 'process_list1';

#    my $refseqs = {};
#    my $inseqs  = [];

   my $msgs = [];
   my $count = {
         MSA=>{Success=>0, Fail=>0, Empty=>0},
         GBK=>{Success=>0, Fail=>0, Empty=>0},
   };
   for (my $ct = 0; $ct<=$#{@$accs}; $ct++) {
        my $msg_msa = '';
        my $msg_gbk = '';
        my $faa = '';
        my $number = $accs->[$ct]->[0];
        my $acc = $accs->[$ct]->[1];

        # get genbank file from MySQL database
        my $result;
        print STDERR "$subname: \$acc='$acc'\n";
        if ($dbh_ref && $dbh_ref->isa('GBKUpdate::Database')) {
            # get genome gbk from MySQL database
            $result = $dbh_ref->get_genbank($acc);
            $result = $result->[0];
        } elsif (-r "$acc") {
            # This $acc is really a filename
            open my $in, '<', "$acc" or croak "$0: Couldn't open $acc: $OS_ERROR";
            $result = do { local $/; <$in>};
            close $in or croak "$0: Couldn't close $acc: $OS_ERROR";
        }
        if (!$result) {
            print STDERR "$subname: #$number:$acc result is empty\n";
#            print STDOUT "$subname: #$number:$acc result is empty\n";
            my $msg = "$acc \tsrc=MSA \tstatus=Fail \tcomment=Empty genome file";
            print STDERR "$subname: \$msg=$msg\n";
            push @$msgs, $msg;
            $msg = "$acc \tsrc=GBK \tstatus=Fail \tcomment=Empty genome file";
            print STDERR "$subname: \$msg=$msg\n";
            push @$msgs, $msg;
#            $count->{MSA}->{Fail}++;
            $count->{MSA}->{Empty}++;
#            $count->{GBK}->{Fail}++;
            $count->{GBK}->{Empty}++;
            next;
        } else {
            print STDERR "$subname: \$result='".substr($result,0,79)."'\n";
        }

        print STDERR "$subname: #$number:\t$acc processing\n";
#        print STDOUT "$subname: #$number:\t$acc processing\n";

        # Now run the annotation by MUSCLE alignment
        my $gbk = $result;
        my $in_file2 = IO::String->new($gbk);
        my $in  = Bio::SeqIO->new( -fh => $in_file2, -format => 'genbank' );
        $acc = $in->next_seq()->accession_number;
        my $outfile = '';
        $outfile = "$dir_path/$acc" . '_matpept_msagbk.faa' if (!$debug);
        print STDERR "$subname: accession='$acc' \$outfile=$outfile.\n";
        
        my ($feats_msa, $comment_msa) = Annotate_Muscle::annotate_1gbk( $gbk, $exe_dir, $aln_fn, $dir_path);
        my $status_msa = $feats_msa ? 'Success' : ($comment_msa eq 'Refseq with mat_peptide annotation from NCBI, skip') ? 'Skip   ' : 'Fail   ';
        if (!$feats_msa) {
            $count->{MSA}->{Fail}++;
            print STDERR "$subname: no result from Annotate_Muscle::annotate_1gbk comment=$comment_msa\n";
        } else {
            $count->{MSA}->{Success}++;
            $debug && print STDERR "$subname: \$feats_msa='\n". Dumper($feats_msa) . "End of \$feats_msa\n\n";
        }
        $msg_msa = "$acc \tsrc=MSA \tstatus=$status_msa \tcomment=$comment_msa";
        push @$msgs, $msg_msa;
        print STDERR "$subname: \$msg=$msg_msa\n";

        my $faa1 = '';
        # Order the resulting mat_peptides according to the start positions in the CDS
        my $refcds_ids = [ keys %$feats_msa ];
        $debug && print STDERR "$subname: \$refcds_ids='@$refcds_ids'\n";
        for my $i (0 .. $#{$refcds_ids}) {
            my $feats = $feats_msa->{$refcds_ids->[$i]};
            if (!$feats->[0]) {
                print STDERR "$subname: ERROR: NULL feature found for \$id=$refcds_ids->[$i]\n";
                next;
            }
            my $start1 = $feats->[0]->location->start;
            for my $j ($i+1 .. $#{$refcds_ids}) {
              $debug && print STDERR "$subname: \$feats_msa='\n". Dumper($feats_msa) . "End of \$feats_msa\n\n";
              my $feats = $feats_msa->{$refcds_ids->[$j]};
              my $start2 = $feats->[0]->location->start;
              $debug && print STDERR "$subname: \$start1=$start1 \$start2=$start2\n";
              if ($start1 > $start2) {
                my $temp = $refcds_ids->[$i];
                $refcds_ids->[$i] = $refcds_ids->[$j];
                $temp = $refcds_ids->[$j] = $temp;
              }
            }
        }
        $debug && print STDERR "$subname: \$refcds_ids='@$refcds_ids'\n";
        for my $id (@$refcds_ids) {
            my $feats = $feats_msa->{$id};
            if (!$feats->[0]) {
                print STDERR "$subname: ERROR: NULL feature found for \$id=$id\n";
                next;
            }
            $debug && print STDERR "$subname: \$feats=\n". Dumper($feats) . "End of \$feats\n\n";

            my $inseq = $feats->[0]->seq;
            # either print to STDERR or fasta file
            $faa1 .= Annotate_misc::generate_fasta( $feats);
            print STDERR "$subname: accession='$acc' CDS=".$feats->[0]->location->to_FTstring."\n";
#            $debug && print STDERR "$subname: \$faa1 = '\n$faa1'\n";

            Annotate_Verify::check_old_annotation( $acc, $faa1);
#            print STDERR "$subname: accession = '".$acc."'\n";
        }
        $debug && print STDERR "$subname: accession=$acc \$faa1 = '\n$faa1'\n";

        # Gets an array of Bio::PrimarySeq objects, directly from genbank file
        my ($feats_gbk, $comment_gbk) = Annotate_gbk::get_matpeptide( $gbk, $feats_msa,$exe_dir);
        my $status_gbk = $feats_gbk ? 'Success' : ($comment_gbk eq 'Not refseq, skip') ? 'Skip   ' : 'Fail   ';
        if (!$feats_gbk) {
            $count->{GBK}->{Fail}++;
            print STDERR "$subname: \$acc=$acc Empty result from Annotate_gbk::get_matpeptide comment=$comment_gbk\n";
        } else {
            $count->{GBK}->{Success}++;
            $debug && print STDERR "$subname: \$feats_gbk='\n$feats_gbk'\nEnd of \$feats_gbk\n\n";
        }
        $msg_gbk = "$acc \tsrc=GBK \tstatus=$status_gbk \tcomment=$comment_gbk";
        if ($feats_gbk && $faa1 && $acc !~ /^NC_/i) {
            $msg_gbk .= ". Not refseq, take MSA instead";
        }
        print STDERR "$subname: \$msg_gbk=$msg_gbk\n";
        push @$msgs, $msg_gbk;
        $debug && print STDERR "$subname: \$msgs=\n". Dumper($msgs) . "End of \$msgs\n\n";

        # $faa1 vs. $feats_gbk,
        # Take gbk if the genome is refseq,
        # Take $faa1 if exists
        # otherwise take gbk.
        $faa = Annotate_gbk::combine_msa_gbk( $acc, $faa1, $feats_gbk);
        print STDERR "$subname: \$acc=$acc \$faa=\n". Dumper($faa) . "End of \$faa\n\n";

        if ( 1 && $faa && $outfile) {
                open my $OUTFH, '>', $outfile
                    or croak "Can't open '$outfile': $OS_ERROR";
                print {$OUTFH} $faa
                    or croak "Can't write to '$outfile': $OS_ERROR";
                close $OUTFH
                    or croak "Can't close '$outfile': $OS_ERROR";
        }

#   exit;

   } # for (my $ct = 0; $ct<=$#{@$accs}; $ct++)
#   $debug && print STDERR "$subname: \$refseqs = \n".Dumper($refseqs)."End of \$refseqs\n\n";
#   $debug && print STDERR "$subname: \$inseqs = \n".Dumper($inseqs)."End of \$inseqs\n\n";

    print STDOUT "accession \tsource \tstatus \tcomment\n";
    for my $msg (@$msgs) {
        print STDOUT "$msg\n";
    }
    my $summary = '';
    $summary = "\nStatistics of this run:\n";
    $summary .= "Total input genomes: ". ($#{@$accs}+1) ."\n";
    for my $key ('MSA','GBK') {
        $summary .= "$key: \t";
        for my $key2 ('Success','Fail','Empty') {
            $summary .= "$key2 \t$count->{$key}->{$key2}\t";
        }
        $summary .= "\n";
    }
    print STDOUT "$summary\n";
    print STDERR "$summary\n";

    print STDERR "$subname: \$msgs=\n". Dumper($msgs) . "End of \$msgs\n\n";
exit;

} # sub process_list1


sub process_list3 {
    my ($accs, $aln_fn, $dbh_ref, $exe_dir, $exe_name, $dir_path, $test1, $refseq_required) = @_;

    my $debug = 0 && $debug_all;
    my $subname = 'process_list3';

#    $debug && print STDERR "$subname: \$accs = \n".Dumper($accs)."End of \$accs\n\n";
    my $refseqs = {};
    my $inseqs  = [];

   # get all CDS for seqs and refseqs, including the mat_peptides
   for (my $ct = 0; $ct<=$#{@$accs}; $ct++) {
        my $refpolyprots = [];
        my $polyprots    = [];
        my $number = $accs->[$ct]->[0];
        my $acc = $accs->[$ct]->[1];

        # get genbank file from MySQL database
        my $result;
        if ($dbh_ref && $dbh_ref->isa('GBKUpdate::Database')) {
            $result = $dbh_ref->get_genbank($acc);
            $result = $result->[0];
        } elsif (-r "$acc") { # This is really a filename
            open my $in, '<', "$acc" or croak "$0: Couldn't open $acc: $OS_ERROR";
            $result = do { local $/; <$in>};
            close $in or croak "$0: Couldn't close $acc: $OS_ERROR";
            print STDERR "$subname: \$result='".substr($result,0,79)."'\n";
        } else {
            print STDERR "$subname: Couldn't find genome \$acc='$acc'\n";
            next;
        }
#        print STDERR "\n$subname: result from database is '@$result'\n";

        if (!$result) {
            print STDERR "$subname: \$ct = $ct \$acc='$acc' result is empty\n";
            next;
        }

        # Now run the annotation
        my $gbk = $result;
#        print STDERR "$subname: \$result=$result->[0]\n";
        my $in_file2 = IO::String->new($gbk);
        my $in  = Bio::SeqIO->new( -fh => $in_file2, -format => 'genbank' );
        # Only take 1st sequence from each gbk (Note: gbk can hold multiple sequences, we ignore all after 1st)
        my $inseq = $in->next_seq();
        my $taxid = $inseq->species->ncbi_taxid;
        print STDERR "\n$subname: \$ct = $ct \$acc='$acc' \$taxid=$taxid processing\n";
#        print STDOUT "$subname: \$ct = $ct \$acc='$acc' \$taxid=$taxid processing\n";

        # determine the refseq, and get the CDS/mat_peptides in refseq
        $refpolyprots = Annotate_Util::get_refpolyprots( $refseqs, $inseq, $exe_dir);
#        $debug && print STDERR "$subname: \$refpolyprots = $#{@$refpolyprots}\n";
#        $debug && print STDERR "$subname: \$refpolyprots = \n".Dumper($refpolyprots)."End of \$refpolyprots\n\n";

        if ($#{@$refpolyprots} <0) {
            print STDERR "$subname: There is a problem getting refseq for ".$acc.". Skipping\n";
            print STDERR "$subname: \$#{\@\$refpolyprots} ".$#{@$refpolyprots}.".\n";
            next;
        }

        # According to refseq, get the CDS/mat_peptides in inseq, use bl2seq to determine if the CDS matches
        my $num_cds;
        ($polyprots, $num_cds) = Annotate_Util::get_polyprots( $inseq, $refpolyprots);
#        $debug && print STDERR "$subname:    \$polyprots = $num_cds\n";
        $debug && print STDERR "$subname: \$polyprots = \n".Dumper($polyprots)."End of \$polyprots\n\n";

        # add refseq to hash, add inseq to array
        if ($#{@$refpolyprots} >=0 && $refpolyprots->[0]->[1]->primary_tag eq 'CDS') {
            if (!exists($refseqs->{$refpolyprots->[0]->[1]->seq->accession_number})) {
                $refseqs->{$refpolyprots->[0]->[1]->seq->accession_number} = $refpolyprots;
            }
        } else {
            print STDERR "$subname: There is a problem with ".$acc.". Skipping\n";
            print STDERR "$subname: \$#{\@\$refpolyprots} ".$#{@$refpolyprots}.". Skipping\n";
            next;
        }

        my $n = [keys %$polyprots];
        $debug && print STDERR "$subname: polyprotein CDS for acc=".$acc." \$n=$#{@$n}\n";
        if ($#{@$n} <0) {
            print STDERR "$subname: Can't find any polyprotein CDS for acc=".$acc." \$n=$#{@$n}. Skipping\n";
#            print STDOUT "$subname: Can't find any polyprotein CDS for acc=".$acc." \$n=$#{@$n}. Skipping\n";
            next;
        }
        push @$inseqs, [$acc, $polyprots];

        if ($test1) {
            $debug && print STDERR "$subname: Processed one file, exit for debug\n";
            last;
        }
   }
   $debug && print STDERR "$subname: \$refseqs = \n".Dumper($refseqs)."End of \$refseqs\n\n";
   $debug && print STDERR "$subname: \$refseqs =".length($refseqs)."\n";
   $debug && print STDERR "$subname: \$inseqs = \n".Dumper($inseqs)."End of \$inseqs\n\n";
   $debug && print STDERR "$subname: \$inseqs =".$#{@$inseqs}."\n";

   $debug && print STDERR "$subname: Finished loading all seqs and refseqs, ready to run MUSCLE\n\n";

   $debug && print STDERR "$subname: Number of input genome is $#{@$inseqs}.\n";
   if ($#{@$inseqs}<0) {
       $debug && print STDERR "$subname: Nothing to process. Return.\n";
       return;
   }

   my $feats_all;
   $feats_all = Annotate_Muscle::muscle_profile( $refseqs, $inseqs, $aln_fn,$exe_dir);
#   $debug && print STDERR "$subname: \$feats_all=\n". Dumper($feats_all->[0]) . "End of \$feats_all\n\n";

   for (my $i=0; $i<=$#{@$feats_all}; $i++) {
       my $feats = $feats_all->[$i];
       for (my $j=0; $j<=$#{@$feats}; $j++) {
            my $feats_new = $feats->[$j];
            $debug && print STDERR "$subname: \$j=$j \$feats_new=\n".Dumper($feats_new)."End of \$feats_new\n\n";
            next if (!$feats_new);

            my $outfile = '';
            my $accession_number = $feats_new->[0]->seq->accession_number;
            $outfile = $accession_number . '_matpept_muscle.faa' if (!$debug);
            my $faa1 = Annotate_misc::generate_fasta( $feats_new, $outfile, '');
            $debug && print STDERR "$subname: \$outfile=$outfile\n";
            print STDERR "$subname: \$j=$j/$#{@$feats} accession = '".$accession_number."'\n";
            print STDERR "$subname: \$faa1 = '\n$faa1'\n";

            if ( 0 && $outfile) {
                open my $OUTFH, '>', $outfile
                    or croak "Can't open '$outfile': $OS_ERROR";
                print {$OUTFH} $faa1
                    or croak "Can't write to '$outfile': $OS_ERROR";
                close $OUTFH
                    or croak "Can't close '$outfile': $OS_ERROR";
            }

            Annotate_Verify::check_old_annotation( $accession_number, $faa1);
       }

   }

} # sub process_list3


=head2 list_dir_files

List the files with given pattern in a folder

=cut

sub list_dir_files {
    my ($list_fn, $ptn) = @_;

    my $debug = 0 && $debug_all;
    my $subname = 'Annotate_misc::list_dir_files';

    $debug && print STDERR "$subname: \$list_fn=$list_fn\n";
    $debug && print STDERR "$subname: \$ptn=$ptn\n";
    my @files = ();
    my $accs = [];
    if (!-d $list_fn) {
        croak("$subname: Couldn't locate accession file/directory: $list_fn: $OS_ERROR");
    }

    # if input -l file is directory
    opendir(DIR, $list_fn)
           or croak("$subname: Couldn't open dir $list_fn: $OS_ERROR");
    @files = sort readdir(DIR)
           or croak("$subname: Couldn't read dir $list_fn: $OS_ERROR");
    closedir(DIR)
           or croak("$subname: Couldn't close dir $list_fn: $OS_ERROR");
    $debug && print STDERR "$subname: \@files='@files'\n";

    for (my $f = 0; $f<=$#files; $f++) {
            my ($number, $acc);
            my $file = $files[$f];
            chomp $file;
            if ($file !~ /$ptn/) { # Keep the gbk files
#            if ($file !~ /^\s*([^\s]+)(\.(gb|gbk|genbank))\s*$/) { # Keep the gbk files
               $debug && print STDERR "$subname: skipping file: '$file'\n";
               next;
            } else {
                $number = $#{@$accs}+1;
                $acc = "$list_fn/$file";
#                print STDERR "\$1=\'$1\'\n";
            }
            push @$accs, [$number, $acc];
    }
    $debug && print STDERR "$subname: \@\$accs=".Dumper($accs)."\n";

    return $accs;
} # sub list_dir_files


=head2 Usage

Print the useage according to calling sub
 msa_annotate.pl:     for deployment
 msa_annotate_dev.pl: for testing

=cut

sub Usage {
    my ($source) = @_;

    my $debug = 0 && $debug_all;
    my $subname = 'Annotate_misc::Usage';

    my $usages = {};
    $usages = {
                'msa_annotate.pl' =>
	"$source Ver$VERSION
Usage:  -d directory to find the input genome file
        -i name of the input genbank file
        -l name of input folder",
                'msa_annotate_dev.pl' =>
"$source Ver$VERSION
Usage:  -d directory to find the input genome file
        -i name of the input genbank file
        -l name of input folder,
           or the name of a text file of accessions, and genbank file to be obtained from a database (MySQL)",
              };

    my $example = "
    E.g.: $source -d ./ -i NC_001477.gb
    E.g.: $source -d ./ -l all_genomes/ \n\n";

    my $usage = "$subname: default, need input for calling program, production or development";
    if (exists($usages->{$source})) {
        $usage = $usages->{$source};
        $usage .= $example;
    }

    return $usage;
} # sub Usage


1;
