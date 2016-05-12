#!/usr/bin/perl

use strict;
use warnings;
use version; our $VERSION = qv('1.1.1'); # Oct 13, 2010
use Getopt::Long;
use Data::Dumper;
use English;
use Carp;

use Bio::SeqIO;
use Bio::Seq;
use Bio::Tools::Run::Alignment::Muscle;
use IO::String;

## Path to the MUSCLE binaries. You need to configure this, if muscle is not in the path already
#	BEGIN {$ENV{MUSCLEDIR} = '/net/home/gsun/prog/muscle/mus37'}
use Annotate_Util;
use Annotate_Verify;
use Annotate_Muscle;

my $debug = 0;

##//README//##
#
# msa_annotate.pl
#
# This script takes a viral genome in genbank format, uses a refseq to annotate any polyproptein within.
# It outputs a file named as <accession>_matpeptide.faa with the annotated mat_peptides in fasta format.
#
# INPUT: directory of input file, genome file name
#
# OUTPUT: fasta file containing the annotated mat_peptide sequence
#
# DEPENDENCIES:
# This script calls perl and uses BioPerl modules.
# Specify the MUSCLE executable location in your environment as directed above!
#
# USAGE:
# ./msa_annotate.pl -d [dir_path] -i [inputFile.fasta]
# e.g.
# ./msa_annotate.pl -d ./ -i NC_001477_test.gb >> out.txt 2>> err.txt
#
#    Authors: Chris Larsen, clarsen@vecna.com; Guangyu Sun, gsun@vecna.com;
#    September 2010
#
#################

## //OPTIONS// ##
# (-d [dir_path] -i [inputFile.gb])

##################

## //EXECUTE// ##

# Get user-defined options
my $infile;
my $dir_path = './';
my $useropts = GetOptions(
                 "d=s" => \ $dir_path,     # Path to data directory
                 "i=s" => \ $infile,       # [inputFile.gbk]
                 );

my $exe_dir  = './';
my $exe_name = $0;
if ($exe_name =~ /^(.*[\/])([^\/]+[.]pl)$/i) {
    $exe_dir  = $1;
    $exe_name = $2;
}

print STDOUT "$exe_name: $0 V$VERSION executing with infile=$infile.\n";
print STDERR "\n$exe_name: $0 V$VERSION executing with infile=$infile.\n";
$debug && print STDERR "$exe_name: \$exe_dir    = $exe_dir\n";
$debug && print STDERR "$exe_name: input path   = $dir_path\n";
$debug && print STDERR "\n";

(!$infile) && croak("$exe_name: Need input genome in genbank format to proceed. Use -i [infile] option.");

# Now go through the input file, all sequences
if ($infile) {

    $debug && print STDERR "$exe_name: Input genbank file '$infile'\n";
    my $gbk;
    if (!-e "$dir_path/$infile") {
        croak "$exe_name: Couldn't locate file '$infile'\n";
    } else {
        open my $in, '<', "$dir_path/$infile" or croak "$0: Couldn't open $infile: $OS_ERROR";
        $gbk = do { local $/; <$in>};
        close $in or croak "$0: Couldn't close $infile: $OS_ERROR";
        0 && print STDERR "$exe_name: \$gbk='$gbk'\n";
    }

    my $feats_new = Annotate_Muscle::annotate_1gbk( $gbk, $exe_dir);

    $debug && print STDERR "$exe_name: \$feats_new=\n". Dumper($feats_new) . "End of \$feats_new\n";
    if (!$feats_new) {
        print STDERR "$exe_name: ERROR: Annotate_Muscle::annotate_1gbk returned empty result.\n";
        print STDERR "$exe_name: ERROR: Exit.\n";
        exit(1);
    }

    print STDERR "\n";
    my $faa1;
    my $outfile = undef;
    for (my $i=0; $i<=$#{@$feats_new}; $i++) {
        my $feats = $feats_new->[$i];

        # following sub call returns the fasta file in a string, for debug
        $faa1 .= Annotate_Util::generate_fasta( $feats);
        $outfile = $dir_path .'/'. $feats->[0]->seq->accession_number . '_matpept.faa' if (!$debug);
        print STDERR "$exe_name: accession=".$feats->[0]->seq->accession_number." CDS=".$feats->[0]->location->to_FTstring."\n";
    }

    # Use following to write fasta file to disk
    open my $OUTFH, '>', $outfile or croak "Can't open '$outfile': $OS_ERROR";
    print {$OUTFH} $faa1 or croak "Can't write to '$outfile': $OS_ERROR";
    close $OUTFH or croak "Can't close '$outfile': $OS_ERROR";

    print STDERR "$exe_name: \$faa1 = '\n".$faa1."'End of \$faa1\n";
}

print "$exe_name: finished.\n";
print STDERR "$exe_name: finished.\n";

exit(0);

1;
