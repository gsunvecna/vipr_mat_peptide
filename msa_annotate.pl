#!/usr/bin/perl

use strict;
use warnings;
use version; our $VERSION = qv('1.1.7'); # Feb 05, 2013
#use lib ("/home/peptide/loader/ext/BioPerl"); # This is for account peptide on Northrop machine 33
use Getopt::Long;
use English;
use Carp;
use Data::Dumper;

use Bio::SeqIO;
use Bio::Seq;
#use Bio::Tools::Run::Alignment::Muscle;
use IO::String;

## Path to the CLUSTALW binaries. You need to configure this, if clustalw is not in the path already
#	BEGIN {$ENV{MUSCLEDIR} = '/net/home/gsun/prog/clustalw/clustalw-2.0.12'}
#Check for CLUSTALW installation on start:
#$ENV{CLUSTALDIR} or croak 'CLUSTALDIR must be defined in your environment';
#(-e $ENV{CLUSTALDIR}) or croak "CLUSTALDIR ($ENV{CLUSTALDIR}) does not exist";

## Path to the MUSCLE binaries. You need to configure this, if muscle is not already in the path
#	BEGIN {$ENV{MUSCLEDIR} = '/net/home/gsun/prog/muscle/mus37'}
use Annotate_Download;
use Annotate_misc;

my $debug = 0;

##//README//##
#
# msa_annotate.pl
#
# This script takes a viral genome in genbank format, uses a refseq to annotate any polyproptein within.
# It outputs a file named as <accession>_matpept_msa.faa with the annotated mat_peptides in fasta format
# if the result comes from alignment. Otherwise, outputs to a file <accession>_matpept_gbk.faa when the
# result comes from genbank.
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
# For single genome
# ./msa_annotate.pl -d [dir_path] -i [inputFile.gb]
# For multiple input genomes within a directory
# ./msa_annotate.pl -d [dir_path] -l [directory]
# e.g.
# ./msa_annotate.pl -d ./ -i NC_001477_test.gb >> out.txt 2>> err.txt
# ./msa_annotate.pl -d ./ -l test >> test/out.txt 2>> test/err.txt
#
#    Authors: Chris Larsen, clarsen@vecna.com; Guangyu Sun, gsun@vecna.com;
#    September 2011
#
#################


## //EXECUTE// ##

# Program locations, may leave blank if the programs are accessible from the prompt
my $progs = {
    # Location of MUSCLE executible, such as /net/home/gsun/prog/muscle/mus37/muscle
#    muscle  => '/home/peptide/loader/ext/muscle3.7/muscle',
    # Location of CLUSTALW executible, such as /net/home/gsun/prog/clustalw/clustalw-2.0.12/clustalw2
#    clustalw  => '/home/peptide/loader/ext/autoCuration/clustalw';
};

# Get user-defined options
my $infile;
my $list_fn = '';
my $dir_path = './';
my $checkRefseq = 0;
my $checkTaxon = 0;

my $exe_dir  = './';
my $exe_name = $0;
($exe_dir, $exe_name) = ($1, $2) if ($exe_name =~ /^(.*[\/])([^\/]+[.]pl)$/i);
print STDERR "$exe_name: $0 $VERSION executing from command='$0 @ARGV' ".POSIX::strftime("%m/%d/%Y %H:%M:%S", localtime)."\n";
my $useropts = GetOptions(
                 "checkrefseq"  => \ $checkRefseq,    # Check any update for RefSeqs from genbank
                 "checktaxon"  => \ $checkTaxon,    # Check any update for taxon from genbank
                 "d=s"  => \ $dir_path,    # Path to directory
                 "i=s"  => \ $infile,      # [inputFile.gbk]
                 "l=s"  => \ $list_fn,     # list of the gbk file
                 );
$dir_path =~ s/[\/]$//;
$list_fn =~ s/[\/]$//;

$debug && print STDERR "$exe_name: \$exe_dir    = $exe_dir\n";
$debug && print STDERR "$exe_name: input path   = $dir_path\n";
$debug && print STDERR "$exe_name: input $checkRefseq = $checkRefseq\n";
$debug && print STDERR "\n";

# Either a genbank file or a folder/list of genbank files is required
if ($checkTaxon) {
#    my $count = Annotate_Def::checkAllTaxon( $exe_dir);
    print STDERR "\n$exe_name: start to check taxon.\n";
    my $count = Annotate_Download::checkAllTaxon( $exe_dir);
    printf STDERR "\n$exe_name: finished checking %d taxon, exit.\n\n", $count;
    exit(1);

} elsif ($checkRefseq) {
#    $count = Annotate_Def::checkAllRefseq( $exe_dir);
    print STDERR "\n$exe_name: start to check RefSeq.\n";
    my $count = Annotate_Download::checkAllRefseq( $exe_dir);
    printf STDERR "\n$exe_name: finished checking %d RefSeqs, exit.\n\n", $count;
    exit(1);

} elsif (!$infile && !$list_fn) {
    print Annotate_misc::Usage( $exe_name, $exe_dir);
    exit(1);
}

# Now, get all accessions (or files) to be processed
my $dbh_ref = undef;
my $aln_fn  = undef;
my $accs = [];
if ($infile) {

    print STDERR "$exe_name: Input genbank file '$dir_path/$infile'\n";
    if ($infile !~ m/[.](gb|gbk|genbank)/i) {
        print STDERR "$exe_name: WARNING: please make sure input genbank file '$infile' is in genbank format\n";
    }

    # Now run the annotation
    my $accs = [];
    push @$accs, [$#{$accs}+1, "$dir_path/$infile"];

    $dbh_ref = undef;
    Annotate_misc::process_list1( $accs, $aln_fn, $dbh_ref, $exe_dir, $exe_name, $dir_path, $progs);

} elsif ("$dir_path/$list_fn") {

    print STDERR "$exe_name: Input accession are in dir/file '$dir_path/$list_fn'\n";
    my $accs = [];
    if (-d "$dir_path/$list_fn") {
        $dir_path = "$dir_path/$list_fn";
        # if input -l file is directory
        my $ptn = '^\s*([^\s]+)(\.(gb|gbk|genbank))\s*$';
        $accs = Annotate_misc::list_dir_files( "$dir_path", $ptn);
        my $n = $#{$accs}+1;
        print STDERR "$exe_name: from directory: $dir_path, found $n gbk files.\n";
        for my $j (0 .. $#{$accs}) {
            print STDERR "$exe_name: \$acc[$j]=$accs->[$j]->[1]\n";
        }
        print STDERR "\n";
#        $debug && print STDERR "$exe_name: \$accs=\n".Dumper($accs)."End of \$accs\n\n";

        $dbh_ref = undef;
    }

    if ( 1 ) {
        # MSA for each genome
        Annotate_misc::process_list1( $accs, $aln_fn, $dbh_ref, $exe_dir, $exe_name, $dir_path, $progs);
    }

}

    print "\n$exe_name: finished.\n\n";
    print STDERR "\n$exe_name: finished.\n\n";

exit(0);

1;
