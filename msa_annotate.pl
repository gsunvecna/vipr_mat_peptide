#!/usr/bin/perl

use strict;
use warnings;
use version; our $VERSION = qv('1.1.5'); # Sep 20, 2012
use Getopt::Long;
use Data::Dumper;
use English;
use Carp;

use Bio::SeqIO;
use Bio::Seq;
use Bio::Tools::Run::Alignment::Muscle;
use IO::String;

## Path to the CLUSTALW binaries. You need to configure this, if clustalw is not in the path already
#	BEGIN {$ENV{MUSCLEDIR} = '/net/home/gsun/prog/clustalw/clustalw-2.0.12'}
#Check for CLUSTALW installation on start:
$ENV{CLUSTALDIR} or croak 'CLUSTALDIR must be defined in your environment';
(-e $ENV{CLUSTALDIR}) or croak "CLUSTALDIR ($ENV{CLUSTALDIR}) does not exist";

## Path to the MUSCLE binaries. You need to configure this, if muscle is not already in the path
#	BEGIN {$ENV{MUSCLEDIR} = '/net/home/gsun/prog/muscle/mus37'}
#use lib qw(/net/home/gsun/northrop/matpeptide/msa_annotate-1.1.2/);
#use Annotate_gbk;		# for annotation from genbank
#use Annotate_Muscle;
#use Annotate_Util;
#use Annotate_Verify;
use Annotate_misc;

my $debug = 0;

##//README//##
#
# msa_annotate.pl
#
# This script takes a viral genome in genbank format, uses a refseq to annotate any polyproptein within.
# It outputs a file named as <accession>_matpept_msagbk.faa with the annotated mat_peptides in fasta format.
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
#    muscle  => '/home/dbadmin/loader/ext/muscle3.7/muscle',
    # Location of CLUSTALW executible, such as /net/home/gsun/prog/clustalw/clustalw-2.0.12/clustalw2
#    clustalw  => '/home/dbadmin/loader/ext/autoCuration/clustalw';
};

# Get user-defined options
my $infile;
my $list_fn = '';
my $dir_path = './';

my $useropts = GetOptions(
                 "d=s"  => \ $dir_path,    # Path to directory
                 "i=s"  => \ $infile,      # [inputFile.gbk]
                 "l=s"  => \ $list_fn,     # list of the gbk file
                 );
$dir_path =~ s/[\/]$//;
$list_fn =~ s/[\/]$//;
my $exe_dir  = './';
my $exe_name = $0;
if ($exe_name =~ /^(.*[\/])([^\/]+[.]pl)$/i) {
    $exe_dir  = $1;
    $exe_name = $2;
}
print STDERR "$exe_name: $0 executing...\n";
print STDERR "$exe_name: command='$0 @ARGV'\n";

$debug && print STDERR "$exe_name: \$exe_dir    = $exe_dir\n";
$debug && print STDERR "$exe_name: input path   = $dir_path\n";
$debug && print STDERR "\n";

# Either a genbank file or a folder/list of genbank file is required
if (!$infile && !$list_fn) {
    print Annotate_misc::Usage( $exe_name);
    exit(1);
}

# Now, get all accessions (or files) to be processed
my $dbh_ref = undef;
my $aln_fn  = undef;
my $accs = [];
if ($infile) {

    print STDERR "$exe_name: Input genbank file '$infile'\n";

    # Now run the annotation
    my $accs = [];
    push @$accs, [$#{$accs}+1, "$dir_path/$infile"];

    $dbh_ref = undef;
    Annotate_misc::process_list1( $accs, $aln_fn, $dbh_ref, $exe_dir, $exe_name, $dir_path, $progs);

} elsif ("$dir_path/$list_fn") {

    print STDERR "$exe_name: Input accession are in file/dir '$dir_path/$list_fn'\n";
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
#        $debug && print STDERR "$exe_name: \$accs = \n".Dumper($accs)."End of \$accs\n\n";

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
