#!/usr/bin/perl

use strict;
use warnings;
use version; our $VERSION = qv('1.1.5'); # Sep 20, 2012
use File::Temp qw/ tempfile tempdir /;
use Getopt::Long;
use English;
use Carp;
use Data::Dumper;

use Bio::SeqIO;
use Bio::Seq;
use Bio::AlignIO;
use Bio::Tools::Run::StandAloneBlast;
use Bio::Tools::Run::Alignment::Muscle;
use IO::String;

#use lib qw(/net/home/gsun/northrop/matpeptide/msa_annotate-1.1.3/);
use GBKUpdate::Configuration;
use GBKUpdate::Database;

## Path to the CLUSTALW binaries. You need to configure this, if clustalw is not in the path already
#	BEGIN {$ENV{MUSCLEDIR} = '/net/home/gsun/prog/clustalw/clustalw-2.0.12'}
#Check for CLUSTALW installation on start:
$ENV{CLUSTALDIR} or croak 'CLUSTALDIR must be defined in your environment';
(-e $ENV{CLUSTALDIR}) or croak "CLUSTALDIR ($ENV{CLUSTALDIR}) does not exist";

## Path to the MUSCLE binaries. You need to configure this, if muscle is not in the path already
#	BEGIN {$ENV{MUSCLEDIR} = '/net/home/gsun/prog/muscle/mus37'}
#Check for MUSCLE installation on start:
#$ENV{MUSCLEDIR} or croak 'MUSCLEDIR must be defined in your environment';
#(-e $ENV{MUSCLEDIR}) or croak "MUSCLEDIR ($ENV{MUSCLEDIR}) does not exist";

use Annotate_gbk;		# for annotation from genbank
use Annotate_Muscle;		# for live MUSCLE run
#use Annotate_Muscleprofile;	# for live MUSCLE profile run
#use Annotate_bl2seq;		# for live bl2seq run
use Annotate_Util;
use Annotate_Verify;
use Annotate_misc;

my $debug_all  = 1; # Used to turn off all debugging code
my $debug      = 1;
my $withDB     = 1; # If the genomes are optionally stored in a MySQL database
my $test1      = 0;

##//README//##
#
# msa_annotate_dev.pl
#
# This script uses a refseq to annotate the polyprotein in a genome file
# It outputs the annotated mat_peptides in fasta named as <accession>_matpeptide.faa with.
#
# INPUT: dir of input file, input genome file name, ticket (used as subfolder), optional refseq file name
#
# OUTPUT: fasta file containing the annotated mat_peptide sequence
#
# DEPENDENCIES:
# This script calls perl and uses BioPerl modules.
# Specify the blast executable location in your environment as directly below!
#
# USAGE:
# For single input genome
# ./msa_annotate.pl -r [refseq] -d [dir_path] -i [inputFile.gb]
# For multiple input genomes within a directory
# ./msa_annotate.pl -d [dir_path] -l [directory]
# For multiple input genomes whose accession is in a list, and gbk files are in MySQL database
# ./msa_annotate.pl -d [dir_path] -l [list_file.txt]
# e.g.
# ./msa_annotate.pl -d ./ -i NC_001477_test.gb >> out.txt 2>> err.txt
# ./msa_annotate.pl -d ./ -l test >> test/out.txt 2>> test/err.txt
# ./msa_annotate.pl -d test -l nuccore_result.txt >> test/out.txt 2>> test/err.txt
#
#	Authors Chris Larsen, clarsen@vecna.com; Guangyu Sun, gsun@vecna.com
#	September 2011
#
#################

## Path to the BLAST binaries. You must configure this!
my $blast_path = '/home/gsun/prog/blast/blast-2.2.20/bin';
my $td = tempdir( CLEANUP => 1 );  # temp dir (threadsafe)

# Program locations, may leave blank if the programs are accessible from the prompt
my $progs = {
    # Location of MUSCLE executible, such as /net/home/gsun/prog/muscle/mus37/muscle
#    muscle  => '/home/dbadmin/loader/ext/muscle3.7/muscle',
    # Location of CLUSTALW executible, such as /net/home/gsun/prog/clustalw/clustalw-2.0.12/clustalw2
#    clustalw  => '/home/dbadmin/loader/ext/autoCuration/clustalw';
};

# Get user-defined options
my $refseq_required = 0;
my $refseq_fn = '';
my $infile    = '';
my $list_fn   = '';
my $dir_path  = './';
my @aln_fn    = ();
my $aln_fn    = \@aln_fn;

my $exe_dir  = './';
my $exe_name = $0;
if ($exe_name =~ /^(.*[\/])([^\/]+[.]pl)$/i) {
    $exe_dir  = $1;
    $exe_name = $2;
}
print STDERR "$exe_name: $0 executing...\n";
print STDERR "$exe_name: command='$0 @ARGV'\n";
my $useropts = GetOptions(
                 "d=s"  => \ $dir_path,    # Path to directory
                 "i=s"  => \ $infile,      # [inputFile.gbk]
                 "l=s"  => \ $list_fn,     # directory with the gbk file, or list of accessions from genbank search
                 "r=s"  => \ $refseq_fn,   # refseq in gbk file
                 "a=s"  => \ @aln_fn,      # alignment file
                 "t=i"  => \ $test1,       # quit after one run
                 );
$dir_path =~ s/[\/]$//;
$list_fn =~ s/[\/]$//;

print STDERR "$exe_name: Directory=$dir_path \talignment file='@$aln_fn'\n";

# Either a genbank file or a folder/list of genbank files is required
if (!$infile && !$list_fn) {
    print Annotate_misc::Usage( $exe_name);
    exit(1);
}

if ($refseq_fn) {
    print STDERR "$exe_name: Refseq supplied in $refseq_fn.\n";
    $refseq_required = 1;
}

##################

## //EXECUTE// ##

# Get refseq object if $refseq_fn is given
  my $dbh_ref = undef;
  my $refseq;
  if ($refseq_fn) {
      print STDERR "$exe_name: refseq read from $refseq_fn\n";
      $refseq = Annotate_Util::get_refseq( $refseq_fn);
  }

# Now go through the input file of sequences
if ("$infile") {

    print STDERR "$exe_name: Input genbank file '$dir_path/$infile'\n";

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
#        $debug && print STDERR "$exe_name: \$accs = \n".Dumper($accs)."End of \$accs\n\n";

        $dbh_ref = undef;

    } elsif (-f "$dir_path/$list_fn") {
        croak("$exe_name: Need to turn on \$withDB to connect to MySQL database.") if (!$withDB);
        my $list_file;
        open $list_file, '<', "$dir_path/$list_fn"
           or croak("$exe_name: Found $list_fn, but couldn't open");
        while (<$list_file>) {
            my ($number, $acc);
            chomp;
#           print $_;
            if ($_ =~ /^\s*(\d+)\s*:\s*([^\s]+)\s*$/) { # '440: NC_002031'
                $number = $1;
                $acc = $2;
                print STDERR "\$1=\'$1\' \t\$2=\'$2\'\n";
            } elsif ($_ =~ /^\s*([^\s.]+)[.]\d+\s*$/) { # FJ888392.1
                $number = $#{$accs} +1;
                $acc = $1;
#                print STDERR "\$1=\'$1\'\n";
            } elsif ($_ =~ /^(\d+)[.] /) { # 1. Norovirus Hu/
                $number = $1;
#                $_ = <$list_file>;
#                $_ = <$list_file>;
#                print STDERR "\$_='$_'\n";
                while (<$list_file>) {
                  if ($_ =~ /^([a-z0-9_]+)[.]\d /i) { # FR695417.1 GI:308232890
                    $acc = $1;
                    print STDERR "$exe_name: \$number=\'$number\' \$acc=$acc\n";
                  }
                  last if ($_ =~ /^$/x)
                }
            } else {
               0 && $debug && print STDERR "$exe_name: skipping line: '$_'\n" if ($_ || $_ ne "/n");
               next;
            }
            push @$accs, [$number, $acc];

        }
        close $list_file or croak "$0: Couldn't close $dir_path/$list_fn: $OS_ERROR";

        print STDERR "\n$exe_name: finished reading list file: $dir_path/$list_fn, found $#{$accs} gbk files.\n\n";
#        $debug && print STDERR "$exe_name: \$accs = \n".Dumper($accs)."End of \$accs\n\n";

        # configuration file used to connect to MySQL
        my $cfg_file = $exe_dir .'vbrc_retrieveGBK_mysql.cfg';
        my $config_ref = GBKUpdate::Configuration->new($cfg_file);
#        $debug && print "$0: \$config_ref=\n". Dumper($config_ref) . "\n";

        $dbh_ref = GBKUpdate::Database->new($config_ref);
#       $debug && print "$0: \$dbh_ref ='". Dumper($dbh_ref) ."'\n";

    } else {
        croak("$exe_name: Couldn't locate accession file/directory: $list_fn: $OS_ERROR");
    }

    if ( 1 ) {
        # MSA for each genome
        $debug && print STDERR "$exe_name: sub Annotate_misc::process_list1 called\n";
        Annotate_misc::process_list1( $accs, $aln_fn, $dbh_ref, $exe_dir, $exe_name, $dir_path, $progs);
    } else {
        # Giant MSA. For large set, requires long time, run out of time/memory, and possibly give wrong result b/c of gaps
        if ( !$debug ) {
            croak("$exe_name: sub Annotate_misc::process_list3 is for debug only. Quit.\n");
        }
        Annotate_misc::process_list3( $accs, $aln_fn, $dbh_ref, $exe_dir, $exe_name, $dir_path, $progs);
    }
}

    print "\n$exe_name: finished.\n\n";
    print STDERR "\n$exe_name: finished.\n\n";

exit(0);


1;
