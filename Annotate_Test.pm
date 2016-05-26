package Annotate_Test;

use strict;
use warnings;
use English;
use Carp;
use Data::Dumper;

#use version;
our $VERSION = qw('1.3.1'); # May 10 2015
use File::Temp qw/ tempfile tempdir /;
use Bio::SeqIO;
use Bio::Seq;
use IO::String;

use Annotate_Align;      # Ran MSA, and get ready for actually cut the mat_peptides
use Annotate_Def;        # Have the approved RefSeqs, also load taxon info and definition of gene symbols
use Annotate_Download;   # Download the RefSeqs and taxon info from NCBI, and check against data stored in file
use Annotate_gbk;	 # for annotation from genbank
use Annotate_Math;       # Get the mat_peptide location, based on alignment with reference
use Annotate_Util;       # assemble the feature for newly annotated mat_peptide, plus other things like checking
use Annotate_Verify;     # Check the quality of the set of annotation

my $debug_all = 0;

####//README//####
#
# Annotate_Test contains the subroutines for testing the installation
#
#    Authors Guangyu Sun, gsun@vecna.com;
#    February 2015
#

############## Subroutines ##############

# Run vipr_mat_peptide.pl script for the genome

sub run_vipr_mat_peptide {
    my ($acc, $dirIn, $exe_dir, $removeFile) = @_;

    my $debug = 0;
    my $subName = 'run_vipr_mat_peptide';

    my $fmatpept = '';
    my $err = '';

    if (!-e "$dirIn/$acc" && ($acc ne 'usage')) {
        $err = "ERROR: Input file '$dirIn/$acc' doesn't exist";
        return ($fmatpept, $err);
    }

    my $dir_path = $dirIn;
    my $result = '';

    # Run vipr_mat_peptide.pl
    my $ferr = "$acc.err";
    my $fout = "$acc.out";
    $debug && print STDERR "$subName: \$fout=$fout \$ferr=$ferr\n";
    $debug && print STDERR "$subName: \$fout=$dir_path/$fout exists, will be overwritten\n" if (-e "$dir_path/$fout");
    $debug && print STDERR "$subName: \$ferr=$dir_path/$ferr exists, will be overwritten\n" if (-e "$dir_path/$ferr");
    my $cmd = '';
    $cmd = "$exe_dir/vipr_mat_peptide.pl"; # Program
#    $cmd = "$exe_dir/vipr_mat_peptide_dev.pl"; # Program
    if ($acc eq 'usage') {
#        $cmd .= " -d $dir_path -i $acc "; # input parameters
        $cmd .= " >$dir_path/$fout 2>$dir_path/$ferr"; # Redicting stdout and stderr to files
    } else {
        $cmd .= " -d $dir_path -i $acc "; # input parameters
        $cmd .= " >$dir_path/$fout 2>$dir_path/$ferr"; # Redicting stdout and stderr to files
    }
    $result = `$cmd`; # run the vipr_mat_peptide.pl script
    $debug && print STDERR "$subName: \$cmd=$cmd\n";
    $debug && print STDERR "$subName: \$result=\n".Dumper($result)."End of \$result\n";

    if (!-e "$dir_path/$ferr") {
        $err = "ERROR: After running vipr_mat_peptide.pl, can't find STDERR file '$dirIn/$ferr'";
        return ($fmatpept, $err);
    }

    # find the name of new outfile from the $ferr
    open my $in, '<', "$dir_path/$ferr" or croak "$0: Couldn't open $dir_path/$ferr: $OS_ERROR";
    while ( <$in> ) {
        if ($_ =~ /outfile=.*[\/]([^\/]+)[.]$/) {
            $fmatpept = $1;
            $debug && print STDERR "$subName: \$fmatpept='$fmatpept'\n";
        } elsif ($_ =~ /.+MSA.+Fail.+comment=(.+)$/) {
            $err = $1;
            last;
        }
    }
    close $in or croak "$0: Couldn't close $dir_path/$ferr: $OS_ERROR";
    $debug && print STDERR "$subName: \$fmatpept='$fmatpept' \$err='$err'\n";
    `rm $dir_path/$fout $dir_path/$ferr` if ($removeFile);

    return ($fmatpept, $err);
} # sub run_vipr_mat_peptide


# Compares the new output with eixsting output

sub cmp_vipr_mat_peptide {
    my ($acc, $fmatpept, $dirIn, $dirOut) = @_;

    my $debug = 0;
    my $subName = 'cmp_vipr_mat_peptide';

    $debug && print STDERR "$subName: \$acc=$acc \$fmatpept='$fmatpept' \$dirIn='$dirIn'\n";
    my $result = '';
    my $err = '';

  if ($acc eq 'usage') {
    $fmatpept = "usage.out";
  } else {
    $acc = $1 if ($acc =~ m/^([^.]+)[.](gb|gbk|genbank)/);
    if (!$acc) {
        $err = "ERROR: Genome '$acc' is invalid";
        return ($result, $err);
    }
    my $list = '';
    $list = `ls $dirOut/${acc}_matpept_???.faa `;
    $debug && print STDERR "$subName: \$list='$list'\n";
    if (!$list) {
            $err = "No standard output for $acc in $dirOut/";
            return ($result, $err);
    }
    if (!$fmatpept) {
            $err = "ERROR: output file '$fmatpept' is invalid";
            return ($result, $err);
    }
    if (!-e "$dirIn/$fmatpept") {
            $err = "ERROR: new output file '$dirIn/$fmatpept' doesn't exist";
            return ($result, $err);
    }
    if (!-e "$dirOut/$fmatpept") {
            $err = "ERROR: standard output file '$dirOut/$fmatpept' doesn't exist";
            return ($result, $err);
    }
  }

    # Compare new <acc>_matpept_msa.faa with standard result
    $result = "";
    my $cmd = "diff $dirOut/$fmatpept $dirIn/$fmatpept";
    if ( 1 ) {
        $result = `$cmd`;
    } else {
        # Place holder for detailed comparison of mat_peptide sequences
        # in case if the defline in fasta file has changed, or others
        
    }
    if ($result) {
        $result = "$cmd\n" . $result;
        $err = "ERROR: Difference with old result";
    }
    return ($result, $err);
} # sub cmp_vipr_mat_peptide


1;

