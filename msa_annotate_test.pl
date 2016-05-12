#!/usr/bin/perl

use strict;
use warnings;
use version; our $VERSION = qv('1.2.0'); # Apr 12, 2013
use Getopt::Long;
use Data::Dumper;
use English;
use Carp;

use Annotate_misc;
use Annotate_Def;

my $debug = 0;

##################//README//#####################################################
#                                                                               #
# msa_annotate_test.pl                                                          #
#                                                                               #
# This script looks for the test/ directory, and find all genbank files in it.  #
# For each genbank file, this script runs the msa_annotate.pl script, and       #
# compares the result with result stored in test/output/.                       #
# Output for each file Success/Fail status, and the summary in the end.         #
#                                                                               #
# INPUT: none. It looks for any genbenk file inside the test/ directory         #
# OUTPUT: for each genbank file, Success/Fail status, and summary when done     #
#                                                                               #
# DEPENDENCIES:                                                                 #
# This script needs msa_annotate.pl package correctly installed                 #
#                                                                               #
# USAGE:                                                                        #
# ./msa_annotate_test.pl                                                        #
# ./msa_annotate_test.pl -f Togaviridae                                         #
# ./msa_annotate_test.pl -f Togaviridae -debug                                  #
# ./msa_annotate_test.pl -i DQ241304.gb (the file needs to be in ./test)        #
#                                                                               #
# Author: Guangyu Sun, Vecna Technologies, Inc, gsun@vecna.com;                 #
# October 2012                                                                  #
#                                                                               #
#################################################################################

    my $families = {
            'Arenaviridae'    => [],
            'Bunyaviridae'    => [],
            'Caliciviridae'   => [],
            'Coronaviridae'   => [],
            'Filoviridae'     => [],
            'Flaviviridae'    => [],
            'Hepeviridae'     => [],
            'Herpesviridae'   => [],
            'Paramyxoviridae' => [],
            'Picornaviridae'  => [],
            'Poxviridae'      => [],
            'Reoviridae'      => [],
            'Rhabdoviridae'   => [],
            'Togaviridae'     => [],
                   };

## //EXECUTE// ##

my $infile = '';
my $FAMILY = '';
my $useropts = GetOptions(
                 "i=s"  => \ $infile,      # [inputFile.gb]
                 "f=s"  => \ $FAMILY,      # viral family to test
                 "debug"=> \ $debug,       # debug
                 );
my $exe_dir  = './';
my $exe_name = $0;
if ($exe_name =~ /^(.*[\/])([^\/]+[.]pl)$/i) {
    $exe_dir  = $1;
    $exe_name = $2;
}
$exe_dir =~ s/[\\\/]$//;
print STDERR "$exe_name: $0 V$VERSION executing ".POSIX::strftime("%m/%d/%Y %H:%M:%S", localtime)."\n";
print STDOUT "$exe_name: $0 V$VERSION executing ".POSIX::strftime("%m/%d/%Y %H:%M:%S", localtime)."\n";
print STDERR "$exe_name: command='$0 @ARGV'\n";
if ($FAMILY && !exists($families->{$FAMILY})) {
    print STDERR "$exe_name: \$FAMILY='$FAMILY' is entered with invalid value, all families will be tested\n";
    $FAMILY = '';
}
print STDERR "$exe_name: requested \$FAMILY='$FAMILY'\n";

$debug && print STDERR "$exe_name: \$exe_dir    = '$exe_dir'\n";
$debug && print STDERR "\n";

# Make sure test/ and test/output directories exist
my $dirIn = "$exe_dir/test";
my $dirOut = "$exe_dir/test/output";
if (!-d "$dirIn") {
    croak("$exe_name: Can't fine input directory '$dirIn', abort.");
} elsif (!-d "$dirOut") {
    croak("$exe_name: Can't find output directory '$dirOut', abort.");
} else {
    print STDERR ("$exe_name: Found both input directory '$dirIn' and output directory '$dirOut'.\n");
}

# Now, get all accessions (or files) to be processed
my $accs = [];
my $dir_path = "$dirIn";
my $removeAnnotationResult = 0;
if ($infile) {
    if (-e "$dir_path/$infile") {
        $accs = [[0, "$dir_path/$infile"]];
        $removeAnnotationResult = 0;
        $debug && print STDERR ("$exe_name: Found required input gbk file at '$dir_path/$infile'\n");
    } else {
        print STDERR ("$exe_name: Can't find required input gbk file at '$dir_path/$infile', abort\n");
        exit(1);
    }
} elsif (-d "$dirIn") {
    my $ptn = '^\s*([^\s]+)(\.(gb|gbk|genbank))\s*$'; # pattern for genbank file name
    $accs = Annotate_misc::list_dir_files( "$dir_path", $ptn); # the returned accs already contain $dir_path/

    unshift @$accs, [$#{$accs}+1, 'usage']; # add usage to test
}
#$debug && print STDERR "$exe_name: \$accs=\n". Dumper($accs) . "End of \$accs\n";

    my $n = $#{$accs}+1;
    print STDERR "$exe_name: from Input directory: $dir_path, found $n gbk files.\n";
    print STDOUT "$exe_name: from Input directory: $dir_path, found $n gbk files.\n";

    my $summ = {total=>0, pass=>0, fail=>0, skip=>0};
    my ($min, $max) = (0, $#{$accs});
    if ( 0 ) {
        $min = ($#{$accs}-5<0) ? 0 : $#{$accs}-5;
        $max = ($#{$accs}-4<0) ? 0 : $#{$accs}-4; # for debugging
    }
    printf("$exe_name: Starting to test %d genomes: %d-%d\n", $max-$min+1, $min, $max);

    # Separate the accessions into different families
#    my $families = {};
    for my $j ($min .. $max) {
        my $acc = $accs->[$j]->[1];
        $acc = $1 if ($acc =~ m/^.*[\/]([^\/]+)$/); # take only the <acc>.gb, delete directory
        $debug && print STDERR "\n";
        $debug && print STDERR "$exe_name: \$acc[$j]=$accs->[$j]->[1] \$acc=$acc\n";

      if ($acc eq 'usage') {
        push @{$families->{usage}}, [$#{$families->{usage}}+1, $acc];
      } else {
        my $in  = Bio::SeqIO->new( -file => "$dirIn/$acc", -format => 'genbank' );
        $debug && print STDERR "$exe_name: \$exe_dir=$exe_dir \$in=$in\n";
        my $inseq = $in->next_seq();
        $debug && print STDERR "$exe_name: \$exe_dir=$exe_dir \$inseq=$inseq\n";
        my $taxid = $inseq->species->ncbi_taxid;
#        $debug && print STDERR "$exe_name: \$taxid=$taxid \$exe_dir=$exe_dir\n";
        my $ti = Annotate_Def::getTaxonInfo( $taxid, $exe_dir);
        $debug && print STDERR "$exe_name: \$acc=$acc \$taxid=$taxid \$ti='@$ti'\n";
        my $fam = (defined($ti->[6])) ? $ti->[6] : 'unknown';

        $debug && print STDERR "$exe_name: \$acc=$acc \$fam=$fam \$FAMILY=$FAMILY\n";
        push @{$families->{$fam}}, [$#{$families->{$fam}}+1, $acc];
      }
    }
for my $fam (sort keys %$families) {
    $debug && print STDERR "$exe_name: \$fam=$fam\n";
    my $accs = $families->{$fam};
    for my $j (0 .. $#{$accs}) {
        $debug && print STDERR "$exe_name: \$fam=$fam \$accs='@{$accs->[$j]}'\n";
    }
}
#    $debug && print STDERR "$exe_name: \$families=\n".Dumper($families)."End of \$families\n";

    # Process each genbank file
    my $count = 0;
for my $fam ( sort keys %$families) {
    my $accs = $families->{$fam};
    my $msg1 = sprintf (" Testing %-15s %d", "$fam,", scalar @$accs);

    print STDERR "$msg1\n";
    print STDOUT "$msg1\n";
    for my $j (0 .. $#{$accs}) {
        $count++;
        my $acc = $accs->[$j]->[1];
        my $result = '';
        my $err = '';
        my $msg1 = sprintf (" Testing #%-4s %-15s %-13s", "$count,", "$fam,", $acc);
        if (!$FAMILY || $fam eq $FAMILY) {
            ($result, $err) = &run_msa_annotate( $acc, $dirIn, $removeAnnotationResult);
        } else {
            $err = "Family not requested. Skip";
        }

        if (!$err) {
            ($result, $err) = &cmp_msa_annotate( $acc, $result, $dirIn, $dirOut);
        }
        $debug && print STDERR "$exe_name: \$acc='$acc' New:OLd \$result='$result' \$err='$err'\n";

        my $msg = '';
        if ($err =~ m/ERROR/) {
            $debug && print STDERR "$exe_name: Problem with \$acc=$acc\n";
            $msg .= sprintf(" Test #%-4s %-15s %-13s fail, comment=$err\n", "$count,", "$fam,", "$acc,");
            $msg .= sprintf(" Test #%-4s %-15s %-13s \$result=$result\n", "$count,", "$fam,", "$acc,");
            $summ->{fail}++;
            $summ->{total}++;
        } elsif ($err =~ m/Skip$/i) {
            $msg .= sprintf(" Test #%-4s %-15s %-13s skip", "$count,", "$fam,", "$acc,");
            $msg .= sprintf(", comment=$err") if $err;
            $msg .= sprintf("\n");
            $summ->{skip}++;
        } else {
#            printf STDERR (" Test #%-4s %-15s %-13s pass", "$count,", "$fam,", "$acc,");
#            printf STDERR (", comment=$err") if $err;
#            printf STDERR ("\n");
            $msg .= sprintf(" Test #%-4s %-15s %-13s pass", "$count,", "$fam,", "$acc,");
            $msg .= sprintf(", comment=$err") if $err;
            $msg .= sprintf("\n");
            $summ->{pass}++;
            $summ->{total}++;
        }
        print STDOUT "$msg";
        print STDERR "$msg";
    }

}
    print STDERR "\n";
#    $debug && print STDERR "$exe_name: \$accs=\n".Dumper($accs)."End of \$accs\n\n";

    print "Tests performed: total=$summ->{total}, fail=$summ->{fail}, pass=$summ->{pass}; skip=$summ->{skip}\n";



print "\n$exe_name: finished.\n\n";
print STDERR "\n$exe_name: finished.\n\n";

exit(0);


#### Subroutines ####

# Run msa_annotate.pl script for the genome

sub run_msa_annotate {
    my ($acc, $dirIn, $removeFile) = @_;

    my $debug = 0;
    my $subName = 'run_msa_annotate';

    my $fmatpept = '';
    my $err = '';

    if (!-e "$dirIn/$acc" && ($acc ne 'usage')) {
        $err = "ERROR: Input file '$dirIn/$acc' doesn't exist";
        return ($fmatpept, $err);
    }

    my $dir_path = $dirIn;
    my $result = '';

    # Run msa_annotate.pl
    my $ferr = "$acc.err";
    my $fout = "$acc.out";
    $debug && print STDERR "$subName: \$fout=$fout \$ferr=$ferr\n";
    $debug && print STDERR "$subName: \$fout=$dir_path/$fout exists, will be overwritten\n" if (-e "$dir_path/$fout");
    $debug && print STDERR "$subName: \$ferr=$dir_path/$ferr exists, will be overwritten\n" if (-e "$dir_path/$ferr");
    my $cmd = '';
    $cmd = "$exe_dir/msa_annotate.pl"; # Program
#    $cmd = "$exe_dir/msa_annotate_dev.pl"; # Program
    if ($acc eq 'usage') {
#        $cmd .= " -d $dir_path -i $acc "; # input parameters
        $cmd .= " >$dir_path/$fout 2>$dir_path/$ferr"; # Redicting stdout and stderr to files
    } else {
        $cmd .= " -d $dir_path -i $acc "; # input parameters
        $cmd .= " >$dir_path/$fout 2>$dir_path/$ferr"; # Redicting stdout and stderr to files
    }
    $result = `$cmd`; # run the msa_annotate.pl script
    $debug && print STDERR "$subName: \$cmd=$cmd\n";
    $debug && print STDERR "$subName: \$result=\n".Dumper($result)."End of \$result\n";

    if (!-e "$dir_path/$ferr") {
        $err = "ERROR: After running msa_annotate.pl, can't find STDERR file '$dirIn/$ferr'";
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
} # sub run_msa_annotate


# Compares the new output with eixsting output

sub cmp_msa_annotate {
    my ($acc, $fmatpept, $dirIn, $dirOut) = @_;

    my $debug = 0;
    my $subName = 'cmp_msa_annotate';

    my $result = '';
    my $err = '';

  if ($acc eq 'usage') {
    $fmatpept = "usage.out";
  } else {
    $acc = $1 if ($acc =~ m/^([^.]+)[.](gb|gbk|genbank)/);
    $debug && print STDERR "$subName: \$acc=$acc\n";
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
} # sub cmp_msa_annotate


1;

