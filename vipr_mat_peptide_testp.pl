#!/usr/bin/perl

use strict;
use warnings;
use version; our $VERSION = qv('1.2.1'); # May 10, 2013
use Getopt::Long;
use Data::Dumper;
use English;
use Carp;
use Parallel::ForkManager; # To install this module: 1)perl -MCPAN -e shell; 2)force install Parallel::ForkManager

use Annotate_misc;
use Annotate_Def;
use Annotate_Test;

my $debug = 0;

##################//README//#####################################################
#                                                                               #
# vipr_mat_peptide_test.pl                                                          #
#                                                                               #
# This script looks for the test/ directory, and find all genbank files in it.  #
# For each genbank file, this script runs the vipr_mat_peptide.pl script, and       #
# compares the result with result stored in test/output/.                       #
# Output for each file Success/Fail status, and the summary in the end.         #
#                                                                               #
# INPUT: none. It looks for any genbenk file inside the test/ directory         #
# OUTPUT: for each genbank file, Success/Fail status, and summary when done     #
#                                                                               #
# DEPENDENCIES:                                                                 #
# This script needs vipr_mat_peptide.pl package correctly installed                 #
#                                                                               #
# USAGE:                                                                        #
# ./vipr_mat_peptide_test.pl                                                        #
# ./vipr_mat_peptide_test.pl -f Togaviridae                                         #
# ./vipr_mat_peptide_test.pl -f Togaviridae -debug                                  #
# ./vipr_mat_peptide_test.pl -i DQ241304.gb (the file needs to be in ./test)        #
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

print STDERR "command='$0 @ARGV'\n";
print STDERR "$0 V$VERSION executing ".POSIX::strftime("%m/%d/%Y %H:%M:%S", localtime)."\n";
print STDOUT "$0 V$VERSION executing ".POSIX::strftime("%m/%d/%Y %H:%M:%S", localtime)."\n";
my $nproc         = 1;
my $infile = '';
my $FAMILY = '';
my $useropts = GetOptions(
                 "i=s"  => \ $infile,      # [inputFile.gb]
                 "f=s"  => \ $FAMILY,      # viral family to test
                 "debug"=> \ $debug,       # debug
                 "nproc=s"  => \$nproc,         # number of processes, for parallel processing
                 );
my $exe_dir  = './';
my $exe_name = $0;
if ($exe_name =~ /^(.*[\/])([^\/]+[.]pl)$/i) {
    $exe_dir  = $1;
    $exe_name = $2;
}
$exe_dir =~ s/[\\\/]$//;
if ($FAMILY && !exists($families->{$FAMILY})) {
    print STDERR "$exe_name: \$FAMILY='$FAMILY' is entered with invalid value, all families will be tested\n";
    $FAMILY = '';
}
print STDERR "$exe_name: requested \$FAMILY='$FAMILY'\n";
$nproc = ($nproc<0) ? 0 : ($nproc>20) ? 20 : $nproc;
print STDERR "$exe_name: \$nproc=$nproc\n";

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
    } elsif ($infile =~ /usage/i) {
        unshift @$accs, [$#{$accs}+1, 'usage']; # add usage to test
        $debug && print STDERR ("$exe_name: added 'usage' to list of accessions\n");
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

    # Set up parallel process, each input gbk/fasta file is one fork
    # Set $max_procs = 0; for debugging, ie. not forking
    # Set $max_procs to (intended number of process)-1 for processing
    my $max_procs = ($nproc>0) ? $nproc : 0;
    my $pm =  Parallel::ForkManager->new( $max_procs);
    $debug && print STDERR "$exe_name: ForkManager started with \$max_procs=$max_procs\n";
    $pm->run_on_start(
        sub { my ($pid, $ident) = @_;
#            print STDERR "** $ident started, pid: $pid\n";
        }
    );
    $pm->run_on_finish(
        sub { my ($pid, $exit_code, $ident) = @_;
          if (!$exit_code) {
#            print STDERR "** $ident finished PID $pid and exit_code: $exit_code\n";
          } else {
#            print STDERR "** $ident ERROR: finished PID $pid and exit_code: $exit_code\n";
#            print STDERR "** $ident ERROR: needs checking\n";
          }
        }
    );

    print STDERR "$msg1\n";
    print STDOUT "$msg1\n";
    my ($timeStart, $timeEnd);
    for my $j (0 .. $#{$accs}) {
        $count++;
        # starts the fork in child, and returns to for loop in main
        my $pid = $pm->start("$j:$accs->[$j]->[1]") and next;

        $timeStart = time();
        my $acc = $accs->[$j]->[1];
        my $result = '';
        my $err = '';
        my $msg1 = sprintf (" Testing #%-4s %-15s %-13s", "$count,", "$fam,", $acc);
        if (!$FAMILY || $fam eq $FAMILY) {
            ($result, $err) = Annotate_Test::run_vipr_mat_peptide( $acc, $dirIn, $exe_dir, $removeAnnotationResult);
        } else {
            $err = "Family not requested. Skip";
        }

        if (!$err) {
            ($result, $err) = Annotate_Test::cmp_vipr_mat_peptide( $acc, $result, $dirIn, $dirOut);
        }
        $debug && print STDERR "$exe_name: \$acc='$acc' New:OLd \$result='$result' \$err='$err'\n";

        my $accOut = $acc;
        $accOut = $1 . "_matpept_msa.faa" if ($acc =~ /^(.+)[.]gb/i);
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
        } elsif ($err && -e "$dirOut/$accOut") {
            $err .= ", but found $dirOut/$accOut";
            $debug && print STDERR "$exe_name: Problem with \$acc=$acc\n";
            $msg .= sprintf(" Test #%-4s %-15s %-13s fail, comment=$err\n", "$count,", "$fam,", "$acc,");
            $msg .= sprintf(" Test #%-4s %-15s %-13s \$result=$result\n", "$count,", "$fam,", "$acc,") if ($result);
            $summ->{fail}++;
            $summ->{total}++;
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
        my $t = time() - $timeStart;
        printf STDOUT ("%3d sec\t$msg", $t);
        print STDERR "$msg";
        $pm->finish($acc); # do the exit in the child process, pass an exit code to finish
    } # acc
    $pm->wait_all_children; # Wait for all child processes to finish

} # family
    print STDERR "\n";
#    $debug && print STDERR "$exe_name: \$accs=\n".Dumper($accs)."End of \$accs\n\n";

    print STDERR "Tests performed: total=$summ->{total}, fail=$summ->{fail}, pass=$summ->{pass}; skip=$summ->{skip}\n";
    print "Tests performed: total=$summ->{total}, fail=$summ->{fail}, pass=$summ->{pass}; skip=$summ->{skip}\n";



print "\n$exe_name: finished.\n\n";
print STDERR "\n$exe_name: finished.\n\n";

exit(0);


