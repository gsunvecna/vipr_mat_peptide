package Annotate_Verify;

use strict;
use warnings;
use English;
use Carp;
use Data::Dumper;

use version; our $VERSION = qv('1.1.1');
use Bio::SeqIO;
use Bio::Seq;
use Bio::AlignIO;
use Bio::Tools::Run::Alignment::Muscle;
use IO::String;

my $debug_all = 1;

####//README//####
#
# Annotate_Util contains the core functions to perform annotation based on both MUSCLE and bl2seq alignment
#
#    Authors Chris Larsen, clarsen@vecna.com; Guangyu Sun, gsun@vecna.com
#    September 2010
#
##################

## //EXECUTE// ##

# Get refseq object if $refseq_fn is given
#  my $dbh_ref;
#  my $refseq = &get_refseq($refseq_fn, $dir_path);


=head2 check_ranges

Checks the new annotations against the CDS to see if all features (in an array) are consecutive without any gap, as a gap (shown as >=?=<) indicates potential problem.

=cut

sub check_ranges {
    my ($feats) = @_;

    my $debug = 0 && $debug_all;
    my $has_gap = 0;
    my $cds = $feats->[0];
    my @str = ('');
#    print "check_ranges: \$feats is a ".$feats."\n";
    $debug && print "check_ranges: \$cds = \n".Dumper($cds)."\n";
    $feats = [sort {$a->location->start <=> $b->location->start} @{$feats}];
#    print "check_ranges: \$feats is a ".$feats."\n";
    my $head_tail = '';
    my $head_tail_print_length = 5;
    for (my $i=0; $i<=$#{@$feats}; $i++) {
        my $feat = $feats->[$i];
        next if ($feat->primary_tag ne 'mat_peptide' && $feat->primary_tag ne 'sig_peptide');
#        print "check_ranges: \$feat is a ".ref($feat)."\n";
        my $loc = $feat->location;
#        print "check_ranges: \$loc is a ".ref($loc)."\n";
#        print "check_ranges: \$loc = \n".Dumper($loc)."\n";
#        print "check_ranges: \$cds_>loc = \n".Dumper($cds->location)."\n";
        my $level = 0;
        my $loc2 = $feats->[$i+$level+1];
        my @add;
        while ($loc2 && ($loc2->start < $loc->end) && ($loc2->end <= $loc->end)) {
          if (!defined($str[$level+1])) {
            $str[$level+1] = $str[0];
            $str[$level+1] =~ s/\S/ /g;
          }
#          $add[$level+1] .= '('.$loc2->start.'..'.$loc2->end.')';
          $add[$level+1] .= sprintf(" (%4s..%4s)", $loc2->start, $loc2->end);
          $level++;
          $loc2 = $feats->[$i+$level+1];
        }

#        my $s = ' ('.$loc->start.'..'.$loc->end.')';
        my $s = sprintf(" (%4s..%4s)", $loc->start, $loc->end);
        my $gap = '';

        # if the first feature doesn't start at 1st residue of CDS
        if ($i == 0 && ($cds->location->isa('Bio::Location::Simple') || $cds->location->isa('Bio::Location::Split')) && $loc->start != $cds->location->start) {
          $has_gap = 1;
          $gap = $loc->start - $cds->location->start;
          $gap = '='.$gap.'=<';
          $s = $gap . $s;
          $gap =~ s/\S/ /g;
          for my $add (1 .. $#add) {
            $add[$add] = $gap . $add[$add];
          }
        }

        if (!$loc2) {   # This indicates $loc is the last feat
#          print "check_ranges: \$loc->end=".$loc->end." \$cds->end=".$cds->location->end."\n";
          if ($loc->end == $cds->location->end-3 || $loc->end == $cds->location->end) {   # If there is overhang
            $gap = '';
          } else {
            $has_gap = 1;
            $gap = $cds->location->end - $loc->end;
            $gap = '>='. $gap .'=';
          }
        } elsif ($loc2->start != $loc->end+1) {   # If there is a gap between subsequent features
          $has_gap = 1;
          $gap = $loc2->start - $loc->end -1;
          $gap = '>='. $gap .'=<';
        } else {
          $gap = ' ';
        }
        $str[0] .= $s.$gap;
        $gap =~ s/\S/ /g;
        for my $j (1 .. $#add) {
            $str[$j] .= $add[$j].$gap;
        }
        $s .= $gap;
        $s =~ s/\S/ /g;
        for my $j ($level+1 .. $#str) {
            $str[$j] .= $s;   # pad other strings
        }

        for my $j (0 .. $#str) {
            $debug && print "check_ranges: \$str[$j]='$str[$j]'\n";
        }

        # check the head/tail of each mat_peptide
        if ( 1 ) {
            $s = [ $feat->get_tag_values('translation') ];
            $s = $s->[0];
            $head_tail .= substr($s,0,$head_tail_print_length).'.';
            $head_tail .= substr($s,length($s)-$head_tail_print_length);
            $head_tail .= '.^.' if ($loc2);
#            print "check_ranges: head_tail=$s";
        }

        $i += $level;
    }
    print "check_ranges: ". $cds->seq->accession_number ."\thead_tail= $head_tail\n";

    if ( 1 ) {
        for my $j (0 .. $#str) {
            print "check_ranges: ". $cds->seq->accession_number ."\t\$str[$j]  =$str[$j]\n";
        }
    }

    return ($has_gap, \@str);
} # sub check_ranges


=head2 print_matpept_ranges

Takes a title, and 2 arrays of features, print the locations on top of each other, leave blanks for extra locations, then indicate the differences with ^

=cut

sub print_matpept_ranges {
    my ($title, $metbod, $ranges, $ranges_new, $cds_start) = @_;

    my $different = 0 && $debug_all;
    foreach my $ct (0 ..$#{@{$ranges_new}}) {
        if (!defined($ranges->[$ct]) || !defined($ranges_new->[$ct]) || $ranges->[$ct]->location->start != $ranges_new->[$ct]->location->start || $ranges->[$ct]->location->end != $ranges_new->[$ct]->location->end) {
              $different = 1;
              last;
        }
    }

    print "\n";
    if (!@{$ranges}) {
            print "!!! New annotation for mat_peptide from $metbod and genbank\n";
    } elsif ($different || $#{@{$ranges}}!=$#{@{$ranges}}) {
            print "!!! Different mat_peptide from $metbod and genbank\n";
    } else {
            print "!!! Identical between $metbod and genbank\n";
    }
    print "$title\n";
    my $str1='';
    my $str2='';
    my $str3='';
    my $len1 = $#{@$ranges};
    my $len2 = $#{@$ranges_new};
    my $i=0;
    my $j=0;
    while ($i <= $len1 || $j <= $len2) {
#        print "i=$i j=$j\n";
        $str1 .= ' ' if ($str1);
        $str2 .= ' ' if ($str2);
        $str3 .= ' ' if ($str3);
        my $range;
        $range = $ranges->[$i]->location if ($ranges->[$i]);
        my $range_new;
        $range_new = $ranges_new->[$j]->location if ($ranges_new->[$j]);
#        print "print_ranges: \$range is a ".ref($range)."\n";
#        print "print_ranges: \$range_new is a ".ref($range_new)."\n";
        my $fm1 = '%d';
        my $fm2 = '%d';
        my $s;
        if ((!$range) || ($range_new && $range->start > $range_new->end)) {
              # insert gap to $range
              $fm1 = get_format(0, $range_new->start);
              $fm2 = get_format(0, $range_new->end);
#              $s = "(".$range_new->start."..".$range_new->end.")";
              $s = sprintf("($fm1..$fm2)",$range_new->start,$range_new->end);
              $str2 .= $s;
              if ($cds_start) {
                $s = sprintf("($fm1..$fm2)",($range_new->start-$cds_start)/3+1,($range_new->end-$cds_start-2)/3+1);
                $str3 .= $s;
              }
              $s =~ s/\S/ /g;
              $str1 .= $s;
 #             print "print_ranges: i=$i j=$j \$str1='$str1'\n";
 #             print "print_ranges: i=$i j=$j \$str2='$str2'\n\n";
              $j++;
              next;
        }
        if ((!$range_new) || ($range && $range_new->start > $range->end)) {
              # insert gap to $range_new
#              my $s;
              $fm1 = get_format($range->start, 0);
              $fm2 = get_format($range->end, 0);
#              $s = "(".$range->start."..".$range->end.")";
              $s = sprintf("($fm1..$fm2)",$range->start,$range->end);
              $str1 .= $s;
              $s =~ s/\S/ /g;
              $str2 .= $s;
              $str3 .= $s;
 #             print "print_ranges: i=$i j=$j \$str1='$str1'\n";
 #             print "print_ranges: i=$i j=$j \$str2='$str2'\n\n";
              $i++;
              next;
        }

#        $range = $ranges->[$i]->location;
#        $range_new = $ranges_new->[$j]->location;
        $fm1 = get_format($range->start, $range_new->start) if ($range && $range_new);
        $fm2 = get_format($range->end, $range_new->end) if ($range && $range_new);
        $str1 .= sprintf("($fm1..$fm2)",$range->start,$range->end) if ($range);
        $str2 .= sprintf("($fm1..$fm2)",$range_new->start,$range_new->end) if ($range_new);
        $str3 .= sprintf("($fm1..$fm2)",($range_new->start-$cds_start)/3+1,($range_new->end-$cds_start-2)/3+1) if ($range_new && $cds_start);
#        print "print_ranges: \$ranges_new is a ".Dumper($ranges_new->[$j])."\n";
  #      print "print_ranges: i2=$i j=$j \$str1='$str1'\n";
  #      print "print_ranges: i2=$i j=$j \$str2='$str2'\n\n";
        $i++; $j++;
    }

#    print "\n";
    print "\$ranges    ='$str1'\n";
    print "\$ranges_new='$str2'\n";
    print " diff      ='". diff_2str($str1,$str2) ."'\n" if (@{$ranges});
    print "\$ranges_new='$str3'\n" if ($cds_start);

    return;
} # sub print_matpept_ranges


# to print out the numbers, padding the smaller number with space according to the larger number
sub get_format {
    my ($n1, $n2) = @_;

    my $debug = 0 && $debug_all;
    my $fmt = undef;

    return $fmt if (!defined($n1) && !defined($n2));

#    print "get_format: ".length($n1)."\n";
    $fmt = ($n1>$n2)? length($n1) : length($n2);
    $fmt = '%'.$fmt.'d';
#    print "get_format: fmt=$fmt\n";

    return $fmt;
} # sub get_format


=head2 check_alignment

Takes a refseq sequence object, and a feature object (CDS), match against each other.
Runs a pairwise bl2seq query looking for a segment of the second sequence
which matches the first string.

Returns nothing if a match of 50% identity was not found, or a hash with four return values:
  match => the matching sequence string
  start => match start point in the ortholog sequence
  ident => %identity between the match and original shortstring (0-100)

=cut

sub check_alignment {
    # not done yet 4/07/2010
    my ($hsp_cds, $feats, $ref_mats) = @_;

    my (@out);
    print STDOUT "check_alignment: \$hsp_cds=\n".Dumper($hsp_cds)."end of \$hsp_cds\n\n";
#    print STDOUT "check_alignment: \$feats=\n".Dumper($feats)."end of \$feats\n\n";

    foreach my $i (0 .. $#{@$feats}) {
        my $feat = $feats->[$i];
        print STDOUT "check_alignment: \$i=$i \$feat=\n".Dumper($feat)."end of \$feat\n\n";
        my $start = $feat->location->start - $hsp_cds->start('HIT');
        my $len   = $feat->location->end - $feat->location->start;

        $out[0] .= substr($hsp_cds->query_string, $start, $len) .'^|^';
        $out[1] .= substr($hsp_cds->homology_string, $start, $len) .'^|^';
        $out[2] .= substr($hsp_cds->hit_string, $start, $len) .'^|^';

        foreach my $s (@out) {
          print "out=$s\n";
        }
    }
    return;
} # sub check_alignment


sub diff_2str {
    my ($str1,$str2) = @_;

    my $debug = 1 && $debug_all;
    my ($diff);
    my (@a1, @a2);
    @a1 = split(//, $str1);
    @a2 = split(//, $str2);
#      print "$#a1=@a1\n";
#      print "$#a2=@a2\n";
    my $len = ($#a1 >= $#a2) ? $#a1 : $#a2;
    for my $i (0 .. $len) {
        if (defined($a1[$i]) && defined($a2[$i]) && $a1[$i] eq $a2[$i]) {
            if ($a1[$i] eq "\n") {
                $diff .= "\n";
            } else {
                $diff .= '.';
            }
        } else {
            my $s1 = '';
            $s1 .= $a1[$i-1] ? $a1[$i-1] : ' ';
            $s1 .= $a1[$i]   ? $a1[$i]   : ' ';
            $s1 .= $a1[$i+1] ? $a1[$i+1] : ' ';
            my $s2 = '';
            $s2 .= $a2[$i-1] ? $a2[$i-1] : ' ';
            $s2 .= $a2[$i]   ? $a2[$i]   : ' ';
            $s2 .= $a2[$i+1] ? $a2[$i+1] : ' ';
#            $debug && print STDERR "diff_2str: $s1 $s2\n";
#            $diff .= ($#a1 > $#a2) ? $a1[$i] : $a2[$i];
            $diff .= ($i<=$#a1) ? $a1[$i] : $a2[$i];
        }
#      print STDERR "\$i=$i \$diff=$diff\n";
    }
    return $diff;
} # sub diff_2str


1;
