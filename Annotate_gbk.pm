package Annotate_gbk;

use strict;
use warnings;

use Bio::SeqIO;
use Bio::Seq;
use Getopt::Long;
use Data::Dumper;

use Annotate_Util;
use Annotate_Verify;

my $debug_all = 0;

##//README//##
#
# gbk_matpeptide.pl
#
# gbk_matpeptide.pl looks in a given directory (default is ./), and for each DNA
# genbank (gb|gbk|gbank|genbank) file, it looks for mat_peptide features, which
# are generated from polyprotein by viral genome. In case there is translation tag
# for the mat_peptide, it's taken as is. If no such translation tag, the DNA
# sequence is then extracted from the genome sequence according to the given
# range, and translated to AA sequencee. The script then writes out all AA
# mat_peptides to a fasta file, named after the input file plus ".faa".
# The id of each mat_peptide is given as following: first, either GI or NP number of
# parent CDS or protein_id, followed by range of the mat_peptide, followed by its
# 'db_xref', 'protein_id', 'product' in this order, whichever appears first.
# The script assumes that mat_peptide only exists in genome genbank files.
#
# INPUT: directory that contains genbank files
# INPUT: name of the genbank file. If absent, it looks for all genbank files in dir
# OUTPUT: a fasta file for each genbank file, containing AA sequences of mat_peptides
#
# USAGE:
# ./gbk_matpeptide.pl [-d <directory>] [-f <infile>]
# The dir input is optional, default value is ./. The -f option is also optional, the script
# looks for all genbank files if the -f option is absent.
#
# Author: Guangyu Sun, gsun@vecna.com
# Oct. 28 2009
#
#################

sub setDebugAll {
    my ($debug) = @_;
    $debug_all = $debug;
} # sub setDebugAll


=head2
 sub combine_msa_gbk combines the fasta string from MSA and GBK into
  1. one set if they are identical, or
  2. Keep both if different, or
  3. Keep which ever is present
   returns the mat_peptide sequences to a string.
=cut

sub combine_msa_gbk {
    my ($acc, $faa1, $faa3) = @_;

    my $debug = 0 && $debug_all;
    my $subname = 'combine_msa_gbk';

    my $faa = '';

    # Add to MSA result any extra mat_peptide in gbk that doesn't show up in MSA result 
    if ($faa1 && $faa3) {
    }

        # $faa1 vs. $feats_gbk,
        # Take $faa1 if exists and not refseq
        # Tame gbk if the genome is refseq, or no $faa1
    if (!$faa1) {
        $debug && print STDERR "$subname: \$acc=$acc has no MSA annotation, take \$faa3.\n";
        $faa = $faa3;
    } elsif ($acc =~ /^NC_\d+$/i) {
        $debug && print STDERR "$subname: \$acc=$acc is a refseq, take \$faa3.\n";
        $faa = ($faa3) ? $faa3 : $faa1;
    } else {
        $debug && print STDERR "$subname: \$acc=$acc is not a refseq, and has MSA annotation, take \$faa1.\n";
    	$faa = $faa1;
    }

   return ($faa);
} # sub combine_msa_gbk


## sub extract_mature_peptides searches a Bio::Seq::RichSeq object for features that are
## named mat_peptide, adds the translation tag to the feature if not present,
## and returns an array of Bio::PrimarySeq containing only id, alphabet, and seq.
sub extract_mature_peptides {
    my ($seq_obj, $feats_msa,$exe_dir) = @_;

    my $debug = 0 || $debug_all;
    my $subname = 'extract_mature_peptides';

    my $acc = $seq_obj->accession_number;
    my $comment = 'No CDS found in gbk';
    my @records = ();
    my ($parent_cds, $count);
    my $count_err = { total=>0, good=>0, mismatched_seq=>0, stop_in_seq=>0, null_CDS=>0};
# Bio::Seq::RichSeq => Bio::SeqFeature::Generic
    my @feat_objs = $seq_obj->get_SeqFeatures;
    for my $i (0 .. $#feat_objs) {
        my $feat_obj = $feat_objs[$i];
        $debug && print STDERR "$subname: $acc feat_obj #$i is '". $feat_obj->primary_tag ."'\n";
        if ($feat_obj->primary_tag eq "CDS") {
            $comment = 'No mat_peptide found in gbk';
            $parent_cds = $feat_obj;
            $count = 0;
            next;
        }
        next if ($feat_obj->primary_tag ne "mat_peptide");
        # Accession AM408911 has a mat_peptide before CDS. Simply ignore it for now
        if (!$parent_cds) {
            $count_err->{total}++;
            $count_err->{null_CDS}++;
            print STDERR "$subname: Found feat_obj #$i is '". $feat_obj->primary_tag ."' without parent CDS, ignored\n";
            next;
        }
        # See if this is a new mat_peptide, some genomes have duplicates, see NC_003899
        my $seen = 0;
        for my $j (0 .. $i-1) {
            my $f = $feat_objs[$j];
            if ($feat_obj->location->start==$f->location->start
                && $feat_obj->location->end==$f->location->end
                && $feat_obj->primary_tag eq $f->primary_tag) {
                $seen = 1;
                last;
            }
        }
        next if ($seen);
        # following process prints the literal DNA sequence, w/o consideration of backups
#	print '$feat_obj->seq='.$feat_obj->seq->seq ."\n";
	# for the join(12332..12358,12358..15117) in DQ848678, this produceds garbage
#	print '$feat_obj->seq='.$feat_obj->seq->translate->seq ."\n";

# The id of mat_peptide is given as: first, either GI (db_xref) or NP (protein_id) of
# parent CDS or protein_id, followed by range of the mat_peptide, followed by its
# 'db_xref', 'protein_id', and 'product' in this order, whichever appears first.
        my (@id);
        my $id = 0;
        # get the id from parent CDS feature, 1) GI, 2) NP
        $debug && print STDERR "$subname: \$parent_cds=\n". Dumper($parent_cds) . "End of \$parent_cds\n\n";
        if ($parent_cds->has_tag('db_xref')) {
                @id = $parent_cds->get_tag_values('db_xref');
                for my $id1 (@id) {
                   if ($id1 =~ /^GI:/i) {
                	$id = 'CDS='. $id1 .'|';
                	last;
                   }
                }
        }
        if (($id !~ /^CDS=GI:/) && $parent_cds->has_tag('protein_id')) {
                @id = $parent_cds->get_tag_values('protein_id');
                $id = 'CDS='. $id[0] .'|';
        }

        # In case $id is still empty
        if (($id !~ /^CDS=/) ) {
             my @tags = ('db_xref', 'protein_id', 'product'); # per Client request
             for my $tag (@tags) {
                if ($parent_cds->has_tag($tag)) {
                   @id = $feat_obj->get_tag_values($tag);
                   $id .= 'CDS='. $id[0] .'|';
                   last;
                }
             }
        }
        $id = "CDS=unknown|" if (($id !~ /^CDS=/) ); # indicates a potential problem
	    $id = "Ver=$acc.".$parent_cds->entire_seq->{'_version'}.'|' .$id;
	    $id = "ACC=$acc|" .$id;
        $id .= '|'; # Save the space for refseq=

        # add the range of mat_peptide
        my $s = '';
        $s = Annotate_Util::get_DNA_loc( $feat_obj);
#        $id .= 'Loc='. $feat_obj->location->to_FTstring ."|";
        $id .= 'Loc='. $s .'|';

        # Take care of codon_start if CDS has it, this is only needed for the first mat_peptide
        if ($feat_obj->location->start == $parent_cds->location->start) {
            my $codon_start = 1;
            if ($parent_cds->has_tag('codon_start')) {
                $codon_start = [$parent_cds->get_tag_values('codon_start')];
                $codon_start = $codon_start->[0];
            }
            if ($codon_start != 1 && !$feat_obj->has_tag('codon_start')) {
                $feat_obj->add_tag_value('codon_start', $codon_start);
                $id .= 'cstart='. $codon_start;
            }
        }
        $id .= '|';
        $s = Annotate_Util::get_AA_loc( $feat_obj, $parent_cds);
        $id .= 'AA='. $s .'|';

        # add gene_symbol
        my $symbol = Annotate_Def::get_gene_symbol( $feat_obj,$exe_dir);
    my $productName = '';
  if (1) {
    my @tags = ('product');
    for my $tag (@tags) {
        if ($feat_obj->has_tag($tag)) {
           @id = $feat_obj->get_tag_values($tag);
           $productName = $id[0];
           last;
        }
    }
    if ((!$symbol || $symbol eq 'unk') && length($productName)<=6 && length($productName)>0) {
#        $symbol = $productName;
    }
  }

        $id .= 'symbol='. $symbol .'|';

        # add the id of mat_peptide
        my @tags = ('db_xref', 'protein_id'); # per Client request
        for my $tag (@tags) {
             if ($feat_obj->has_tag($tag)) {
                   @id = $feat_obj->get_tag_values($tag);
                   $id .= 'mat_pept='.$id[0] .'|';
                   last;
             }
        }
        @tags = ('product'); # per Client request
        for my $tag (@tags) {
             if ($feat_obj->has_tag($tag)) {
                   @id = $feat_obj->get_tag_values($tag);
                   $id .= 'product='.$id[0];
                   last;
             }
        }
        $debug && print STDERR "$subname: feat_obj #$i id=$id\n";

        # get the translation, either from DNA or directly from the feature
        my $translation = '';
        if ( ! $feat_obj->has_tag('translation') ) {
           my $s;
           if ( 0 ) {
               $s = get_feature_nuc($feat_obj, $seq_obj->seq);
           } else {
               $s = Annotate_Math::get_dna_byloc( $feat_obj, $seq_obj->seq);
           }
           $debug && print STDERR "$subname: feat_obj #$i \n\$s=$s\n";

           my $f = Bio::PrimarySeq->new(-seq => $s, -desc => $id, -alphabet => 'dna');
           $f->revcom() if ($feat_obj->strand == -1);

           # only keep the mat_peptide that confirms to CDS
           $translation = $f->translate()->seq;
           $debug && print STDERR "$subname: \$f=\n".Dumper($f)."end of \$f\n\n";

           # add the translation to the feature
           $feat_obj->add_tag_value('translation', $f->translate()->seq) if ( 0 );
           $debug && print STDERR "$subname: #$i translation obtained from DNA\n";

        } else {
           print "$subname: WARNING: found sequence for mat_peptide from gb file.\n";
           my @s = $feat_obj->get_tag_values('translation');
           $translation = $s[0];
        }
        # Get translation for CDS
        my @parent_cds_seq = ();
        if ($parent_cds->has_tag('translation')) {
            @parent_cds_seq = $parent_cds->get_tag_values('translation');
        } else {
            @parent_cds_seq = (Annotate_Util::get_new_translation( $parent_cds, $parent_cds));
        }
        if ($translation =~ /[*]/) {
#            $count_err->{total}++;
#            $count_err->{stop_in_seq}++;
=head1
            print STDERR "$subname: ERROR: Found stop* in $acc #$i mat_peptide AA seq=".$feat_obj->location->to_FTstring."\n";
            print STDERR "$subname: ERROR: indicating problem in GBK file, such as AM408911\n";
            print STDERR "$subname: feat_obj #$i \n\$translation      =$translation\n";
            print STDERR "$subname: feat_obj #$i \n\$parent_cds_seq[0]=$parent_cds_seq[0]\n";
=cut
            if (0) {
                next;
            } else {
                # This section changes * to ., so that the translation from CDS can be taken as that of the new mat_peptide
                my $t = $translation;
                $t =~ s/[*]/./g;
                $translation = $t;
                $debug && print STDERR "$subname: feat_obj #$i \$t=$t\n";
                $debug && print STDERR "$subname: feat_obj #$i \$parent_cds_seq[0]=$parent_cds_seq[0]\n";
                if ($parent_cds_seq[0] =~ /($t)/) {
                    $debug && print STDERR "$subname: \$1=$1\n";
                    $translation = $1;
                }
            }
        }
        if ($translation && ($parent_cds_seq[0] !~ /$translation/i)) {
            $count_err->{total}++;
            $count_err->{mismatched_seq}++;
            print STDERR "$subname: WARNING: $acc translation for mat_peptide=".$feat_obj->location->to_FTstring." doesn't match CDS=".$parent_cds->location->to_FTstring.".\n";
            $debug && print STDOUT "$subname: WARNING: $acc translation for mat_peptide=".$feat_obj->location->to_FTstring." doesn't match CDS=".$parent_cds->location->to_FTstring.".\n";
            print STDERR "$subname: \$parent_cds & mat_peitide: $id\n";
            print STDERR "$subname: \$translation=$translation\n";
            next;
        }
        $count_err->{total}++;

        # see if MSA result contains same mat_peptide. If so, add MSA to source
        my $source = 'src=GBK';
        my $str1 = $parent_cds->location->start .'..'. $parent_cds->location->end;
        $debug && print STDERR "$subname: \$str1=$str1\n";
        # each CDS and corresponding mat_peptides
#        CDS: for (my $i=0; $i<=$#{$feats_msa}; $i++) {
#            my $feats = $feats_msa->[$i];
        CDS: for my $i (keys %$feats_msa) {
            my $feats = $feats_msa->{$i};
            $debug && print STDERR "$subname: \$acc=$acc \$feats=\n". Dumper($feats) . "End of \$feats\n\n";

            my $cds = $feats->[0];
            if (!$cds) {
                print STDERR "$subname: ERROR: \$acc=$acc NULL feature found at \$i=$i\n";
                next;
            }
            my $str2 = $cds->location->start .'..'. $cds->location->end;
            $debug && print STDERR "$subname: $acc \$str1=$str1 \$str2=$str2\n";
#            next if ($str1 ne $str2);

            foreach my $ifeat (0 .. $#{$feats}) {
                my $feat = $feats->[$ifeat];
                $debug && print STDERR "$subname: \$ifeat=$ifeat \$feat=".$feat->primary_tag.":".$feat->location->to_FTstring."\n";
                next if ($feat->primary_tag eq 'CDS'); # Exclude CDS
                # get the translation
                my @values = $feat->get_tag_values('translation');
                my $s = $values[0];
                $s =~ s/X/./i;
                $debug && print STDERR "$subname: $acc \$s          =$s\n";
                $debug && print STDERR "$subname: $acc \$translation=$translation\n";
                # if ($s eq $translation) {
                if ($translation =~ m/^$s$/i) {
                    $source = 'src=MSA,GBK';
                    $debug && print STDERR "$subname: $acc \$source=$source\n";

                    # update the source in note for the MSA feature
                    my $notes = [ $feat->get_tag_values('note') ];
                    for my $note (@$notes) {
                        next if ($note !~ /Desc:/i);
                        $note =~ s/src=MSA[|]/src=MSA,GBK|/i;
                        $debug && print STDERR "$subname: \$acc=$acc \$note=$note\n";
                        $feat->remove_tag('note');
                        $feat->add_tag_value('note', @$notes);
                    }
                    $debug && print STDERR "$subname: \$acc=$acc \$feat=\n". Dumper($feat) . "End of \$feat\n\n";
                    $debug && print STDERR "$subname: \$acc=$acc \$feats=\n". Dumper($feats) . "End of \$feats\n\n";

                    last CDS;
                } elsif ($feat->location->to_FTstring eq $feat_obj->location->to_FTstring) {
                    $debug && print STDERR "$subname: $acc         diff=".Annotate_Verify::diff_2str($s, $translation)."\n";
                }
            }
        }
        $id = "$source|".$id;

        my $f = Bio::PrimarySeq->new(
                                      -seq => $translation,
#                                      -id => $id, # id in fasta can't have whitespace, desc can
                                      -desc => $id, # id in fasta can't have whitespace, desc can
                                      -alphabet => 'protein'
                                    );
        push @records, $f;

    } # end of for my $feat_obj
    my $n = $count_err->{total};
    for my $key (keys %$count_err) {
    	$n -= $count_err->{$key} if ($key ne "total");
    }
    if ($count_err->{total}) {
        $comment = "In GBK, $count_err->{total} total";
        $comment .= "; $n good" if ($n);
    } else {
        $comment = "No mat_peptides in gbk";
    }
    if ($count_err->{mismatched_seq} || $count_err->{stop_in_seq} ||$count_err-> {null_CDS}) {
#        my $count_err = { total=>0, mismatched_seq=>0, stop_in_seq=>0, null_CDS=>0};
        my $comm = '';
        for my $k ('mismatched_seq','stop_in_seq','null_CDS') {
            next if (!$count_err->{$k});
            $comm .= ", " if ($comm);
            $comm .= "$count_err->{$k} $k" ;
        }
        $comment .= "; problem: $comm";
    }
    $debug && print STDERR "$subname: \@records=\n".Dumper(@records)."end of \@records\n\n";

    return (\@records, $comment);
} # sub extract_mature_peptides


## sub extract_mature_peptides0 searches a Bio::Seq::RichSeq object for features
## named mat_peptide, adds the translation tag to the feature if not present,
## and returns an array of Bio::PrimarySeq containing only id, alphabet, and seq.
## Backup version. Feb 14, 2011
sub extract_mature_peptides0 {
    my ($seq_obj, $feats_msa,$exe_dir) = @_;

    my $debug = 0 && $debug_all;
    my $subname = 'extract_mature_peptides0';

    my $comment = 'No CDS found in gbk';
    my @records = ();
    my ($parent_cds, $count);
    # Bio::Seq::RichSeq => Bio::SeqFeature::Generic
    my @feat_objs = $seq_obj->get_SeqFeatures;
    for my $i (0 .. $#feat_objs) {
        my $feat_obj = $feat_objs[$i];
        if ($feat_obj->primary_tag eq "CDS") {
            $comment = 'No mat_peptide found in gbk';
            $parent_cds = $feat_obj;
            $count = 0;
            next;
        }
        next if ($feat_obj->primary_tag ne "mat_peptide");
        my $n = $#records+1;
        $comment = "Found $n mat_peptides in gbk";
        $debug && print STDERR "$subname: feat_obj #$i is '". $feat_obj->primary_tag ."'\n";
        # following process prints the literal DNA sequence, w/o consideration of backups
#    print '$feat_obj->seq='.$feat_obj->seq->seq ."\n";
    # for the join(12332..12358,12358..15117) in DQ848678, this produceds garbage
#    print '$feat_obj->seq='.$feat_obj->seq->translate->seq ."\n";

# The id of mat_peptide is given as: first, either GI (db_xref) or NP (protein_id) of
# parent CDS or protein_id, followed by range of the mat_peptide, followed by its
# 'db_xref', 'protein_id', and 'product' in this order, whichever appears first.
           my (@id);
           my $id = 0;
           # get the id from parent CDS feature, 1) GI, 2) NP
           if ($parent_cds->has_tag('db_xref')) {
                @id = $parent_cds->get_tag_values('db_xref');
                for my $id1 (@id) {
                   if ($id1 =~ /^GI:/i) {
                	$id = 'CDS='. $id1 ."|";
                	last;
                   }
                }
           }
           if (($id !~ /^CDS=GI:/) && $parent_cds->has_tag('protein_id')) {
                @id = $parent_cds->get_tag_values('protein_id');
                $id = 'CDS='. $id[0] ."|";
           }

           # In case $id is still empty
           if (($id !~ /^CDS=/) ) {
             my @tags = ('db_xref', 'protein_id', 'product'); # per Chris request
             for my $tag (@tags) {
                if ($parent_cds->has_tag($tag)) {
                   @id = $feat_obj->get_tag_values($tag);
                   $id .= 'CDS='. $id[0] ."|";
                   last;
                }
             }
           }
           $id = "CDS=unknown|" if (($id !~ /^CDS=/) ); # indicates a potential problem
           $id = "ACC=". $parent_cds->seq->accession_number .'.'.$parent_cds->entire_seq->{'_version'}.'|' .$id;
           
           my $source = 'src=GBK';
            for (my $i=0; $i<=$#{$feats_msa}; $i++) {
                my $feats = $feats_msa->[$i];
                $debug && print STDERR "$subname: \$feats=\n". Dumper($feats) . "End of \$feats\n\n";

                my $inseq = $feats->[0]->seq;
                # either print to STDERR or fasta file
                my $acc = $inseq->accession_number;
#                my $outfile = $acc . '_matpept_muscle.faa' if (!$debug);
#                my $faa1 .= Annotate_misc::generate_fasta( $feats, $outfile, '');
#                print STDERR "$subname: accession='$acc' CDS=".$feats->[0]->location->to_FTstring."\n";
#                $debug && print STDERR "$subname: \$faa1 = '\n$faa1'\n";

#                Annotate_Verify::check_old_annotation( $acc, $faa1);
                print STDERR "$subname: accession = '".$acc."'\n";
            }
           
           $id = "$source|".$id;

           # add the range of mat_peptide
           my $s = '';
           $s = Annotate_Util::get_DNA_loc( $feat_obj);
#           $id .= 'Loc='. $feat_obj->location->to_FTstring ."|";
           $id .= 'Loc='. $s ."|";

           $s = Annotate_Util::get_AA_loc( $feat_obj, $parent_cds);
           $id .= 'AA='. $s .'|';

           # add gene_symbol
           my $symbol = Annotate_Util::get_gene_symbol( $feat_obj,$exe_dir);
           $id .= 'symbol='. $symbol ."|";

           # add the id of mat_peptide
           my @tags = ('db_xref', 'protein_id', 'product'); # per Chris request
           for my $tag (@tags) {
                if ($feat_obj->has_tag($tag)) {
                   @id = $feat_obj->get_tag_values($tag);
                   $id .= 'mat_pept='.$id[0];
                   last;
                }
           }
           $debug && print STDERR "$subname: feat_obj #$i id=$id\n";

        # get the translation, either from DNA or directly from the feature
        if ( ! $feat_obj->has_tag('translation') ) {
           my $s = get_feature_nuc($feat_obj, $seq_obj->seq);
           $debug && print STDERR "$subname: feat_obj #$i \n\$s=$s\n";

           my $f = Bio::PrimarySeq->new(-seq => $s, -desc => $id, -alphabet => 'dna');
           $f->revcom() if ($feat_obj->strand == -1);

           # only keep the mat_peptide that confirms to CDS
           $s = $f->translate()->seq;
           my @parent_cds_seq = $parent_cds->get_tag_values('translation');
           if ($parent_cds_seq[0] !~ /$s/i) {
                print STDERR "WARNING: translation for mat_peptide doesn't match CDS.\n";
                print STDOUT "WARNING: translation for mat_peptide doesn't match CDS.\n";
                print STDERR "\$parent_cds & mat_peitide: $id\n";
                print STDERR "\$s=$s\n";
                next;
           }
           push @records, $f;
           $debug && print STDERR "$subname: \$f=\n".Dumper($f)."end of \$f\n\n";

           # add the translation to the feature
           $feat_obj->add_tag_value('translation', $f->translate()->seq) if ( 0 );

        } else {
           print "WARNING: gbk_matpeptide.pl found sequence for mat_peptide from gb file.\n";
           my @s = $feat_obj->get_tag_values('translation');
#          print "s[0]=$s[0]\n";
           my $f = Bio::PrimarySeq->new(-seq => $s[0], -id => $id, -alphabet => 'protein');
           push @records, $f;
        }
   } # end of for my $feat_obj
   $debug && print STDERR "$subname: \@records=\n".Dumper(@records)."end of \@records\n\n";

   return (\@records, $comment);
} # sub extract_mature_peptides0


## sub get_feature_nuc fetches the sequence directly from seq string
## according to the location range
sub get_feature_nuc {
    my ($feat_obj,$seq_obj_seq) = @_;

    my $debug = 0 && $debug_all;
    my $subname = 'get_feature_nuc';

    $debug && print STDERR "$subname: \$feat_obj=\n".Dumper($feat_obj)."end of \$feat_obj\n\n";
    my $s = '';
      if ($feat_obj->location->isa('Bio::Location::Simple')) {
#	   print "  Simple:DNA    :seq=". $feat_obj->seq->seq ."\n";
	   $s = $feat_obj->seq->seq;
              $debug && print STDERR "$subname: Simple \n\$s=$s\n";
#	   print "  Simple:protein:seq=". $feat_obj->seq->translate()->seq ."\n";
      } elsif ( $feat_obj->location->isa('Bio::Location::Fuzzy')) {
#	   print "  Simple:DNA    :seq=". $feat_obj->seq->seq ."\n";
	   $s = $feat_obj->seq->seq;
              $debug && print STDERR "$subname: Fuzzy  \n\$s=$s\n";
#	   print "  Simple:protein:seq=". $feat_obj->seq->translate()->seq ."\n";
      } elsif ($feat_obj->location->isa('Bio::Location::Split')) {
	   my $s0 = $seq_obj_seq;
#	   print "location=". $feat_obj->location->to_FTstring . "\n";
	   for my $loc ($feat_obj->location->sub_Location) {
	      $s .= substr($s0, $loc->start -1, $loc->end +1 - $loc->start);
#	      print "start=". $loc->start ."\tend=". $loc->end ."\n";
              $debug && print STDERR "$subname: Split  \n\$s=$s\n";
	   }
      }
    $debug && print STDERR "$subname: \n\$s=$s\n";

    return $s;
} # sub get_feature_nuc


=head2
 sub get_matpeptide takes a genome genbank file, searches for the mat_peptide in it
   returns the mat_peptide sequences to a FASTA string.
=cut

sub get_matpeptide {
#    my ($gbk, $feats_msa,$exe_dir) = @_;
    my ($seq_obj, $feats_msa,$exe_dir) = @_;

    my $debug = 0 && $debug_all;
    my $subname = 'get_matpeptide';

    my $comment = '';
=head2
    my $in_file2 = IO::String->new($gbk);
    my $seqio_obj = Bio::SeqIO->new( -fh => $in_file2, -format => 'genbank' );
#   my $seqio_obj = Bio::SeqIO->new(-file=> $infile);
   my $seq_obj = $seqio_obj->next_seq;
=cut

   my @records = ();
    # Only take 1st sequence from each gbk (Note: gbk can hold multiple sequences, we ignore all after 1st)
     $debug && print STDERR "$subname: \$seq_obj=\n".Dumper($seq_obj)."end of \$seq_obj\n\n";

     my $r1;
     # Gets an array of Bio::PrimarySeq objects, directly from genbank file
     ($r1, $comment) = extract_mature_peptides($seq_obj, $feats_msa,$exe_dir);
     $debug && print STDERR "$subname: \$r1=\n".Dumper($r1)."end of \$r1\n\n";
     push @records, @$r1;
     $debug && print STDERR "$subname: sub extract_mature_peptides returned $#{$r1} mat_peptides\n";

   # write out the entire genbank file
     if ( 0 ) {
	my $outgbfile = $seq_obj->accession_number."_matpeptide".$2;
        my $seq_out = Bio::SeqIO->new('-file' => ">$outgbfile", '-format' => 'genbank');
        $seq_out->write_seq($seq_obj);
     }

#    if (my $seq_obj1 = $seqio_obj->next_seq) {
#        print STDERR "$subname: There are additional sequences (".$seq_obj1->accession_number.") in input file, ignored\n";
#    }

        # write all mat_peptides to fasta file
        my $faa_gbk = '';
        if (@records) {
          my $vfile = IO::String->new($faa_gbk);
          my $seq_out = Bio::SeqIO->new(
                  '-fh'     => $vfile,
                  '-format' => 'fasta'
                              );
          my $i = 0;
          for my $seq (@records) {
	     $i++;
	     if ($seq->alphabet eq 'dna') {
		$seq_out->write_seq($seq->translate());
	     } elsif($seq->alphabet eq 'protein') {
		$seq_out->write_seq($seq);
	     }
          }
        }
       $debug && print STDERR "$subname: \$faa_gbk=\n". Dumper($faa_gbk) . "End of \$faa_gbk\n\n";

   if ($seq_obj->accession_number !~ /^NC_/i) {
#      $faa_gbk = '';
#      $comment = "Not refseq, skip";
   }

   return ($faa_gbk, $comment);
} # sub get_matpeptide


1;
