package Annotate_Def;

use strict;
use warnings;
use English;
use Carp;
use Data::Dumper;

use LWP::Simple;
use File::Temp qw/ tempfile tempdir /;
use Bio::SeqIO;
use Bio::Seq;
use Bio::AlignIO;
use Bio::Tools::Run::StandAloneBlast;
use IO::String;

use version; our $VERSION = qv('1.1.7'); # Feb 05 2013

my $debug_all = 1;

####//README//####
#
# Annotate_Def contains reference data used in the msa_annotate.pl script
#
#    Authors Chris Larsen, clarsen@vecna.com; Guangyu Sun, gsun@vecna.com
#    November 2012
#
##################

## //Global variables// ##

#loadTaxonTable: $TAXON=
#$VAR1 = {
#          'taxon_fn' => 'Annotate_taxon_records.txt',
#          'taxon' => {
#             '336486' => ['336486', '336486', '10260', '10240', 'Turkeypox virus', 'Avipoxvirus', 'Poxviridae'],
#             '463719' => ['463719', '357231', '142786', '11974', 'Murine norovirus', 'Norovirus', 'Caliciviridae'],
#             '389230' => ['389230', '389230', '694002', '11118', 'Bat coronavirus (BtCoV/133/2005)', 'Betacoronavirus', 'Coronaviridae'],
#                     },
#          'taxon_loaded' => 1
#        };

my $TAXON = {
               taxon_loaded => 0,
               taxon_fn     => "Annotate_taxon_records.txt",
            };

#get_refseq_acc: $REFSEQS=
#$VAR1 = {
#          'refseq_loaded' => 1,
#          'refs' => {
#                      'Arenaviridae' => {
#                                          '11619' => 'NC_005081'
#                                        },
#                      'Flaviviridae' => {
#                                          '11080' => 'NC_007580',
#                                          '31655' => 'NC_004102',
#                                          '11070' => 'NC_002640',
#                                        },
#                    }
#        };
my $REFSEQS = {
               refseq_loaded => 0,
               refs => {},
            };

#get_gene_symbol: $GENE_SYM=
#$VAR1 = {
#          'symbol' => {
#                        'NC_004718' => {
#                                         '8485..9984' => 'nsp4',
#                                         '13372..13410' => 'nsp11',
#                                         '17970..19550' => 'nsp14',
#                                         '805..2718' => 'nsp2',
#                                       }
#                      },
#          'symbol_loaded' => 1,
#          'symbol_fn' => 'Annotate_symbol_records.txt'
#        };
my $GENE_SYM = {
               symbol_loaded => 0,
               symbol_fn     => "Annotate_symbol_records.txt",
            };


## //subroutines// ##

=head2
 sub downloadRefseq searches genbank for all refseqs for a family, then downloads all the refseqs
=cut
=head2
sub downloadRefseq {
    my ($fam, $exe_dir, $download_REFSEQ) = @_;

    my $debug = 0 && $debug_all;
    my $subn = 'downloadRefseq';

    my $download_SUMMARY = $download_REFSEQ;
#    $download_REFSEQ = 0;

    $exe_dir = './' if (!$exe_dir);
    $debug && print STDERR "$subn: \$exe_dir='$exe_dir'\n";

    # Search for NC_* for each family from genbank
    my $count = 0;
    my $web = '';
    my $key = '';
    my $base = 'http://eutils.ncbi.nlm.nih.gov/entrez/eutils/';
#    my $outFileName = undef;
    my $outFileName = '';
    $debug && print STDERR "\n$subn: \$fam=$fam\n";

    if ( $download_SUMMARY ) {
        # modeled after http://www.ncbi.nlm.nih.gov/books/NBK25498/#chapter3.Application_3_Retrieving_large -11/28/2012

        #assemble the esearch URL
        my $query = "${fam}[orgn]+AND+srcdb_refseq[prop]";
        $debug && print STDERR "\n$subn: \$fam=$fam \$query='$query'\n";
        my $url = $base . "esearch.fcgi?db=nucleotide&term=$query&usehistory=y&retmax=100";
        $debug && print STDERR "$subn: \$fam=$fam \$url=$url\n";

        # post the esearch URL
        my $output = '';
        for my $ntry (0 .. 2) {
            sleep(10);
            $output = get($url);
            $debug && print STDERR "$subn: \$fam=$fam \$output='$output'\n";
            last if ($output !~ /Unable to obtain query/);
        }
        return $outFileName if ($output =~ /Unable to obtain query/);
        # parse WebEnv, QueryKey and Count (# records retrieved)
        $web = $1 if ($output =~ /<WebEnv>(\S+)<\/WebEnv>/);
        $key = $1 if ($output =~ /<QueryKey>(\d+)<\/QueryKey>/);
        $count = $1 if ($output =~ /<Count>(\d+)<\/Count>/);
        my $time = POSIX::strftime("%m/%d/%Y %H:%M:%S", localtime);
        print STDERR "$subn: \$fam=$fam \$count=$count \$web=$web \$key=$key at $time\n";
    }

    # Have to get the summary before downloading gbk files
    $outFileName = "RefSeq_${fam}.gbk";
    if ( $download_REFSEQ  && $download_SUMMARY) {
        # open output file for writing
        my $bakFileName = "$outFileName.bak";
        if (-e "$exe_dir/$outFileName") {
                for my $i (reverse 1 .. 8) {
                    $bakFileName = sprintf("$outFileName.bak$i");
                    my $bakFileName1 = sprintf("$outFileName.bak%d", $i+1);
                    if (-e "$exe_dir/$bakFileName") {
                        `mv $exe_dir/$bakFileName $exe_dir/$bakFileName1`;
                    }
                }
                my $result = `mv $exe_dir/$outFileName $exe_dir/$bakFileName`;
                print STDERR "$subn: \$fam=$fam existing file $outFileName moved to $bakFileName \$result=$result\n";
        }
        open (my $OUT, ">$exe_dir/$outFileName") || die "$subn: Can't open file'$outFileName'!\n";
        # Download the NC_* as a set
        # retrieve data in batches of 100
        my $retmax = 100;
        for (my $retstart = 0; $retstart < $count; $retstart += $retmax) {
            $debug && print STDERR "$subn: \$fam=$fam Downloading \$retstart=$retstart\n";
            ($retstart>0) && sleep(10);
            my $efetch_url = $base;
            $efetch_url .= "efetch.fcgi?";
            $efetch_url .= "db=nucleotide";
            $efetch_url .= "&WebEnv=$web";
            $efetch_url .= "&query_key=$key";
            $efetch_url .= "&retstart=$retstart";
            $efetch_url .= "&retmax=$retmax";
            $efetch_url .= "&rettype=gb";
            $efetch_url .= "&retmode=text";
            my $efetch_out = '';
            for my $ntry (0 .. 2) {
                $efetch_out = get($efetch_url);
                last if ($efetch_out !~ /Unable to obtain query/);
            }
            return undef if ($efetch_out =~ /Unable to obtain query/);
            print $OUT "$efetch_out";
        }
        close $OUT;
        $debug && print STDERR "$subn: Downloaded $count genomes for $fam, saved to $exe_dir/$outFileName\n";
        # check if there is any difference between new and existing genbank file
        if (-e "$bakFileName") {
              my $result = `diff $exe_dir/$bakFileName $exe_dir/$outFileName`;
              if ($result) {
                print STDERR "$subn: ERROR: \$fam=$fam difference between $bakFileName and $outFileName \$result='$result'\n";
              } else {
                print STDERR "$subn: \$fam=$fam identical $bakFileName and $outFileName \$result='$result'\n";
              }
        }
    }
    if (!-e "$exe_dir/$outFileName") {
        $outFileName = '';
    }

    return $outFileName;
} # sub downloadRefseq
=cut

=head2
 sub downloadTaxon for a family downloads the the taxonnomy information for all strains/species
=cut

=head2
sub downloadTaxon {
    my ( $fam, $exe_dir, $download_TAXON) = @_;

    my $debug = 0 && $debug_all;
    my $subn = 'downloadTaxon';

    my $count = 0;
    $exe_dir = './' if (!$exe_dir);
    $debug && print STDERR "$subn: \$exe_dir='$exe_dir'\n";

    # Search for NC_* for each family from genbank
    my $web = '';
    my $key = '';
    my $base = 'http://eutils.ncbi.nlm.nih.gov/entrez/eutils/';
    my $outFileName = '';
#    $outFileName = "RefSeq_${fam}.gbk";
#    print STDERR "\n$subn: \$fam=$fam \$outFileName=$outFileName\n";

    my $taxonFileName = "taxon_${fam}.xml";
    $download_TAXON = 0 if (!$download_TAXON);
    if ( $download_TAXON ) {
        #assemble the esearch URL
        my $query = "${fam}[orgn]";
        print STDERR "\n$subn: \$fam=$fam \$outFileName=$outFileName \$query='$query'\n";
        my $url = $base . "esearch.fcgi?db=taxonomy&term=$query&usehistory=y&retmax=10&rettype=xml";
        $debug && print STDERR "$subn: \$fam=$fam \$url='$url'\n";

        my $delay = 1;
        # post the esearch URL
        sleep( $delay );
        my $output = get($url);
        $debug && print STDERR "$subn: \$fam=$fam \$output='$output'\n";
        # parse WebEnv, QueryKey and Count (# records retrieved)
        $web = $1 if ($output =~ /<WebEnv>(\S+)<\/WebEnv>/);
        $key = $1 if ($output =~ /<QueryKey>(\d+)<\/QueryKey>/);
        $count = $1 if ($output =~ /<Count>(\d+)<\/Count>/);
        my $time = POSIX::strftime("%m/%d/%Y %H:%M:%S", localtime);
        print STDERR "$subn: \$fam=$fam \$count=$count \$web=$web \$key=$key at $time\n";

        # Download the TaxonSet for each id, from efetch.fcgi
        open (my $OUT, ">$exe_dir/$taxonFileName") || die "$subn: Can't open file='$taxonFileName'!\n";
#        for (my $i=0; $i<$count; $i++) {
        my $retmax = 5000;
        # download the taxon
        for (my $retstart = 0; $retstart < $count; $retstart += $retmax) {
            $debug && print STDERR "$subn: \$fam=$fam \$count=$count Downloading \$retstart=$retstart\n";
            sleep( $delay );
            my $efetch_url = $base;
            $efetch_url .= "efetch.fcgi?";
            $efetch_url .= "db=taxonomy";
            $efetch_url .= "&WebEnv=$web";
            $efetch_url .= "&query_key=$key";
            $efetch_url .= "&retstart=$retstart";
            $efetch_url .= "&retmax=$retmax";
            $efetch_url .= "&rettype=xml";
            $efetch_url .= "&retmode=xml";
            my $efetch_out = '';
            for my $ntry (0 .. 2) {
                ($ntry>0) && sleep( $delay );
                $efetch_out = get($efetch_url);
                print STDERR "$subn: \$ntry=$ntry \$efetch_out=".substr($efetch_out, 0, 300)."\n";
                last if ($efetch_out !~ /Unable to obtain query/);
            }
            if ($efetch_out =~ /Unable to obtain query/) {
                print STDERR "$subn: ERROR: \$count=$count \$retstart=$retstart \$efetch_url='$efetch_url'\n";
                print STDERR "$subn: \$efetch_out=$efetch_out\n";
                return undef ;
            } else {
                print $OUT "$efetch_out";
            }
        }

        close $OUT;
        $debug && print STDERR "$subn: Downloaded $count genomes for $fam, saved to $exe_dir/$taxonFileName\n";
    }
    if (!-e "$exe_dir/$taxonFileName") {
        print STDERR "$subn: ERROR: \$fam=$fam file '$exe_dir/$taxonFileName' non-existant, and no live download performed\n";
        $taxonFileName = '';
    }

    return $taxonFileName;
} # sub downloadTaxon
=cut

=head2
 sub loadXmlTaxon reads an XML file, searchs for each taxid, and its name, speciesid, species,
 genus, genusid, family, and familyid, to be compared with exsting data
=cut
=head2

sub loadXmlTaxon {
    my ($fam, $taxonFileName, $exe_dir) = @_;

    my $debug = 0 && $debug_all;
    my $subn = 'loadXmlTaxon';

    my $newTax = {};
    $exe_dir = './' if (!$exe_dir);
    $debug && print STDERR "$subn: \$exe_dir='$exe_dir'\n";
    if (!-e "$exe_dir/$taxonFileName") {
        print STDERR "$subn: ERROR: \$fam=$fam file '$exe_dir/$taxonFileName' non-existant, and no live download performed\n";
        return $newTax;
    }

    open (my $taxFile, "<$exe_dir/$taxonFileName") || die "$subn: Can't open file='$taxonFileName'!\n";

    my @set = ();
    while (<$taxFile>) {
        chomp;
        next if ($_ !~ /^(<TaxaSet><Taxon>|<Taxon>)/);
        # First, identify the <Taxon></Taxon>
        $_ =~ s/(<TaxaSet>|<\/TaxaSet>)//;
        @set = ();
        push @set, $_;
        while (<$taxFile>) {
            chomp;
            $_ =~ s/(<TaxaSet>|<\/TaxaSet>)//;
            push @set, $_;
            $debug && print STDERR "$subn: \$_='$_'\n";
            last if ($_ =~ /^<\/Taxon>/);
        }
        $debug && print STDERR "$subn: \@set=".Dumper(@set)."\n";

        # Parse the <Taxon></Taxon>
        # turn xml to hash
        my $taxid = '';
        my $tax = xml2hash( @set);

        my $taxon = getTaxonString( $tax);
        my $t = [
                  $taxon->{'TaxId'},
                  $taxon->{'speciesid'},
                  $taxon->{'genusid'},
                  $taxon->{'familyid'},
                  $taxon->{'species'},
                  $taxon->{'genus'},
                  $taxon->{'family'},
                ];

        $newTax->{$taxon->{'TaxId'}} = $t;
        $debug && print STDERR "$subn: \$tax=\n".Dumper($tax)."\n";

    }
    $debug && print STDERR "$subn: \$newTax=\n".Dumper($newTax)."\n";

    return $newTax;
} # sub loadXmlTaxon
=cut

=head2
 sub getTaxonString takes a hash of a taxon, and finds the relative info, and returns an array
=cut

=head2
sub getTaxonString {
    my ($tax) = @_;

    my $debug = 0 && $debug_all;
    my $subn = 'getTaxonString';

    $debug && print STDERR "$subn: Just got in the sub: \$tax=\n".Dumper($tax)."\n";
    my $taxInfo = {
                    'TaxId' => -1,
                    'speciesid' => -1,
                    'species' => '',
                    'genusid' => -1,
                    'genus' => '',
                    'familyid' => -1,
                    'family' => '',
                  };
    my @starts = ();
    my @tags = ();
    @tags = keys %{$tax->{'Taxon'}};

    $taxInfo->{'TaxId'} = $tax->{'Taxon'}->{'TaxId'};
    if ($tax->{'Taxon'}->{'Rank'} eq 'species') {
        $taxInfo->{'speciesid'} = $tax->{'Taxon'}->{'TaxId'};
        $taxInfo->{'species'} = $tax->{'Taxon'}->{'ScientificName'};
    }
    if ($tax->{'Taxon'}->{'Rank'} eq 'genus') {
            $taxInfo->{'genusid'} = $tax->{'Taxon'}->{'TaxId'};
            $taxInfo->{'genus'} = $tax->{'Taxon'}->{'ScientificName'};
    }
    if ($tax->{'Taxon'}->{'Rank'} eq 'family') {
            $taxInfo->{'familyid'} = $tax->{'Taxon'}->{'TaxId'};
            $taxInfo->{'family'} = $tax->{'Taxon'}->{'ScientificName'};
    }

    my $ts = $tax->{'Taxon'}->{'LineageEx'}->{'Taxon'};
    for my $t (@$ts) {
        if ($t->{'Rank'} eq 'species') {
            $taxInfo->{'speciesid'} = $t->{'TaxId'};
            $taxInfo->{'species'} = $t->{'ScientificName'};
        }
        if ($t->{'Rank'} eq 'genus') {
            $taxInfo->{'genusid'} = $t->{'TaxId'};
            $taxInfo->{'genus'} = $t->{'ScientificName'};
        }
        if ($t->{'Rank'} eq 'family') {
            $taxInfo->{'familyid'} = $t->{'TaxId'};
            $taxInfo->{'family'} = $t->{'ScientificName'};
        }
    }

    $debug && print STDERR "$subn: right before returning from subroutine \$taxInfo=\n".Dumper($taxInfo)."\n";

    return $taxInfo;
} # sub getTaxonString
=cut

=head2
 sub xml2hash takes an array of strings of XML, turns it to a hash
=cut

=head2
sub xml2hash {
    my (@set) = @_;

    my $debug = 0 && $debug_all;
    my $subn = 'xml2hash';

    $debug && print STDERR "$subn: Just got in the sub: \@set=\n".Dumper(@set)."\n";
    my $tax = {};
    my @starts = ();
    my $start = -1;
    my @ends = ();
    my @tags = ();

    for (my $i=0; $i<=$#set; $i++) {
        my $line = $set[$i];
        $debug && print STDERR "$subn: \$i=$i \$line='$line'\n";
        my $tag = '';
#        $tag = '' if (!$tag);
        $debug && print STDERR "$subn: \$tag='$tag' \$start=$start \@tags='@tags' \$line='$line'\n";
        if ($#tags<0) {
            if ($line =~ /<([^\/]+)>/) {
                $tag = $1;
                if ($line =~ /<$tag>(.*)<\/$tag/) {
                    $set[$i] =~ s/<$tag>(.*)<\/$tag//;
                    $tax->{$tag} = $1;
                    $debug && print STDERR "$subn: 1. \$1=$1 \$start=$start \@tags='@tags' \$tax=\n".Dumper($tax)."\n";
                    $line =~ s/<$tag>$tax->{$tag}<\/$tag>//;
                    next;
                }
                push @tags, $tag;
                $set[$i] =~ s/<([^\/]+)>//;
                $start = $i;
                next;
            }
            $debug && print STDERR "$subn: 2. \$start=$start \@tags='@tags' \$tax=\n".Dumper($tax)."\n";
        } else {
            if ($line =~ /<([^\/]+)>/) {
                $tag = $1;
                if ($line =~ /<$tag>(.*)<\/$tag/) {
                    next;
                }
                push @tags, $tag;
                next;
            } elsif ($line =~ /<\/$tags[$#tags]>/) {
                my @set1 = ( @set[$start .. $i] );
                $debug && print STDERR "$subn: 3. Found last tags: \$start=$start \@tags='@tags'  \$set1=\n".Dumper(@set1)."\n";
                $tag = pop @tags;
                if ($#tags<0) {
                    $set[$i] =~ s/<\/$tag>//;
                    my $t = xml2hash( @set[$start .. $i]);
                    for my $k ($start .. $i) { $set[$k] = ''; } # Clear up those lines already processed
                    if (exists($tax->{$tag}) && ref($tax->{$tag}) eq 'HASH') {
                        $tax->{$tag} = [$tax->{$tag}, $t];
                    } elsif (exists($tax->{$tag}) && ref($tax->{$tag}) eq 'ARRAY') {
                        push @{$tax->{$tag}}, $t;
                    } else {
                        $tax->{$tag} = $t;
                    }
                }
                $debug && print STDERR "$subn: 4. \$tag=$tag \$start=$start \@tags='@tags' \$tax->{$tag}=".ref($tax->{$tag})."\n";
            } else {
                $debug && print STDERR "$subn: 1. \$i=$i problem with line='$line'\n";
            }
        }
        $debug && print STDERR "$subn: \$i=$i \$tax=\n".Dumper($tax)."\n";
    }

    $debug && print STDERR "$subn: right before returning from subroutine \$tax=\n".Dumper($tax)."\n";
    return $tax;
} # sub xml2hash
=cut

=head2
 sub checkAllTaxon looks over each family and see if there is any update in its taxonomy from NCBI
=cut

=head2
sub checkAllTaxon {
    my ( $exe_dir) = @_;

    my $debug = 0 && $debug_all;
    my $subn = 'checkAllTaxon';

    # Load the $REFSEQS if not already
    if ( !$TAXON->{taxon_loaded} ) {
        Annotate_Def::loadTaxonTable( $exe_dir);
        my $ctLines = scalar(@{$TAXON->{taxon}});
        $debug && print STDERR "$subn: loaded $ctLines lines from file '$exe_dir/$TAXON->{taxon_fn}'\n";
    }
    $debug && print STDERR "$subn: \$TAXON=\n".Dumper($TAXON)."end of \$TAXON\n\n";
    if ( !$REFSEQS->{refseq_loaded} ) {
        my $nstrain = Annotate_Def::initRefseq();
        my @fam = sort keys %{$REFSEQS->{refs}};
        printf STDERR ("$subn: loaded RefSeqs for $nstrain strains in %d families: @fam\n", $#fam+1);
    }
    $debug && print STDERR "$subn: \$REFSEQS=\n".Dumper($REFSEQS)."end of \$REFSEQS\n\n";

    my $count = 0;
    $exe_dir = './' if (!$exe_dir);
    $debug && print STDERR "$subn: \$exe_dir='$exe_dir'\n";

    my $taxonFileName = "Annotate_taxon_records.txt"; # plain text file storing the taxon info for all families
    my $taxonTempName = "$taxonFileName.0";
    my $taxonFile = undef;
    open($taxonFile, ">$exe_dir/$taxonTempName") || die "$subn: Can't open file='$taxonTempName'!\n";

    # Check strains for each family from genbank
    my $DOWNLOAD_TAXON = 0;
    my $genera = {};
    my $TAXON_UPDATED = 0;
    for my $fam (sort keys %{$REFSEQS->{refs}}){
        my $web = '';
        my $key = '';
        my $base = 'http://eutils.ncbi.nlm.nih.gov/entrez/eutils/';
        my $outFileName = '';

        $outFileName = Annotate_Def::downloadTaxon( $fam, $exe_dir, $DOWNLOAD_TAXON);
        $outFileName = Annotate_Def::downloadTaxon( $fam, $exe_dir, $DOWNLOAD_TAXON) if (!$outFileName); # Effectively try again to download
        if ($DOWNLOAD_TAXON) {
            print STDERR "$subn: \$fam=$fam, performed live download\n";
        } else {
            print STDERR "$subn: \$fam=$fam, existing file used, no live download\n";
        }
        print STDERR "$subn: \$outFileName=$outFileName\n";

        if (!$outFileName) {
            print STDERR "$subn: \$fam=$fam\n";
            print STDERR "$subn: ERROR: \$fam=$fam \$taxonFileName='$outFileName' is empty.\n";
            print STDERR "$subn: \$fam=$fam, perhaps turn on the download options in sub Annotate_Def::downloadTaxon\n\n";
            next;
        }

        my $newTax = loadXmlTaxon( $fam, $outFileName, $exe_dir);
        $debug && print STDERR "$subn: \$newTax=".scalar(keys %$newTax)."\n".Dumper($newTax)."end of \$newTax\n\n";

        my $status = { total => 0, strain => 0, identical => 0, new => 0, };
        my $msg = '';
        my $str = sprintf("# Family=$fam taxon loaded from file=$outFileName \n");
        my $time = POSIX::strftime("%m/%d/%Y %H:%M:%S", localtime);
        $str .= sprintf("# Saved from Annotate_Dev.pm Ver$VERSION on %s \n", $time);
        my $str1 = '';
        my $ct = -1;
        for my $taxid (sort {$newTax->{$a}->[2]<=>$newTax->{$b}->[2] || $newTax->{$a}->[1]<=>$newTax->{$b}->[1]} keys %$newTax) {
            next if ($newTax->{$taxid}->[1]<0);
            $ct++;
            $status->{total}++;
            my $newStr = '';
            $newStr = sprintf("| %7d | %7d |", $newTax->{$taxid}->[0], $newTax->{$taxid}->[1]);
            $newStr .= sprintf(" %7d | %7d |", $newTax->{$taxid}->[2], $newTax->{$taxid}->[3]);
            $newStr .= sprintf(" %-28s |", $newTax->{$taxid}->[4]);
            $genera->{$newTax->{$taxid}->[6]}->{$newTax->{$taxid}->[5]}->{$newTax->{$taxid}->[4]}++;
            $newStr .= sprintf(" %-16s |", $newTax->{$taxid}->[5]);
            $genera->{$newTax->{$taxid}->[6]}->{$newTax->{$taxid}->[5]}->{total}++;
            $newStr .= sprintf(" %-15s |", $newTax->{$taxid}->[6]);
            $str .= $newStr . "\n";
            $debug && print STDERR "$subn: \$ct=$ct \$taxid=$taxid \$newStr='$newStr'\n";
#            if ($newTax->{$taxid}->[1] == -1) { # speciesid
#                $msg .= "#$ct\t---\t$newStr\n";
#                next;
#            }
            my $oldStr = "";
            if (exists($TAXON->{taxon}->{$taxid})) {
                my $oldTax = $TAXON->{taxon}->{$taxid};
                $oldStr = sprintf("| %7d | %7d |", $oldTax->[0], $oldTax->[1]);
                $oldStr .= sprintf(" %7d | %7d |", $oldTax->[2], $oldTax->[3]);
                $oldStr .= sprintf(" %-28s |", $oldTax->[4]);
                $oldStr .= sprintf(" %-16s |", $oldTax->[5]);
                $oldStr .= sprintf(" %-15s |", $oldTax->[6]);
            }
            
            $debug && print STDERR "$subn: \$taxid=$taxid\n";
            $status->{strain}++;
            if ($newStr eq $oldStr) {
                $status->{identical}++;
                $msg .= "#$ct\t---\t$newStr\n";
                $debug && print STDERR "$subn: \$ct=$ct \$taxid=$taxid identical Taxon info found. Skip...\n";
            } else {
                $TAXON_UPDATED = 1;
                $status->{new}++;
                $msg .= "#$ct\tnew\t$newStr\n";
                $debug && print STDERR "$subn: \$ct=$ct \$taxid=$taxid found new Taxon info from NCBI download for $fam:\n";
            }
            $debug && print STDERR "$subn: \$ct=$ct \$taxid=$taxid \$oldStr='$oldStr'\n";
            $debug && print STDERR "$subn: \$ct=$ct \$taxid=$taxid \$newStr='$newStr' \$TAXON_UPDATED=$TAXON_UPDATED\n";
        } # for my $taxid (sort {}) {

        my $msg1 = "$fam found taxon from NCBI: $status->{total} total:";
        $msg1 .= " $status->{strain} strains,";
        $msg1 .= " $status->{identical} identical,";
        $msg1 .= " $status->{new} new\n";
        $msg = $msg1 . $msg;
        # Save to Anntate_taxon_records.txt
        $str .= sprintf("# Family=$fam taxon loaded from file=$outFileName, total=%d\n", $status->{total});
        $count += $status->{total};
        print $taxonFile "$str\n" || die "$subn: Can't write to file='$taxonFileName'!\n";
#        print STDERR "$subn: \$str='\n$str'\n";
        print STDERR "$subn: \$msg='\n$msg'\n";

    } # for my $fam (sort keys %{$REFSEQS->{refs}}){
    my @fams = ( sort keys %{$REFSEQS->{refs}} );
    printf $taxonFile ("# Updated $count strains in %d families: @fams\n", scalar(@fams));
    close $taxonFile || die "$subn: Can't close file='$taxonFileName'!\n";

if ( 0 ) {
    $debug && print STDERR "$subn: \$genera=\n".Dumper($genera)."\n";
    for my $fam (sort keys %$genera) {
      my $family = $genera->{$fam};
      print STDERR "\n$subn: family=$fam\n";
      for my $gen (sort keys %$family) {
        my $genus = $family->{$gen};
        print STDERR "$fam  '$gen'\n";
        printf STDERR ("$fam  %-12s  total  %3d\n", "'".$gen."'", $genus->{total});
        for my $spe (sort keys %$genus) {
          next if ($spe eq 'total');
          printf STDERR ("$fam  %-12s  %4d  '$spe'\n", "'".$gen."'", $genus->{$spe});
        }
      }
      print STDERR "\n";
    }
}

    # save and report the result
    printf STDERR ("$subn: # Checked $count strains in %d families: \n%s\n\n", scalar(@fams), join("\n", @fams));
    if ($DOWNLOAD_TAXON) {
        print STDERR "$subn: performed live download\n";
    } else {
        print STDERR "$subn: existing file used, no live download\n";
    }
    if ( !$TAXON_UPDATED ) {
        print STDERR "$subn: Message: No update is found between taxon from NCBI and stored in '$exe_dir/$taxonFileName'\n";
        print STDERR "$subn: Message: No change is made to the current file '$exe_dir/$taxonFileName'\n";
    } else {
        print STDERR "$subn: Found update in taxon from NCBI, full info will be saved to file='$exe_dir/$taxonFileName'\n";
        my $bak1 = "$taxonFileName.1";
        my $result = `diff $exe_dir/$taxonFileName $exe_dir/$bak1`;
        if (-e "$exe_dir/$taxonFileName" && $result) {
            for my $i (reverse 1 .. 8) { # back up upto 8 earlier copies
                $bak1 = sprintf("Annotate_taxon_records.txt.%d", $i);
                next if (!-e "$exe_dir/$bak1");
                my $bak2 = sprintf("Annotate_taxon_records.txt.%d", $i+1);
                $result = `mv $exe_dir/$bak1 $exe_dir/$bak2`;
                print STDERR "$subn: \$result='$result'" if ($result);
            }
            $result = `mv $exe_dir/$taxonFileName $exe_dir/$bak1`;
            print STDERR "$subn: backing up $taxonFileName to $bak1 \$result='$result'\n" if ($result);
        }
        $result = `mv $exe_dir/$taxonTempName $exe_dir/$taxonFileName`;
        print STDERR "$subn: \$result='$result'\n";
        print STDERR "$subn: Moved $taxonTempName to $taxonFileName\n";
    }

    return $count;
} # sub checkAllTaxon
=cut

=head2
 sub checkAllRefseq look over each RefSeq (NC_??????) in $REFSEQS and see if there is any update in the genbank file
=cut

=head2
sub checkAllRefseq {
    my ( $exe_dir) = @_;

    my $debug = 0 && $debug_all;
    my $subn = 'checkAllRefseq';

    # Load the $REFSEQS if not already
    if ( !$REFSEQS->{refseq_loaded} ) {
        my $nstrain = Annotate_Def::initRefseq();
        my @fam = sort keys %{$REFSEQS->{refs}};
        printf STDERR ("$subn: loaded RefSeqs for $nstrain strains in %d families: \n%s\n\n", $#fam+1, join("\n", @fam));
    }
    $debug && print STDERR "$subn: \$REFSEQS=\n".Dumper($REFSEQS)."end of \$REFSEQS\n\n";

    my $count = 0;
    $exe_dir = './' if (!$exe_dir);
    $debug && print STDERR "$subn: \$exe_dir='$exe_dir'\n";

    my $refs = {};
    # Collect all the refseqs
    for my $fam (keys %{$REFSEQS->{refs}}){
        $debug && print STDERR "$subn: \$fam=$fam\n";
        my $list = $REFSEQS->{refs}->{$fam};
        for my $taxid (keys %$list) {
            $debug && print STDERR "$subn: \$taxid=$taxid refseq=$list->{$taxid}\n";
            $refs->{$list->{$taxid}}++;
        }
    }
    $debug && print STDERR "$subn: \$refs=\n".Dumper($refs)."end of \$refs\n\n";

    # Change to 1 in order to actually download from genbank
    my $DOWNLOAD_REFSEQ = 0;
    # Search for NC_* for each family from genbank
    for my $fam (sort keys %{$REFSEQS->{refs}}){
        my $web = '';
        my $key = '';
        my $count = 0;
        my $base = 'http://eutils.ncbi.nlm.nih.gov/entrez/eutils/';
        my $outFileName = '';
#        $outFileName = "RefSeq_${fam}.gbk";
#        print STDERR "\n$subn: \$fam=$fam \$outFileName=$outFileName\n";

        $outFileName = Annotate_Def::downloadRefseq( $fam, $exe_dir, $DOWNLOAD_REFSEQ);
        if (!$outFileName) {
            print STDERR "$subn: \$fam=$fam\n";
            print STDERR "$subn: ERROR: \$fam=$fam \$outFileName='$outFileName' is empty.\n";
            print STDERR "$subn: \$fam=$fam Please turn on the download options in sub Annotate_Def::downloadRefseq\n\n";
            next;
        }

        # Check each genbank file, if it has mat_peptide, compare with stored file, alert if different
        my $seqio = Bio::SeqIO->new( -file => "$exe_dir/$outFileName");
        my $nseq = 0;
        my $result = {}; # one line for each RefSeq
        my $status = {
                       total => 0,
                       no_file_on_disk => 0,
                       has_matpeptide => 0,
                       change_in_matpeptide => 0,
                     };
        while (my $new = $seqio->next_seq()) {
            $nseq++;
            my $acc = $new->accession_number;
            my $taxid = $new->species->ncbi_taxid;
            my $speciesid;
            $speciesid = Annotate_Def::getTaxonInfo( $taxid);
            $debug && print STDERR "$subn: \$nseq=$nseq \$acc=$acc \$taxid=$taxid \$speciesid='@$speciesid'\n";
            my $species = '';
            $species = ($#{$speciesid}>-1) ? $speciesid->[4] : '';
            $speciesid = ($#{$speciesid}>-1) ? $speciesid->[1] : '-1';
            $debug && print STDERR "$subn: \$nseq=$nseq \$acc=$acc\n";
            my $errcode = check1Refseq($new, $exe_dir);
            $debug && print STDERR "$subn: \$nseq=$nseq \$acc=$acc \$errcode=\n".Dumper($errcode)."\n";

#            print STDERR "$subn: \$nseq=$nseq \$acc=$acc\n";
            # print out those new mat_peptides
            if ($errcode->{change_in_matpeptide}>0) {
              for my $newf (@{$errcode->{new_mat_peptides}}) {
                print STDERR "$subn: \$nseq=$nseq \$acc=$acc newf='$newf'\n";
              }
            }

            $status->{total}++;
            $status->{no_file_on_disk}++ if ($errcode->{no_file_on_disk});
            $status->{has_matpeptide}++ if ($errcode->{has_matpeptide});
            $status->{change_in_matpeptide}++ if ($errcode->{change_in_matpeptide});
            my $msg1 = "#$nseq\t";
            $msg1 .= "$fam\t";
            $msg1 .= "$acc\t";
            $msg1 .= ($errcode->{has_matpeptide}) ? "$errcode->{has_matpeptide}\t" : "-\t";
            $msg1 .= ($errcode->{no_file_on_disk}) ? "$errcode->{no_file_on_disk}\t" : "-\t";
            $msg1 .= ($errcode->{change_in_matpeptide}) ? "$errcode->{change_in_matpeptide}\t" : "-\t";
            $msg1 .= "$taxid\t";
            $msg1 .= "$speciesid\t";
            $msg1 .= "'$species'\n";

            $result->{$speciesid}->{$acc} = $msg1;
        }

        # collect all summaries
        my $msg = "#\tfamily  \taccession\thas_mat\tno_file\tchanged\ttaxid\tspeciesid\tspecies\n";
        for my $speciesid (sort {$a<=>$b} keys %{$result}) {
          for my $acc (sort keys %{$result->{$speciesid}}) {
            $msg .= $result->{$speciesid}->{$acc};
          }
        }
        $msg .= "$subn: ERROR: total RefSeq $status->{total} doesn't match \$count=$count" if ($count && $count!=$status->{total});

        $msg .= "family=$fam\ttotal=$status->{total}\t";
        $msg .= "$status->{has_matpeptide}\t";
        $msg .= "$status->{no_file_on_disk}\t";
        $msg .= "$status->{change_in_matpeptide}\n";

        printf STDERR ("$msg\n");
        printf STDERR ("$subn: family=$fam total RefSeqs checked:  $status->{total}\n");
        printf STDERR ("$subn: family=$fam Refseqs w/ mat_peptide: $status->{has_matpeptide}/$status->{total}\n");
        printf STDERR ("$subn: family=$fam Refseqs w/o gbk file:   $status->{no_file_on_disk}/$status->{total}\n");
        printf STDERR ("$subn: family=$fam Refseqs w/ changes:     $status->{change_in_matpeptide}/$status->{total}\n\n");

    }

    return $count;
} # sub checkAllRefseq
=cut

=head2
 sub check1Refseq Bio::Seq object based on genbank file, checks if it has same feature annotation
 as anyfile in refseq directory. returs following values:
0: no difference, or existing gbk file
1. found difference in term of mat_peptide annotation, alert the user
=cut

=head2
sub check1Refseq {
    my ( $new, $exe_dir) = @_;

    my $debug = 0 && $debug_all;
    my $subn = 'check1Refseq';

    $exe_dir = './' if (!$exe_dir);
    my $errcode = {
                    no_file_on_disk => 0,
                    has_matpeptide => 0,
                    change_in_matpeptide => 0,
                    change_in_file => 0,
                    mat_peptides =>[],
                    new_mat_peptides =>[],
                  };

    my $acc = $new->accession_number;
    my $newfeats = [ $new->get_SeqFeatures ];
    $debug && print STDERR "$subn: $acc has $#{$newfeats} features in total\n";

    # Gather all mat_peptides
    my $has_matpeptide = 0;
    for (my $newct = 0; $newct<=$#{$newfeats}; $newct++) {
        my $newfeat = $newfeats->[$newct];
#        $debug && print STDERR "$subn: \$newct=$newct \$newfeat=".$newfeat->primary_tag.':'.$newfeat->location->to_FTstring."\n";
        next if ($newfeat->primary_tag ne 'CDS'); # Looking for the first CDS
        my $cds = $newfeat;
        my $cdsid = '';
        for my $tag ('db_xref', 'protein_id') {
            if ($cds->has_tag($tag)) {
               my @id = $cds->get_tag_values($tag);
               $cdsid = 'CDS='. $id[0];
               last;
            }
        }
        while ($newfeats->[$newct+1] && ($newfeats->[$newct+1]->primary_tag eq 'mat_peptide'
                 || $newfeats->[$newct+1]->primary_tag eq 'sig_peptide'
                 || $newfeats->[$newct+1]->primary_tag eq 'misc_feature') ) {
            $newfeat = $newfeats->[++$newct];
            next if ($newfeat->primary_tag ne 'mat_peptide');
            $errcode->{has_matpeptide} += 1;
            my $str = $acc."|$cdsid|Loc=".$newfeat->location->to_FTstring;
            push @{$errcode->{mat_peptides}}, $str;
        }
    }

    ($errcode->{has_matpeptide}==0) && return $errcode;

    my $old = undef;
    my $oldfeats = [ ];
    # See if genbank file exists on disk
    if (!-e "$exe_dir/refseq/$acc.gb") {
        $errcode->{no_file_on_disk} = 1;
#        return $errcode;
    } else {
        # Check between new and old RefSeq
        $old = Bio::SeqIO->new( -file => "$exe_dir/refseq/$acc.gb")->next_seq();
        $oldfeats = [ $old->get_SeqFeatures ];
    }
    my $msg = '';
    for (my $newct = 0; $newct<=$#{$errcode->{mat_peptides}}; $newct++) {
        my $newfeat = $errcode->{mat_peptides}->[$newct];
        my $newfeat1 = $newfeat;
        $newfeat1 =~ s/[<>]//g;
        $debug && print STDERR "$subn: \$newct=$newct \$newfeat='$newfeat'\n";

        my $seen = 0;
        for (my $oldct = 0; $oldct<=$#{$oldfeats}; $oldct++) {
            my $oldfeat = $oldfeats->[$oldct];
#            $debug && print STDERR "$subn: \$oldct=$oldct \$oldfeat=".$oldfeat->primary_tag.':'.$oldfeat->location->to_FTstring."\n";
            next if ($oldfeat->primary_tag ne 'CDS'); # Looking for the first CDS
            my $cds = $oldfeat;
            my $cdsid = '';
            for my $tag ('db_xref', 'protein_id') {
                if ($cds->has_tag($tag)) {
                    my @id = $cds->get_tag_values($tag);
                    $cdsid = 'CDS='. $id[0];
                    last;
                }
            }
            while ($oldfeats->[$oldct+1] && ($oldfeats->[$oldct+1]->primary_tag eq 'mat_peptide'
                    || $oldfeats->[$oldct+1]->primary_tag eq 'sig_peptide'
                    || $oldfeats->[$oldct+1]->primary_tag eq 'misc_feature')) {
                $oldfeat = $oldfeats->[++$oldct];
                my $oldfeat = $acc."|$cdsid|Loc=".$oldfeat->location->to_FTstring;
                my $oldfeat1 = $oldfeat;
                $oldfeat1 =~ s/[<>]//g;
                if ($newfeat eq $oldfeat) {
                  $debug && print STDERR "$subn: \$newct=$newct \$newfeat='$newfeat' \$oldfeat='$oldfeat'\n";
                  $seen = 1;
                  last;
                } elsif ($newfeat1 eq $oldfeat1) {
                  $debug && print STDERR "$subn: \$newct=$newct \$newfeat='$newfeat' \$oldfeat='$oldfeat' are similar\n";
                }
            }
        }
        if (!$seen) {
            $errcode->{change_in_matpeptide} += 1;
            push @{$errcode->{new_mat_peptides}}, $newfeat;
        }
    }
    
    return $errcode;
} # sub check1Refseq
=cut

=head2
 sub getTaxonInfo takes a taxid, returns the array of [taxid, speciesid, genusid, familyid, species, genus, family]
=cut

sub getTaxonInfo {
    my ( $taxid, $exe_dir) = @_;

    my $debug = 0 && $debug_all;
    my $subn = 'getTaxonInfo';

    $exe_dir = './' if (!$exe_dir);
    if ( !$TAXON->{taxon_loaded} ) {
        Annotate_Def::loadTaxonTable( $exe_dir);
        my $ctLines = scalar(keys %{$TAXON->{taxon}});
        $debug && print STDERR "$subn: loaded $ctLines lines from file '$exe_dir/$TAXON->{taxon_fn}'\n";
    }

    $debug && print STDERR "$subn: \$taxid=$taxid\n";
    my @taxInfo = ();
    if (exists($TAXON->{taxon}->{$taxid})) {
        @taxInfo = @{$TAXON->{taxon}->{$taxid}};
        $debug && print STDERR "$subn: \$taxid=$taxid \@taxInfo='@taxInfo'\n";
    } else {
        # See if $taxid is species
        for my $t (keys %{$TAXON->{taxon}}) {
            $debug && print STDERR "$subn: \$taxid=$taxid \$t='$t' species=$TAXON->{taxon}->{$t}->[1]\n";
            next if ($taxid != $TAXON->{taxon}->{$t}->[1]);
            @taxInfo = @{$TAXON->{taxon}->{$t}};
            last;
        }
    }
    return [ @taxInfo ];
} # sub getTaxonInfo


=head2
sub loadTaxonTable loads the taxonomy info from a text file into a global variable
=cut

sub loadTaxonTable {
    my ($exe_dir) = @_;

    my $debug = 0 && $debug_all;
    my $subn = 'loadTaxonTable';

    my $ctLines = 0;
    # Load the taxon information from an ASCII file
    $debug && print STDERR "$subn: \$TAXON=$TAXON->{taxon_loaded} file=$exe_dir/$TAXON->{taxon_fn}\n";
    if ( !-f "$exe_dir/$TAXON->{taxon_fn}" ) {
        croak "$subn: Couldn't find file=$exe_dir/$TAXON->{taxon_fn}\n";
        return 0;
    } else {
        $debug && print STDERR "$subn: found taxon_fn=$exe_dir/$TAXON->{taxon_fn}\n";

        open my $tfile, '<', "$exe_dir/$TAXON->{taxon_fn}"
               or croak("$subn: found '$exe_dir/$TAXON->{taxon_fn}', but couldn't open: $OS_ERROR");
        while (<$tfile>) {
            s/(^[|]*\s+|\s+$)//x;
            my $words = [ split(/\s*\|\s*/x) ];
#            $debug && print STDERR "$subn: \$words='@$words'\n";
            next if (!$words->[0] || $words->[0] !~ /^\d+$/x);
            next if (!$words->[1] || $words->[1] !~ /^\d+$/x); # ignore species=-1

            my $strainid = $words->[0];
            $ctLines++;
            $TAXON->{'taxon'}->{$strainid} = $words;

        }
        close $tfile or croak "$subn: Couldn't close $exe_dir/$TAXON->{taxon_fn}: $OS_ERROR";
        my @keys = keys %{$TAXON->{'taxon'}};
        $TAXON->{taxon_loaded} = 1 if ($#keys>1);

        $debug && print STDERR "$subn: finished reading $ctLines lines from file: '$exe_dir/$TAXON->{taxon_fn}'.\n";
        $debug && print STDERR "$subn: \$TAXON=\n".Dumper($TAXON)."End of \$TAXON\n\n";
    }
    $debug && print STDERR "$subn: \$TAXON->{taxon_loaded}=$TAXON->{taxon_loaded}\n";
    print STDERR "$subn: loaded $ctLines lines from file '$exe_dir/$TAXON->{taxon_fn}'\n";

    return $TAXON;
} # sub loadTaxonTable


=head2
sub printRefseqList prints the RefSeq in following format
Family=Flaviviridae
Family=Flaviviridae     species=12637   strain=11053    refseq=NC_001477
Family=Flaviviridae     species=12637   strain=11060    refseq=NC_001474
Family=Flaviviridae     species=11072   strain=         refseq=NC_001437
=cut

sub printRefseqList {
    my ($exe_dir) = @_;

    my $debug = 0 && $debug_all;
    my $subn = 'printRefseqList';

    $debug && print STDERR "$subn: Got here inside sub printRefseqList\n";
    $debug && print STDERR "$subn: \$REFSEQS=\n".Dumper($REFSEQS)."End of \$REFSEQS\n\n";
    $debug && print STDERR "$subn: \$TAXON=\n".Dumper($TAXON)."End of \$TAXON\n\n";
#    print a table of family, species, strain, RefSeq
    my $nfamily = 0;
    my $nspecies = 0;
    my $nstrain = 0;
    my $t = "";
    for my $fam (sort keys %{$REFSEQS->{refs}}) {
        my $t1 = "Family=$fam\n";
        $debug && print STDERR "$subn: \$fam=$fam\n";
        my $strains = $REFSEQS->{refs}->{$fam};
        if (!%$strains) {
            $t1 .= "Family=$fam \tNo species found for this family\n";
            $debug && print STDERR "$subn: \$t1='$t1'\n";
#            next;
        }
        $nfamily++;
        my $istrain = 0;
        for my $id (sort keys %$strains) {
            $istrain++;
            $t1 .= "$fam\t";
            $t1 .= "#$istrain\t";
            # Find the speciesid
            my $found = 0;
            my $speciesid = '';
            if (exists($TAXON->{'taxon'}->{$id})) {
                $found = 1;
                $speciesid = $TAXON->{'taxon'}->{$id}->[1];
                my $words = $TAXON->{'taxon'}->{$id};
                $t1 .= "species=$words->[1] \t";
                if ($words->[0] eq $words->[1]) {
                    $t1 .= "strain=----- \t";
                    $nspecies++;
                } else {
                    $t1 .= "strain=$words->[0] \t";
                }
                $nstrain++;
                $t1 .= "ref=$strains->{$id} \t\"$words->[4]\"\n";
            } else {
                for my $s (keys %{$TAXON->{'taxon'}}) {
                    next if ($TAXON->{'taxon'}->{$s}->[1] ne $id);
                    $found = 1;
                    $t1 .= "species=$id \tstrain=----- \t";
                    $t1 .= "ref=$strains->{$id} \t\"$TAXON->{'taxon'}->{$s}->[4]\"\n";
                    $nspecies++;
                    $nstrain++;
                    last;
                }
            }
            if (!$found) {
                print STDERR "$subn: ERROR: Species/strain=$id doesn't exist in \$TAXON, abort\n";
                exit(1);
            }
            $debug && print STDERR "$subn: Families: $nfamily, species: $nspecies, strains: $nstrain\n";
        }
        $t .= ($t) ? "\n".$t1 : $t1;
    }
    my $t1 = "Viral families and species covered in Ver$VERSION are:\n";
    $t1 .= "Families: $nfamily, species: $nspecies, strains: $nstrain\n";
    return $t1 . "\n". $t;
} # sub printRefseqList


=head2 get_refpolyprots
Takes a sequence object,
 return the [CDS_id, CDS, mat_peptide 1, mat_peptide 2, ...] in the refseq. Empty if there is any problem
=cut

sub get_refpolyprots {
    my ( $refseqs, $inseq, $exe_dir) = @_;

    my $debug = 0 && $debug_all;
    my $subn = 'get_refpolyprots';

    my $refpolyprots = [];
    my $refseq;

    my @ann = $inseq->annotation->get_Annotations('date_changed');
    $debug && print STDERR "$subn: accession=".$inseq->accession_number."\tdate=$ann[0]->{'value'}\n";

    if (!$refseq || !$refseq->isa('Bio::Seq::RichSeq') || !check_refseq($refseq, $inseq)) {

        $refseq = Annotate_Def::get_refseq( $inseq, $exe_dir);
        my $acc = $inseq->accession;
        my $taxid  = $inseq->species->ncbi_taxid;
        if (!$refseq) {
#            my $species = Annotate_Def::getSpecies( $taxid);
            my $species = Annotate_Def::getTaxonInfo( $taxid);
            $species = $species->[4];
            print STDERR "$subn: ERROR genome=$acc w/ taxid=$taxid($species) is not covered in V$VERSION.";
            print STDERR " Please contact script author for any update.\n";
            return undef;
        } else {
#            print "$subn: \$refseq=$refseq\n\n";
            $debug && print STDOUT "$subn: Refseq changed to \$refseq=".$acc.' w/ taxid='.$taxid."\n";
            $debug && print STDERR "$subn: Refseq changed to \$refseq=".$acc.' w/ taxid='.$taxid."\n";
        }
    }
    if (exists($refseqs->{$refseq->accession_number})) {
        return $refseqs->{$refseq->accession_number};
    }

    # For a new refseq, go through all features in refseq, look for polyproptein CDS
    my $reffeats = [ $refseq->get_SeqFeatures ];
    $debug && print STDERR "$subn: ".$refseq->accession_number." has $#{$reffeats} features\n";

    for (my $refct = 0; $refct<=$#{$reffeats}; $refct++) {
        my $reffeat = $reffeats->[$refct];
        $debug && print STDERR "$subn: \$refct=$refct \$reffeat=".$reffeat->primary_tag.':'.$reffeat->location->to_FTstring."\n";
        next if ($reffeat->primary_tag ne 'CDS'); # Looking for the first CDS

        my $acc = $reffeat->seq->accession_number;
        my $is_poly = 0;
        $is_poly = Annotate_Def::is_polyprotein( $reffeat, ['product', 'note'], $reffeats->[$refct+1]);
        if (!$is_poly) {
            $debug && print STDERR "$subn: \$refct=$refct \$is_poly=$is_poly is not a polyprotein, skipping\n";
            next; # only run for CDS labeled as polyprotein in refseq
        } else {
#            $debug && print STDERR "$subn: \$refct=$refct \$reffeat=".$reffeat->primary_tag."\n";
        }

        my $refcds = $reffeat;
        $debug && print STDERR "$subn: Found polyprotein at \$refct=$refct\n";

        # Found a polyprotein CDS, now collect all mat_peptide and sig_peptide following it
        my $refmatps = [];
        my %allowed_tags = (
                             'mat_peptide' => 1,
                             'sig_peptide' => 1,
                             'misc_feature' => 1,
                             'proprotein' => 1, # This is needed for refseq NC_003988 of Simian enterovirus A (310907)
                           );
        while ($refct<$#{$reffeats} && $allowed_tags{$reffeats->[$refct+1]->primary_tag}) {

            my $reffeat = $reffeats->[$refct+1];
            # Special cases to correct refseq location of preM of NC_001809
            if ($acc eq "NC_001809" && $reffeat->location->to_FTstring eq "466..972") {
                $reffeat->location->end(969); # The annotation in gbk has 1 extra codon at the end
                print STDERR "$subn: \$reffeat=\n".Dumper($reffeat)."End of $reffeat\n";
            }
            # Special cases to correct refseq location of preM of NC_003635
            if ($acc eq "NC_003635" && $reffeat->location->to_FTstring eq "<440..>700") {
                $reffeat->location->end(925); # The annotation stops at the pre part
                print STDERR "$subn: \$reffeat=\n".Dumper($reffeat)."End of $reffeat\n";
            }
            
            $debug && print STDERR "$subn: \$refct=$refct+1 \$reffeat=".$reffeats->[$refct+1]->primary_tag.' '.$reffeats->[$refct+1]->location->to_FTstring."\n";
            if ($reffeats->[$refct+1]->primary_tag eq 'sig_peptide') {
                ++$refct;
                next;
            };
            if ($reffeats->[$refct+1]->primary_tag eq 'misc_feature') {
                ++$refct;
                next;
            };
            if ($reffeats->[$refct+1]->primary_tag eq 'proprotein') {
                ++$refct;
                next;
            };
            
            push @$refmatps, $reffeats->[++$refct];
        }

        my $cds_id = '';
        my @allowed_tags = ( 'db_xref', 'protein_id');
        foreach my $tag (@allowed_tags) {
            next if (!$refcds->has_tag($tag));
            my @values = $refcds->get_tag_values($tag);
            foreach my $v (@values) {
                if ($v =~ /^(GI:.+)$/i) {
                    $cds_id = $1;
                } elsif ($v =~ /^(NP_\d+)[.]\d+$/i) {
                    $cds_id = $1;
                }
                $debug && $cds_id && print STDERR "$subn: \$cds_id=$cds_id\n";
                last if ($cds_id);
            }
            last if ($cds_id);
        }
        $cds_id = 'unknown' if (!$cds_id);


        # check if this is a new CDS, e.g. NC_001547.1 has a 2nd CDS that's part of 1st one
        # If so, choose the longer CDS, and combine the list of mat_peptides
        my $seen=0;
if (0) {
#        foreach my $j (keys %$polyprots) {
        foreach my $k (0 .. $#{$refpolyprots}) {
            my $old_cds = $refpolyprots->[$k]->[1];
            $debug && print STDERR "$subn: \$k=$k \$old_cds=\n".Dumper($old_cds)."end of \$old_cds\n\n";
            $debug && print STDERR "$subn: \$k=$k \$refcds=\n".Dumper($refcds)."end of \$refcds\n\n";
            my $s1 = Annotate_Def::get_new_translation( $old_cds, $old_cds);
            my $s2 = Annotate_Def::get_new_translation( $refcds, $refcds);
            $debug && print STDERR "$subn: \$s1=$s1\n";
            $debug && print STDERR "$subn: \$s2=$s2\n";
            if ($s1 =~ /$s2/i || $s2 =~ /$s1/i) { # $s1 is related to $s2
              $seen=1;
              printf STDERR "$subn: old CDS=%d..%d is related to CDS=%d..%d, skip the shorter one\n", $old_cds->location->start, $old_cds->location->end, $refcds->location->start, $refcds->location->end;
              $old_cds = $refcds if ($s2 =~ /$s1/i);
              # See if each mat_peptide is a duplicate, leave out duplicates. E.g. NC_003899
              foreach my $refmatp (@$refmatps) {
                $debug && print STDERR "\n$subn: \$k=$k \$refmatp=".$refmatp->location->to_FTstring."\n";
                my $seen1 = 0;
                for my $n (2 .. $#{$refpolyprots->[$k]}) {
                  $debug && print STDERR "$subn: \$k=$k \$n=$n \$refmatp=".$refpolyprots->[$k]->[$n]->location->to_FTstring."\n";
                  my $refpolyprot = $refpolyprots->[$k]->[$n]->location;
                  if ($refmatp->location->start == $refpolyprot->start && $refmatp->location->end == $refpolyprot->end) {
                     $seen1 = 1;
                  }
                }
                $debug && print STDERR "$subn: \$k=$k \$refmatp=".$refmatp->location->to_FTstring." \$seen1=$seen1\n";
                push @{$refpolyprots->[$k]}, $refmatp if (!$seen1);
              }
#              last;
            } else {
              # if $s1 is not related to $s2
              $debug && print STDERR "$subn: \$s1 is not related to \$s2\n";
            }
        }
#        }
}
        if (!$seen && $#{$refmatps}>=0) {
            $debug && print STDERR "$subn: \$refct=$refct $#{$refmatps}+1 mat_peptides in ".$refseq->accession_number."|$cds_id\n";
            push @$refpolyprots, [
                               $cds_id,
                               $refcds,
                               @$refmatps,
                                 ];
        }
        $debug && print STDERR "$subn: \$refct=$refct \$#refpolyprots=$#{$refpolyprots}+1\n";
    } # for (my $refct = 0; $refct<=$#{$reffeats}; $refct++) {
    $debug && print STDERR "$subn: \$refpolyprots=\n".Dumper($refpolyprots)."end of \$refpolyprots\n\n";

    return $refpolyprots;
} # sub get_refpolyprots


sub is_polyprotein {
    my ($feat, $allowed_tags, $feat2) = @_;
#        $is_poly = Annotate_Def::is_polyprotein( $reffeat, ['product', 'note']);

    my $debug = 0 && $debug_all;
    my $subn = 'is_polyprotein';

    my $is_poly = 0;
    foreach my $tag (@$allowed_tags) {
#        $debug && print STDERR "$subn: \$refct=$refct \$tag=$tag\n";
        next if (!$feat->has_tag($tag));
#        $debug && print STDERR "$subn: \$refct=$refct Found \$tag=$tag\n";
        my @tag = $feat->get_tag_values($tag);
        foreach my $t (@tag) {
#            $debug && print STDERR "$subn: \$refct=$refct \$t=$t\n";
            if ($t =~ /polyprotein/i) {
                $is_poly = 1;
                last;
            }
        }
        last if $is_poly;
    }

    if (!$is_poly && $feat2) { {
        $debug && print STDERR "$subn: \$is_poly=$is_poly \$feat2=$feat2\n";
        last if ($feat2->primary_tag ne 'mat_peptide');
        last if ($feat2->location->start < $feat->location->start);
        last if ($feat->location->end   < $feat2->location->end);
        $is_poly = 1;
    } }

    return $is_poly;
} # sub is_polyprotein


=head2 get_refseq
Takes either a file name, or a Bio::Seq object based on genbank file, returns refseq in Bio::Seq.
 For a file name, simply load the genbank file
 For a Bio::Seq object, look at its taxon. Load the refseq either from file system (or MySQL)
 Return the refseq in Bio::Seq object.
=cut

sub get_refseq {
    my ($inseq, $exe_dir) = @_;

    my $debug = 0 && $debug_all;
    my $subn = 'get_refseq';
    my $refseq = undef;
    return $refseq if (!$inseq);

    $debug && print STDERR "$subn: \$inseq=$inseq is a ".ref($inseq)."\n";
    my $speciesid = '-1';
    my $species = 'unknown';
    my $taxid = '';
    my $refseq_fn;

    if ($inseq->isa('Bio::Seq::RichSeq')) {
            my $refacc = undef;
            $taxid = $inseq->species->ncbi_taxid;
            ($refacc, $speciesid, $species) = Annotate_Def::get_refseq_acc( $inseq, $exe_dir);
            return undef if (!$refacc);
            if (-e "${exe_dir}refseq/${refacc}_matpeptide.gb") {
                $refseq_fn = $refacc.'_matpeptide.gb';
            } elsif (-e "${exe_dir}refseq/${refacc}_msaa.gb") {
                $refseq_fn = $refacc.'_msaa.gb';
            } elsif (-e "${exe_dir}refseq/${refacc}.gb") {
                $refseq_fn = $refacc.'.gb';
            } elsif (-e "${exe_dir}refseq/${refacc}.gbk") {
                $refseq_fn = $refacc.'.gbk';
            } elsif (-e "${exe_dir}refseq/${refacc}.genbank") {
                $refseq_fn = $refacc.'.genbank';
            } else {
                $refseq_fn = $refacc.'.gb';
            }
            $debug && print STDERR "$subn: REFSEQ is ${exe_dir}refseq/$refseq_fn\n";

    } elsif ($inseq) {
        $refseq_fn = $inseq;
    }

    # Get refseq object if $refseq_fn is given, search in ./refseq, ./, and $dir_path/$ticket/
    # First look for refseq in file system
    $exe_dir = './' if (!$exe_dir);
    if ($refseq_fn) {
        $debug && print STDERR "$subn: REFSEQ is $refseq_fn\n";
        if (-e "$exe_dir/refseq/$refseq_fn") {
            $refseq = Bio::SeqIO->new( -file => "$exe_dir/refseq/$refseq_fn")->next_seq();
            print STDERR "$subn: Found REFSEQ at '$exe_dir/refseq/$refseq_fn' for taxid=$taxid, species=$speciesid ($species)\n";
        } elsif (-e "$exe_dir/$refseq_fn") {
            $refseq = Bio::SeqIO->new( -file => "$exe_dir/$refseq_fn")->next_seq();
            print STDERR "$subn: Found REFSEQ at '$exe_dir/$refseq_fn' for taxid=$taxid, species=$speciesid ($species)\n";
#        } elsif ($debug && -e "./$refseq_fn") {
#            $refseq = Bio::SeqIO->new( -file => "./$refseq_fn")->next_seq();
#            print STDERR "$subn: Found REFSEQ at './$refseq_fn'\n";
        } else {
            print STDERR "$subn: ERROR: Can't find required refseq: $refseq_fn in '$exe_dir/refseq' for taxid=$taxid, species=$speciesid ($species).\n";
        }
    }

    return $refseq;
} # sub get_refseq


=head2 get1RefseqAcc
Takes a taxid, determine the RefSeq from $REFSEQS
 Return the refseq accession
=cut

sub get1RefseqAcc {
    my ($taxid) = @_;

    my $debug = 0 && $debug_all;
    my $subn = 'get1RefseqAcc';

    my $refseq_acc = '';
    $debug && print STDERR "$subn: \$taxid=$taxid\n";
    return $refseq_acc if (!$taxid);

    # Going through $REFSEQS and find the RefSeq for a given taxid
    my $speciesid = '-1';
    my $species = 'unknown';
    for my $fam (keys %{$REFSEQS->{refs}}) {
        my $refseq_list1 = $REFSEQS->{refs}->{$fam};
        $debug && print STDERR "$subn: \$fam=$fam \$refseq_list1=\n".Dumper($refseq_list1)."\n";
        next if (!($refseq_list1));
        $debug && print STDERR "$subn: \$refseq_list1=".keys(%$refseq_list1)."\n";
        if (exists($refseq_list1->{$taxid})) {
            $refseq_acc = $refseq_list1->{$taxid};
            my $taxinfo = Annotate_Def::getTaxonInfo( $taxid);
            $speciesid = $taxinfo->[1] if (scalar($taxinfo));
            $species = $taxinfo->[4] if (scalar($taxinfo));
            $debug && print STDERR "$subn: Got from \$REFSEQS \$refseq_acc=$refseq_acc for \$taxid=$taxid species=$speciesid ($species)\n";
        }
    }
    # If $REFSEQS doesn't contain $taxid, try find the species and take the refseq there
    if (!$refseq_acc && $TAXON->{taxon_loaded}) {
        my $taxinfo = Annotate_Def::getTaxonInfo( $taxid);
        $speciesid = $taxinfo->[1] if (scalar($taxinfo));
        $species = $taxinfo->[4] if (scalar($taxinfo));
        if ($speciesid>0) {
            $speciesid = $TAXON->{'taxon'}->{$taxid}->[1];
            $species = $TAXON->{'taxon'}->{$taxid}->[4];
    for my $fam (keys %{$REFSEQS->{refs}}) {
        my $refseq_list1 = $REFSEQS->{refs}->{$fam};
            if (exists($refseq_list1->{$speciesid})) {
                $refseq_acc = $refseq_list1->{$speciesid};
                $debug && print STDERR "$subn: Got from \$REFSEQS \$refseq_acc=$refseq_acc for \$taxid=$taxid species=$speciesid ($species)\n";
            }
        }
    }
=head2
        $refseq_acc = Annotate_Def::get1RefseqAcc($speciesid);
            next if (!exists($TAXON->{'taxon'}->{$speciesid}->{$taxid}));
            # Could use 0 to signal the strains not covered
#            next if ($TAXON->{'taxon'}->{$speciesid}->{$taxid}!=1);
            my $refacc = '';
            $refacc = $TAXON->{'taxon'}->{$speciesid}->{$speciesid} if (exists($TAXON->{'taxon'}->{$speciesid}->{$speciesid}));
            $debug && print "$subn: \$taxid=$taxid \$speciesid=$speciesid \$refacc=$refacc\n";
            $refseq_acc = $refacc if ($refacc && $refacc =~ /^NC_0/i);
            print STDERR "$subn: \$taxid=$taxid \$speciesid=$speciesid \$refacc=$refacc \$refseq_acc=$refseq_acc\n";
            last;
#        }
=cut
        $debug && print STDERR "$subn: Determined from \$REFSEQS & \$TAXON \$refseq_acc=$refseq_acc for \$taxid=$taxid species=$speciesid ($species)\n";

    }

    return ($refseq_acc, $speciesid, $species);
} # sub get1RefseqAcc


=head2 get_refseq_acc

Takes either a file name, or a Bio::Seq object based on genbank file, returns refseq in Bio::Seq.
 For a file name, simply load the genbank file
 For a Bio::Seq object, look at its taxon. Load the refseq either from file system (or MySQL)
 Return the refseq in Bio::Seq object.
=cut

sub get_refseq_acc {
    my ($inseq, $exe_dir) = @_;

    my $debug = 0 && $debug_all;
    my $subn = 'get_refseq_acc';

    # Initialize the $REFSEQS if needed
    if ( !$REFSEQS->{refseq_loaded} ) {
        my $nstrain = Annotate_Def::initRefseq();
        my @fam = sort keys %{$REFSEQS->{refs}};
        $debug && print STDERR "$subn: loaded $nstrain strain for $#fam families: @fam\n";
    }
    $debug && print STDERR "$subn: \$REFSEQS=\n".Dumper($REFSEQS)."end of \$REFSEQS\n\n";

    # Load taxon information from an ASCII file
    $debug && print STDERR "$subn: \$TAXON=$TAXON->{taxon_loaded} file=$exe_dir/$TAXON->{taxon_fn}\n";
    if ( !$TAXON->{taxon_loaded} ) {
        Annotate_Def::loadTaxonTable( $exe_dir);
    $debug && print STDERR "$subn: \$TAXON=\n".Dumper($TAXON)."end of \$TAXON\n\n";
        my $ctLines = scalar(keys %{$TAXON->{taxon}});
        $debug && print STDERR "$subn: loaded $ctLines lines from file '$exe_dir/$TAXON->{taxon_fn}'\n";
    }
    $debug && print STDERR "$subn: \$TAXON->{taxon_loaded}=$TAXON->{taxon_loaded}\n";
    $debug && print STDERR "$subn: \$TAXON=\n".Dumper($TAXON)."end of \$TAXON\n\n";

#    If no inseq, print a table of family, species, strain, RefSeq
    if (!$inseq) {
        $debug && print STDERR "$subn: Going into sub printRefseqList\n";
        return Annotate_Def::printRefseqList();
    }

    my $speciesid = '-1';
    my $species = 'unknown';
    my $refseq_acc = '';
    my $taxid = $inseq->species->ncbi_taxid;
    $debug && print STDERR "$subn: \$taxid=$taxid \$inseq=$inseq is a ".ref($inseq)."\n";


    # Determine the refseq by a 2-step process
    # First, see if the $REFSEQS contains $taxid, take the refseq if so
if (0) {
    for my $fam (keys %{$REFSEQS->{refs}}) {
        my $refseq_list1 = $REFSEQS->{refs}->{$fam};
        $debug && print STDERR "$subn: \$fam=$fam \$refseq_list1=\n".Dumper($refseq_list1)."\n";
        next if (!($refseq_list1));
        $debug && print STDERR "$subn: \$refseq_list1=".keys(%$refseq_list1)."\n";
        if (%$refseq_list1 && exists($refseq_list1->{$taxid})) {
            $refseq_acc = $refseq_list1->{$taxid};
            $debug && print STDERR "$subn: Determined from \$REFSEQS \$taxid=$taxid \$refseq_acc=$refseq_acc\n";
        }
    }
} else {
    ($refseq_acc, $speciesid, $species) = Annotate_Def::get1RefseqAcc($taxid);
}

    return ($refseq_acc, $speciesid, $species);
} # sub get_refseq_acc


=head2
sub initRefseq stores the approved RefSeq for its species/strain in a global variable, which is not accessible from outside.
=cut

sub initRefseq {
    my ($inseq, $exe_dir) = @_;

    my $debug = 0 && $debug_all;
    my $subn = 'initRefseq';

#    my $refseq_acc = undef;
#    return $refseq_acc if (!$inseq);

        # list of refseqs
        $REFSEQS->{refseq_loaded} = 1;
        $REFSEQS->{refs} = {
        # strain_id => 'accession_refseq',
           # Family=Flaviviridae
           'Flaviviridae' => {
           11079 => 'NC_000943', # Murray Valley encephalitis virus; Ver1.1.3
           11072 => 'NC_001437', # Japanese encephalitis virus
##           11073 => 'NC_001437', # Japanese encephalitis virus. Refseq has 13 mat_peptides, w/ fuzzy ranges
##           11076 => 'NC_001437', # Japanese encephalitis virus. Refseq has 13 mat_peptides, w/ fuzzy ranges
           11060 => 'NC_001474', # Dengue virus 2
           11069 => 'NC_001475', # Dengue virus 3
           11053 => 'NC_001477', # Dengue virus 1
           12637 => 'NC_001477', # Dengue virus
           11082 => 'NC_001563', # West Nile virus (lineage II strain 956)
           31658 => 'NC_001564', # Cell fusing agent virus; Ver1.1.3
           11084 => 'NC_001672', # Tick-borne encephalitis virus. Refseq has 13 mat_peptides, w/ fuzzy ranges
##           11092 => 'NC_001672', # Tick-borne encephalitis virus. Refseq has 13 mat_peptides, w/ fuzzy ranges
##           11094 => 'NC_001672', # Tick-borne encephalitis virus. Refseq has 13 mat_peptides, w/ fuzzy ranges
##           47300 => 'NC_001672', # Tick-borne encephalitis virus. Refseq has 13 mat_peptides, w/ fuzzy ranges
##          638787 => 'NC_001672', # Tick-borne encephalitis virus. Refseq has 13 mat_peptides, w/ fuzzy ranges
           11086 => 'NC_001809', # Louping ill virus; Ver1.1.3; Special treatment to correct end point of preM
           11089 => 'NC_002031', # Yellow fever virus (YFV). This refseq has 13 mat_peptides, w/ fuzzy ranges
           11070 => 'NC_002640', # Dengue virus 4
           64300 => 'NC_003635', # Modoc virus; Ver1.1.3
           64285 => 'NC_003675', # Rio Bravo virus; Ver1.1.3
           64280 => 'NC_003676', # Apoi virus; Ver1.1.3
           11083 => 'NC_003687', # Powassan virus; Ver1.1.3
##           11085 => 'NC_003690', # Langat virus
##          161675 => 'NC_003996', # Tamana bat virus
           11103 => 'NC_004102', # Hepatitis C virus genotype 1
#           11103 => 'AF211032', # Used to test complement location in reffeat, need to comment out in production
           31646 => 'NC_004102', # Hepatitis C virus genotype 1a
           31647 => 'NC_004102', # Hepatitis C virus genotype 1b
           31649 => 'NC_004102', # Hepatitis C virus genotype 2a
           31650 => 'NC_004102', # Hepatitis C virus genotype 2b
           31655 => 'NC_004102', # Hepatitis C virus genotype 6a
          356426 => 'NC_004102', # Hepatitis C virus genotype 3a
           64312 => 'NC_004119', # Montana myotis leukoencephalitis virus; Ver1.1.3
          172148 => 'NC_004355', # Alkhurma hemorrhagic fever virus; Ver1.1.3
##           64294 => 'NC_005039', # Yokose virus
##           12542 => 'NC_005062', # Omsk hemorrhagic fever virus
##          218849 => 'NC_005064', # Kamiti River virus
           64286 => 'NC_006551', # Usutu virus; Ver1.1.3
##           64287 => 'NC_006947', # Karshi virus
           11080 => 'NC_007580', # St. Louis encephalitis virus
##           11080 => 'NC_001563', # St. Louis encephalitis => West Nile virus (lineage II strain 956)
##          390844 => 'NC_008604', # Culex flavivirus
##           64283 => 'NC_008718', # Entebbe bat virus
##           44026 => 'NC_008719', # Sepik virus
##           64303 => 'NC_009026', # Bussuquara virus
##           59563 => 'NC_009028', # Ilheus virus
##           44024 => 'NC_009029', # Kokobera virus
#           40271 => 'NC_009823', # Hepatitis C virus genotype 2
#          356114 => 'NC_009824', # Hepatitis C virus genotype 3
#           33745 => 'NC_009825', # Hepatitis C virus genotype 4
#           33746 => 'NC_009826', # Hepatitis C virus genotype 5
#           42182 => 'NC_009827', # Hepatitis C virus genotype 6
#           11082 => 'NC_009942', # West Nile virus (lineage I strain NY99), missing 2k
           },

           # Family=Caliciviridae
           'Caliciviridae' => {
#mysql> select taxid,accession,speciesid,species from genome left join taxon on taxid=taxon.id where accession like "NC_%" and taxid in (select id from taxon where family="Caliciviridae") order by speciesid,taxid,accession;
#+---------+-----------+-----------+---------------------------------------+
#| taxid   | accession | speciesid | species                               |
#+---------+-----------+-----------+---------------------------------------+
  11976 => 'NC_001543', #   11976 | Rabbit hemorrhagic disease virus      | Good refseq w/ 8 mat_peptides, V1.1.6
  11978 => 'NC_001481', #   11978 | Feline calicivirus                    | Good refseq w/ 8+2 mat_peptides, V1.1.6
  11983 => 'NC_001959', #   11983 | Norwalk virus                         | Good refseq w/ 6 mat_peptides, V1.1.1
  33756 => 'NC_002615', #   33756 | European brown hare syndrome virus    | Good refseq w/ 7 mat_peptides, V1.1.6
  35612 => 'NC_002551', #   35612 | Vesicular exanthema of swine virus    | Good refseq w/ 7+2 mat_peptides, V1.1.6
  74724 => 'NC_004542', #   74724 | Canine calicivirus                    | Good refseq w/ 7+2 mat_peptides, V1.1.6
## 106333 => 'NC_000940', #   95342 | Sapporo virus                         | No mat_peptide in refseq
## 234601 => 'NC_010624', #   95342 | Sapporo virus                         | No mat_peptide in refseq
# 290314 => 'NC_006554', #   95342 | Sapporo virus                         | Good refseq w/ 1 mat_peptides
## 291175 => 'NC_006269', #   95342 | Sapporo virus                         | No mat_peptide in refseq
 146073 => 'NC_004541', #  146073 | Walrus calicivirus                    | Good refseq w/ 7+2 mat_peptides, V1.1.6
## 303317 => 'NC_008580', #  303317 | Rabbit vesivirus                      | No mat_peptide in refseq
 # NC_008311 seems good enough to be used as RefSeq for species 357231
 223997 => 'NC_008311', #  357231 | Murine norovirus                      | Good refseq w/ 6 mat_peptides, V1.1.6
 357231 => 'NC_008311', #  357231 | Murine norovirus                      | Good refseq w/ 6 mat_peptides, V1.1.6
## 436911 => 'NC_011050', #  436911 | Steller sea lion vesivirus            | No mat_peptide in refseq
## 576948 => 'NC_011704', #  576948 | Rabbit calicivirus Australia 1 MIC-07 | No mat_peptide in refseq
## 520973 => 'NC_012699', #  646294 | St-Valerien swine virus               | No mat_peptide in refseq
## 190239 => 'NC_004064', #  696856 | Newbury-1 virus                       | No mat_peptide in refseq
# 331642 => 'NC_007916', #  696856 | Newbury-1 virus                       | Good refseq w/ 2 mat_peptides
##1185359 => 'NC_017936', # 1185359 | Bat sapovirus TLC58/HK                | No mat_peptide in refseq
#+---------+-----------+-----------+---------------------------------------+
#19 rows in set (0.24 sec)
           },

           # Family=Togaviridae
           'Togaviridae' => {
  # Alphavirus
  11036 => 'NC_001449', # 11036 |     11036 | Venezuelan equine encephalitis virus; Ver1.1.3
  11027 => 'NC_001512', # 11027 |     11027 | O'nyong-nyong virus; Ver1.1.3
  11029 => 'NC_001544', # 11029 |     11029 | Ross River virus; Ver1.1.3
  11034 => 'NC_001547', # 11034 |     11034 | Sindbis virus; Ver1.1.3
  11020 => 'NC_001786', # 11020 |     11020 | Barmah Forest virus; Ver1.1.3
  11033 => 'NC_003215', # 11033 |     11033 | Semliki forest virus; Ver1.1.3
  59301 => 'NC_003417', # 59301 |     59301 | Mayaro virus; Ver1.1.3
  78540 => 'NC_003433', # 78540 |     84589 | Salmon pancreas disease virus; Ver1.1.3
  11021 => 'NC_003899', # 11021 |     11021 | Eastern equine encephalitis virus; Ver1.1.3
  44158 => 'NC_003900', # 44158 |     44158 | Aura virus; Ver1.1.3
  11039 => 'NC_003908', # 11039 |     11039 | Western equine encephalomyelitis virus; Ver1.1.3
  84589 => 'NC_003930', # 84589 |     84589 | Salmon pancreas disease virus; Ver1.1.3
#  37124 => 'NC_001512', # 37124 |     37124 | Chikungunya virus, refseq=O'nyong-nyong virus
  37124 => 'NC_004162', # 37124 |     37124 | Chikungunya virus; Ver1.1.3
#  59300 => 'NC_001544', # 59300 |     59300 | Getah virus, refseq=Ross River virus
  59300 => 'NC_006558', # 59300 |     59300 | Getah virus; Ver1.1.3
  11024 => 'NC_012561', # 11024 |     11024 | Highlands J virus; Ver1.1.3
  48544 => 'NC_013528', # 48544 |     48544 | Fort Morgan virus; Ver1.1.3
  # Rubivirus
  11041 => 'NC_001545', # 11041 |     11041 | Rubella virus; Ver1.1.3
            },

           # Family=Coronaviridae
           'Coronaviridae' => {
#+---------+-----------+-----------+---------------------------------------------------------+
#| taxid   | accession | speciesid | species                                                 |
#+---------+-----------+-----------+---------------------------------------------------------+
#  11137 => 'NC_002645', #     11137 | Human coronavirus 229E          | No mat_peptide in refseq
  28295 => 'NC_003436', #     28295 | Porcine epidemic diarrhea virus | Has 13 mat_peptides, but w/ gaps b/w them, species
#  46473 => 'NC_007447', #     74501 | Bovine torovirus                | No mat_peptide in refseq
# 277944 => 'NC_005831', #    277944 | Human coronavirus NL63          | No mat_peptide in refseq
 290028 => 'NC_006577', #    290028 | Human coronavirus HKU1          | Has 15 mat_peptides; Ver1.1.4, species
# 389230 => 'NC_008315', #    389230 | Bat coronavirus (BtCoV/133/2005)| No mat_peptide in refseq
# 393767 => 'NC_010437', #    393767 | Bat coronavirus 1A              | No mat_peptide in refseq
# 393768 => 'NC_010436', #    393768 | Bat coronavirus 1B              | No mat_peptide in refseq
# 405554 => 'NC_008516', #    405554 | White bream virus               | No mat_peptide in refseq
# 502105 => 'NC_012949', #    502105 | Bovine respiratory coronavirus bovine/US/OH-440-TC/1996 | No mat_peptide in refseq
# 502108 => 'NC_012948', #    502108 | Bovine respiratory coronavirus AH187 | No mat_peptide in refseq
# 627439 => 'NC_012950', #    627439 | Human enteric coronavirus strain 4408| No mat_peptide in refseq
#  11135 => 'NC_002306', #    693997 | Alphacoronavirus 1              | No mat_peptide in refseq
# 693998 => 'NC_009988', #    693998 | Rhinolophus bat coronavirus HKU2| No mat_peptide in refseq
# 693999 => 'NC_009657', #    693999 | Scotophilus bat coronavirus 512 | No mat_peptide in refseq
# 694001 => 'NC_010438', #    694001 | Miniopterus bat coronavirus HKU8| No mat_peptide in refseq
# The species Bovine coronavirus (11128) has been lumped into species Betacoronavirus 1 (694003), and is no longer a species
# Alignment of all polyproteins in 694003 showed no gaps around the cleavage sites
  11128 => 'NC_003045', #    694003 | Betacoronavirus 1               | has mat_peptide; Ver1.1.4
 694003 => 'NC_003045', #    694003 | Betacoronavirus 1               | Has 15 mat_peptides; Ver1.1.5, species
#  31631 => 'NC_005147', #    694003 | Betacoronavirus 1               | No mat_peptide in refseq
#  42005 => 'NC_007732', #    694003 | Betacoronavirus 1               | No mat_peptide in refseq
# 136187 => 'NC_010327', #    694003 | Betacoronavirus 1               | No mat_peptide in refseq
# Murine hepatitis virus (11138) that used to be species has been moved to under Murine coronavirus (694005)
  11142 => 'NC_001846', #    694005 | Murine coronavirus              | has mat_peptide; Ver1.1.4
 694005 => 'NC_001846', #    694005 | Murine coronavirus              | has mat_peptide; Ver1.1.4, species
  11144 => 'NC_006852', #    694005 | Murine coronavirus              | has mat_peptide; Ver1.1.4
# 502102 => 'NC_012936', #    694005 | Murine coronavirus              | No mat_peptide in refseq
# 694006 => 'NC_009021', #    694006 | Rousettus bat coronavirus HKU9  | No mat_peptide in refseq
# 694007 => 'NC_009019', #    694007 | Tylonycteris bat coronavirus HKU4    | No mat_peptide in refseq
# 694008 => 'NC_009020', #    694008 | Pipistrellus bat coronavirus HKU5    | No mat_peptide in refseq
 227859 => 'NC_004718', #    694009 | Severe acute respiratory syndrome-related coronavirus| Has 15 mat_peptides; Ver1.1.4, species
 694009 => 'NC_004718', #    694009 | Severe acute respiratory syndrome-related coronavirus| Has 15 mat_peptides; Ver1.1.4, species
  11120 => 'NC_001451', #    694014 | Avian coronavirus               | Has 14 mat_peptides; Ver1.1.4
  11152 => 'NC_010800', #    694014 | Avian coronavirus               | Has 15 mat_peptides; Ver1.1.4
 694014 => 'NC_010800', #    694014 | Avian coronavirus               | Has 15 mat_peptides; Ver1.1.4, species
# 572289 => 'NC_011550', #    694014 | Avian coronavirus               | No mat_peptide in refseq
# 572290 => 'NC_011549', #    694014 | Avian coronavirus               | No mat_peptide in refseq
# 694015 => 'NC_010646', #    694015 | Beluga Whale coronavirus SW1    | No mat_peptide in refseq
# 864596 => 'NC_014470', #    864596 | Bat coronavirus BM48-31/BGR/2008| No mat_peptide in refseq
#1159902 => 'NC_016996', #   1159902 | Common-moorhen coronavirus HKU21| No mat_peptide in refseq
#1159903 => 'NC_016993', #   1159903 | Magpie-robin coronavirus HKU18  | No mat_peptide in refseq
#1159904 => 'NC_016994', #   1159904 | Night-heron coronavirus HKU19   | No mat_peptide in refseq
#1159905 => 'NC_016990', #   1159905 | Porcine coronavirus HKU15       | No mat_peptide in refseq
#1159906 => 'NC_016992', #   1159906 | Sparrow coronavirus HKU17       | No mat_peptide in refseq
#1159907 => 'NC_016991', #   1159907 | White-eye coronavirus HKU16     | No mat_peptide in refseq
#1159908 => 'NC_016995', #   1159908 | Wigeon coronavirus HKU20        | No mat_peptide in refseq
#1160968 => 'NC_017083', #   1160968 | Rabbit coronavirus HKU14        | No mat_peptide in refseq
#+---------+-----------+-----------+---------------------------------------------------------+
           },

           # Family=Arenaviridae
           'Arenaviridae' => {
#  11623 => 'NC_004291', # Lymphocytic choriomeningitis virus |     11623 | No mat_peptide in refseq
#  11623 => 'NC_004294', # Lymphocytic choriomeningitis virus |     11623 | No mat_peptide in refseq
#  11631 => 'NC_004292', # Tacaribe virus                     |     11631 | No mat_peptide in refseq
#  11631 => 'NC_004293', # Tacaribe virus                     |     11631 | No mat_peptide in refseq
#  11620 => 'NC_004296', # Lassa virus                        |     11620 | No mat_peptide in refseq
#  11620 => 'NC_004297', # Lassa virus                        |     11620 | No mat_peptide in refseq
#  45219 => 'NC_005077', # Guanarito virus                    |     45219 | No mat_peptide in refseq
#  45219 => 'NC_005082', # Guanarito virus                    |     45219 | No mat_peptide in refseq
#  11628 => 'NC_005078', # Machupo virus                      |     11628 | No mat_peptide in refseq
#  11628 => 'NC_005079', # Machupo virus                      |     11628 | No mat_peptide in refseq
#  11619 => 'NC_005080', # Junin virus                        |     11619 | No mat_peptide in refseq
  11619 => 'NC_005081', # Junin virus                        |     11619 | Good refseq w/ 2 mat_peptides; Ver1.1.4
#  49891 => 'NC_005894', # Pirital virus                      |     49891 | No mat_peptide in refseq
#  49891 => 'NC_005897', # Pirital virus                      |     49891 | No mat_peptide in refseq
#  45709 => 'NC_006313', # Sabia virus                        |     45709 | No mat_peptide in refseq
#  45709 => 'NC_006317', # Sabia virus                        |     45709 | No mat_peptide in refseq
#  11630 => 'NC_006439', # Pichinde virus                     |     11630 | No mat_peptide in refseq
#  11630 => 'NC_006447', # Pichinde virus                     |     11630 | No mat_peptide in refseq
# 300180 => 'NC_006572', # Lassa virus                        |     11620 | No mat_peptide in refseq
# 300180 => 'NC_006573', # Lassa virus                        |     11620 | No mat_peptide in refseq
# 300175 => 'NC_006574', # Mopeia virus                       |     11629 | No mat_peptide in refseq
# 300175 => 'NC_006575', # Mopeia virus                       |     11629 | No mat_peptide in refseq
#  55097 => 'NC_007903', # Mobala virus                       |     55097 | No mat_peptide in refseq
#  55097 => 'NC_007904', # Mobala virus                       |     55097 | No mat_peptide in refseq
#  55096 => 'NC_007905', # Ippy virus                         |     55096 | No mat_peptide in refseq
#  55096 => 'NC_007906', # Ippy virus                         |     55096 | No mat_peptide in refseq
#  45218 => 'NC_010247', # Amapari virus                      |     45218 | No mat_peptide in refseq
#  45218 => 'NC_010251', # Amapari virus                      |     45218 | No mat_peptide in refseq
#  42764 => 'NC_010248', # Oliveros virus                     |     42764 | No mat_peptide in refseq
#  42764 => 'NC_010250', # Oliveros virus                     |     42764 | No mat_peptide in refseq
# 144752 => 'NC_010249', # Allpahuayo virus                   |    144752 | No mat_peptide in refseq
# 144752 => 'NC_010253', # Allpahuayo virus                   |    144752 | No mat_peptide in refseq
# 208899 => 'NC_010252', # Cupixi virus                       |    208899 | No mat_peptide in refseq
# 208899 => 'NC_010254', # Cupixi virus                       |    208899 | No mat_peptide in refseq
# 192848 => 'NC_010255', # Bear Canyon virus                  |    192848 | No mat_peptide in refseq
# 192848 => 'NC_010256', # Bear Canyon virus                  |    192848 | No mat_peptide in refseq
# 499556 => 'NC_010562', # Chapare virus                      |    499556 | No mat_peptide in refseq
# 499556 => 'NC_010563', # Chapare virus                      |    499556 | No mat_peptide in refseq
#  46919 => 'NC_010700', # Whitewater Arroyo virus            |     46919 | No mat_peptide in refseq
#  46919 => 'NC_010703', # Whitewater Arroyo virus            |     46919 | No mat_peptide in refseq
#  45223 => 'NC_010701', # Tamiami virus                      |     45223 | No mat_peptide in refseq
#  45223 => 'NC_010702', # Tamiami virus                      |     45223 | No mat_peptide in refseq
#  45222 => 'NC_010756', # Parana virus                       |     45222 | No mat_peptide in refseq
#  45222 => 'NC_010761', # Parana virus                       |     45222 | No mat_peptide in refseq
#  45220 => 'NC_010757', # Flexal virus                       |     45220 | No mat_peptide in refseq
#  45220 => 'NC_010759', # Flexal virus                       |     45220 | No mat_peptide in refseq
#  45221 => 'NC_010758', # Latino virus                       |     45221 | No mat_peptide in refseq
#  45221 => 'NC_010760', # Latino virus                       |     45221 | No mat_peptide in refseq
# 649188 => 'NC_012776', # Lujo virus                         |    649188 | No mat_peptide in refseq
# 649188 => 'NC_012777', # Lujo virus                         |    649188 | No mat_peptide in refseq
# 573900 => 'NC_013057', # Morogoro virus                     |    573900 | No mat_peptide in refseq
# 573900 => 'NC_013058', # Morogoro virus                     |    573900 | No mat_peptide in refseq
            },

  # Family=Bunyaviridae
  'Bunyaviridae' => {
# 992212 => 'NC_018136', # Severe fever with thrombocytopenia syndrome virus | 1003835 | No mat_peptide in refseq
# 992212 => 'NC_018137', # Severe fever with thrombocytopenia syndrome virus | 1003835 | No mat_peptide in refseq
 992212 => 'NC_018138', # Severe fever with thrombocytopenia syndrome virus | 1003835 | Good refseq, w/ 2 mat_peptide + 2 sig_peptide, V1.1.7
#  35304 => 'NC_001925', # Bunyamwera virus                                |  35304 | No mat_peptide in refseq
#  35304 => 'NC_001926', # Bunyamwera virus                                |  35304 | No mat_peptide in refseq
#  35304 => 'NC_001927', # Bunyamwera virus                                |  35304 | No mat_peptide in refseq
#  11588 => 'NC_002043', # Rift Valley fever virus                         |  11588 | Replaced by NC_014397
#  11588 => 'NC_002044', # Rift Valley fever virus                         |  11588 | Replaced by NC_014396
#  11588 => 'NC_002045', # Rift Valley fever virus                         |  11588 | Replaced by NC_014395
#  46607 => 'NC_003466', # Andes virus                                     |  46607 | No mat_peptide in refseq
#  46607 => 'NC_003467', # Andes virus                                     |  46607 | No mat_peptide in refseq
#  46607 => 'NC_003468', # Andes virus                                     |  46607 | No mat_peptide in refseq
#  11577 => 'NC_004108', # California encephalitis virus                   |  11577 | No mat_peptide in refseq
#  11577 => 'NC_004109', # California encephalitis virus                   |  11577 | No mat_peptide in refseq
#  11577 => 'NC_004110', # California encephalitis virus                   |  11577 | No mat_peptide in refseq
#  11595 => 'NC_004157', # Dugbe virus                                     |  11595 | No mat_peptide in refseq
#  11595 => 'NC_004158', # Dugbe virus                                     |  11595 | Has only 1 mat_peptide, ignore for V1.1.4
#  11595 => 'NC_004159', # Dugbe virus                                     |  11595 | No mat_peptide in refseq
#  11591 => 'NC_005214', # Uukuniemi virus                                 |  11591 | No mat_peptide in refseq
  11591 => 'NC_005220', # Uukuniemi virus                                 |  11591 | Good refseq, w/ 2 mat_peptide + 2 sig_peptide, V1.1.4
#  11591 => 'NC_005221', # Uukuniemi virus                                 |  11591 | No mat_peptide in refseq
#  37705 => 'NC_005215', # Sin Nombre virus                                |  37705 | No mat_peptide in refseq
#  37705 => 'NC_005216', # Sin Nombre virus                                |  37705 | No mat_peptide in refseq
#  37705 => 'NC_005217', # Sin Nombre virus                                |  37705 | No mat_peptide in refseq
#  11599 => 'NC_005218', # Hantaan virus                                   |  11599 | No mat_peptide in refseq
  11599 => 'NC_005219', # Hantaan virus                                   |  11599 | Good refseq, w/ 2 mat_peptide + 1 sig_peptide, V1.1.4
#  11599 => 'NC_005222', # Hantaan virus                                   |  11599 | No mat_peptide in refseq
#  11604 => 'NC_005223', # Puumala virus                                   |  11604 | No mat_peptide in refseq
#  11604 => 'NC_005224', # Puumala virus                                   |  11604 | No mat_peptide in refseq
#  11604 => 'NC_005225', # Puumala virus                                   |  11604 | No mat_peptide in refseq
#  37133 => 'NC_005226', # Tula virus                                      |  37133 | No mat_peptide in refseq
#  37133 => 'NC_005227', # Tula virus                                      |  37133 | No mat_peptide in refseq
#  37133 => 'NC_005228', # Tula virus                                      |  37133 | No mat_peptide in refseq
#  12506 => 'NC_005233', # Dobrava-Belgrade virus                          |  12506 | No mat_peptide in refseq
#  12506 => 'NC_005234', # Dobrava-Belgrade virus                          |  12506 | No mat_peptide in refseq
#  12506 => 'NC_005235', # Dobrava-Belgrade virus                          |  12506 | No mat_peptide in refseq
#  11608 => 'NC_005236', # Seoul virus                                     |  11608 | No mat_peptide in refseq
#  11608 => 'NC_005237', # Seoul virus                                     |  11608 | No mat_peptide in refseq
#  11608 => 'NC_005238', # Seoul virus                                     |  11608 | No mat_peptide in refseq
#  11593 => 'NC_005300', # Crimean-Congo hemorrhagic fever virus           |  11593 | No mat_peptide in refseq
#  11593 => 'NC_005301', # Crimean-Congo hemorrhagic fever virus           |  11593 | No mat_peptide in refseq
#  11593 => 'NC_005302', # Crimean-Congo hemorrhagic fever virus           |  11593 | No mat_peptide in refseq
# 118655 => 'NC_005775', # Oropouche virus                                 | 118655 | No mat_peptide in refseq
# 118655 => 'NC_005776', # Oropouche virus                                 | 118655 | No mat_peptide in refseq
# 118655 => 'NC_005777', # Oropouche virus                                 | 118655 | No mat_peptide in refseq
# 206160 => 'NC_006318', # Sandfly fever Naples virus                      | 206160 | No mat_peptide in refseq
# 206160 => 'NC_006319', # Sandfly fever Naples virus                      | 206160 | No mat_peptide in refseq
# 206160 => 'NC_006320', # Sandfly fever Naples virus                      | 206160 | No mat_peptide in refseq
#  93830 => 'NC_006433', # Z10                                             |  93830 | No mat_peptide in refseq
#  93830 => 'NC_006435', # Z10                                             |  93830 | No mat_peptide in refseq
#  93830 => 'NC_006437', # Z10                                             |  93830 | No mat_peptide in refseq
#  70566 => 'NC_009894', # Akabane virus                                   |  70566 | No mat_peptide in refseq
#  70566 => 'NC_009895', # Akabane virus                                   |  70566 | No mat_peptide in refseq
#  70566 => 'NC_009896', # Akabane virus                                   |  70566 | No mat_peptide in refseq
# 262967 => 'NC_010704', # Thottapalayam virus                             | 262967 | No mat_peptide in refseq
# 262967 => 'NC_010707', # Thottapalayam virus                             | 262967 | No mat_peptide in refseq
# 262967 => 'NC_010708', # Thottapalayam virus                             | 262967 | No mat_peptide in refseq
# 267407 => 'NC_013106', # European mountain ash ringspot-associated virus | 267407 | No mat_peptide in refseq
# 267407 => 'NC_013107', # European mountain ash ringspot-associated virus | 267407 | No mat_peptide in refseq
# 267407 => 'NC_013108', # European mountain ash ringspot-associated virus | 267407 | No mat_peptide in refseq
         },

  # Family=Filoviridae
#  11269 => 'NC_001608', # Lake Victoria marburgvirus |  11269 | No mat_peptide in refseq
# 186538 => 'NC_002549', # Zaire ebolavirus           | 186538 | No mat_peptide in refseq
# 186539 => 'NC_004161', # Reston ebolavirus          | 186539 | No mat_peptide in refseq
# 186540 => 'NC_006432', # Sudan ebolavirus           | 186540 | No mat_peptide in refseq

  # Family=Paramyxoviridae
#mysql> select taxid,accession,species,taxid from genome left join taxon on taxid=taxon.id where accession like "NC_%" and taxid in (select id from taxon where family="Paramyxoviridae") ;
#  11234 => 'NC_001498', # Measles virus                      |  11234 | No mat_peptide in refseq
#  11191 => 'NC_001552', # Sendai virus                       |  11191 | No mat_peptide in refseq
#  11250 => 'NC_001781', # Human respiratory syncytial virus  |  11250 | No mat_peptide in refseq
#  11216 => 'NC_001796', # Human parainfluenza virus 3        |  11216 | No mat_peptide in refseq
#  12814 => 'NC_001803', # Respiratory syncytial virus        |  12814 | No mat_peptide in refseq
#  63330 => 'NC_001906', # Hendra virus                       |  63330 | No mat_peptide in refseq
#  11232 => 'NC_001921', # Canine distemper virus             |  11232 | No mat_peptide in refseq
#  11246 => 'NC_001989', # Bovine respiratory syncytial virus |  11246 | No mat_peptide in refseq
#  11215 => 'NC_002161', # Bovine parainfluenza virus 3       |  11215 | No mat_peptide in refseq
#  92129 => 'NC_002199', # Tupaia paramyxovirus               |  92129 | No mat_peptide in refseq
#  11161 => 'NC_002200', # Mumps virus                        |  11161 | No mat_peptide in refseq
# 139270 => 'NC_002617', # Newcastle disease virus            | 139270 | No mat_peptide in refseq
# 121791 => 'NC_002728', # Nipah virus                        | 121791 | No mat_peptide in refseq
# 157619 => 'NC_003043', # Avian paramyxovirus 6              | 157619 | No mat_peptide in refseq
#  11212 => 'NC_003443', # Human parainfluenza virus 2        |  11212 | No mat_peptide in refseq
#  12730 => 'NC_003461', # Human parainfluenza virus 1        |  12730 | No mat_peptide in refseq
# 162013 => 'NC_004074', # Tioman virus                       | 162013 | No mat_peptide in refseq
# 162145 => 'NC_004148', # Human metapneumovirus              | 162145 | No mat_peptide in refseq
# 204987 => 'NC_005036', # Goose paramyxovirus SF02           | 204987 | No mat_peptide in refseq
# 122203 => 'NC_005084', # Fer-de-lance virus                 | 122203 | No mat_peptide in refseq
#  37131 => 'NC_005283', # Cetacean morbillivirus             |  37131 | No mat_peptide in refseq
# 241630 => 'NC_005339', # Mossman virus                      | 241630 | No mat_peptide in refseq
#  11242 => 'NC_006296', # Rinderpest virus                   |  11242 | No mat_peptide in refseq
#  31604 => 'NC_006383', # Peste-des-petits-ruminants virus   |  31604 | No mat_peptide in refseq
#  11228 => 'NC_006428', # Simian virus 41                    |  11228 | No mat_peptide in refseq
#  11207 => 'NC_006430', # Simian virus 5                     |  11207 | No mat_peptide in refseq
# 270473 => 'NC_006579', # Murine pneumonia virus             | 270473 | No mat_peptide in refseq
# 322067 => 'NC_007454', # J-virus                            | 322067 | No mat_peptide in refseq
# 152219 => 'NC_007620', # Menangle virus                     | 152219 | No mat_peptide in refseq
#  38525 => 'NC_007652', # Avian metapneumovirus              |  38525 | No mat_peptide in refseq
# 341053 => 'NC_007803', # Beilong virus                      | 341053 | No mat_peptide in refseq
#  43140 => 'NC_009489', # Mapuera virus                      |  43140 | No mat_peptide in refseq
#  53179 => 'NC_009640', # Porcine rubulavirus                |  53179 | No mat_peptide in refseq

  # Family=Hepeviridae
  'Hepeviridae' => {
#mysql> select taxid,accession,species,taxid from genome left join taxon on taxid=taxon.id where accession like "NC_%" and taxid in (select id from taxon where family="Hepeviridae") ;
#+---------+-----------+-----------------------+---------+
#| taxid   | accession | species               | taxid   |
#+---------+-----------+-----------------------+---------+
#    12461 => 'NC_001434', # Hepatitis E virus     |   12461 | No mat_peptide in refseq
#  1016879 => 'NC_015521', # Cutthroat trout virus | 1016879 |  Good refseq w/ 4 mat_peptides
#+---------+-----------+-----------------------+---------+
#2 rows in set (0.01 sec)
    },

  # Family=Picornaviridae # V1.1.6
  'Picornaviridae' => {
#mysql> select taxid,accession,speciesid,species from genome left join taxon on taxid=taxon.id where accession like "NC_%" and taxid in (select id from taxon where family="Picornaviridae") order by speciesid,taxid,accession;
#+---------+-----------+-----------+----------------------------------+
#| taxid   | accession | speciesid | species                          |
#+---------+-----------+-----------+----------------------------------+
   12064 => 'NC_001859', #   12064 | Bovine enterovirus               | Good refseq w/ 11 mat_peptides, V1.1.6
##  269638 => 'NC_017676', #   12064 | Bovine enterovirus               |  No mat_peptide in refseq
   12092 => 'NC_001489', #   12092 | Hepatitis A virus                | Good refseq w/ 12 mat_peptides, V1.1.6
   12104 => 'NC_001479', #   12104 | Encephalomyocarditis virus       | Good refseq w/ 11 mat_peptides, V1.1.6
##   12111 => 'NC_011450', #   12110 | Foot-and-mouth disease virus     | No mat_peptide in refseq
 # The VP1-2A cleavage site doesn't look good in 12116 => 'NC_002554', still used for subspecies 12116
   12116 => 'NC_002554', #   12110 | Foot-and-mouth disease virus     | Good refseq w/ 14 mat_peptides, V1.1.6
 # All cleavage sites seem good in 12118 => 'NC_004004', taking as refseq for entire species 12110
   12110 => 'NC_004004', #   12110 | Foot-and-mouth disease virus     | Good refseq w/ 12 mat_peptides, V1.1.6
   12118 => 'NC_004004', #   12110 | Foot-and-mouth disease virus     | Good refseq w/ 12 mat_peptides, V1.1.6
##   12122 => 'NC_011451', #   12110 | Foot-and-mouth disease virus     | No mat_peptide in refseq
##   12123 => 'NC_011452', #   12110 | Foot-and-mouth disease virus     | No mat_peptide in refseq
##   35292 => 'NC_003992', #   12110 | Foot-and-mouth disease virus     | No mat_peptide in refseq
 # Problem: 5 proteins (P1-2A) are lumped together in 110195 => 'NC_004915'
##  110195 => 'NC_004915', #   12110 | Foot-and-mouth disease virus     | Good refseq w/ 7 mat_peptides
   47000 => 'NC_003982', #   47000 | Equine rhinitis A virus          | Good refseq w/ 12 mat_peptides, V1.1.6
   70796 => 'NC_003990', #   70796 | Avian encephalomyelitis virus    | Good refseq w/ 11 mat_peptides, V1.1.6
   72149 => 'NC_001918', #   72149 | Aichi virus                      | Good refseq w/ 10 mat_peptides, V1.1.6
  106966 => 'NC_004441', #  106966 | Porcine enterovirus B            | Good refseq w/ 11 mat_peptides, V1.1.6
  118140 => 'NC_003985', #  118140 | Porcine teschovirus              | Good refseq w/ 12 mat_peptides, V1.1.6
  138948 => 'NC_001612', #  138948 | Human enterovirus A              | Good refseq w/ 11 mat_peptides, V1.1.6
  442855 => 'NC_010412', #  138948 | Human enterovirus A              | Good refseq w/ 11 mat_peptides, V1.1.6
  442858 => 'NC_010413', #  138948 | Human enterovirus A              | Good refseq w/ 11 mat_peptides, V1.1.6
  138949 => 'NC_001472', #  138949 | Human enterovirus B              | Good refseq w/ 11 mat_peptides, V1.1.6
  172038 => 'NC_010411', #  138949 | Human enterovirus B              | Good refseq w/ 11 mat_peptides, V1.1.6
  325447 => 'NC_009887', #  138949 | Human enterovirus B              | Good refseq w/ 11 mat_peptides, V1.1.6
  582384 => 'NC_013114', #  138949 | Human enterovirus B              | Good refseq w/ 11 mat_peptides, V1.1.6
 # Following 10 are required to avoid gap at cleavage site for those species
   12060 => 'NC_013114', #  138949 | Human enterovirus B              | Good refseq w/ 11 mat_peptides, V1.1.6
   12062 => 'NC_013114', #  138949 | Human enterovirus B              | Good refseq w/ 11 mat_peptides, V1.1.6
   12067 => 'NC_013114', #  138949 | Human enterovirus B              | Good refseq w/ 11 mat_peptides, V1.1.6
   47506 => 'NC_013114', #  138949 | Human enterovirus B              | Good refseq w/ 11 mat_peptides, V1.1.6
  103915 => 'NC_013114', #  138949 | Human enterovirus B              | Good refseq w/ 11 mat_peptides, V1.1.6
  222887 => 'NC_013114', #  138949 | Human enterovirus B              | Good refseq w/ 11 mat_peptides, V1.1.6
  222889 => 'NC_013114', #  138949 | Human enterovirus B              | Good refseq w/ 11 mat_peptides, V1.1.6
  318562 => 'NC_013114', #  138949 | Human enterovirus B              | Good refseq w/ 11 mat_peptides, V1.1.6
  318563 => 'NC_013114', #  138949 | Human enterovirus B              | Good refseq w/ 11 mat_peptides, V1.1.6
  325446 => 'NC_013114', #  138949 | Human enterovirus B              | Good refseq w/ 11 mat_peptides, V1.1.6
  582385 => 'NC_013115', #  138949 | Human enterovirus B              | Good refseq w/ 11 mat_peptides, V1.1.6
#  138950 => 'NC_001428', #  138950 | Human enterovirus C              | Good refseq w/ 11 mat_peptides, "PROVISIONAL"
  138950 => 'NC_002058', #  138950 | Human enterovirus C              | Good refseq w/ 11 mat_peptides, "REVIEWED", V1.1.6
  861519 => 'NC_014336', #  138950 | Human enterovirus C              | Good refseq w/ 11 mat_peptides, V1.1.6
  138951 => 'NC_001430', #  138951 | Human enterovirus D              | Good refseq w/ 15 mat_peptides, V1.1.6
  147711 => 'NC_001617', #  147711 | Human rhinovirus A               | Good refseq w/ 11 mat_peptides, V1.1.6
   12131 => 'NC_001490', #  147712 | Human rhinovirus B               | Good refseq w/ 11 mat_peptides, V1.1.6
  147712 => 'NC_001490', #  147712 | Human rhinovirus B               | Good refseq w/ 11 mat_peptides, V1.1.6
  172314 => 'NC_003976', #  172314 | Ljungan virus                    | Good refseq w/ 11 mat_peptides, V1.1.6
  184753 => 'NC_010384', #  184753 | Simian picornavirus strain N125  | Good refseq w/ 11 mat_peptides, V1.1.6
  184754 => 'NC_013695', #  184754 | Simian picornavirus strain N203  | Good refseq w/ 11 mat_peptides, V1.1.6
  194965 => 'NC_004421', #  194965 | Bovine kobuvirus                 | Good refseq w/ 11 mat_peptides, V1.1.6
  195054 => 'NC_001897', #  195054 | Human parechovirus               | Good refseq w/ 10 mat_peptides, V1.1.6
  204711 => 'NC_001366', #  204711 | Theilovirus                      | Good refseq w/ 12 mat_peptides, V1.1.6
  434309 => 'NC_009448', #  204711 | Theilovirus                      | Good refseq w/ 12 mat_peptides, V1.1.6
  511755 => 'NC_010810', #  204711 | Theilovirus                      | Good refseq w/ 12 mat_peptides, V1.1.6
##  263531 => 'NC_008714', #  263531 | Possum enterovirus W1            | No mat_peptide in refseq
##  263532 => 'NC_008715', #  263532 | Possum enterovirus W6            | No mat_peptide in refseq
  310907 => 'NC_003988', #  310907 | Simian enterovirus A             | Good refseq w/ 11 mat_peptides, V1.1.6
 # There are problems with following refseq (NC_003983) for Equine rhinitis B virus, so the species is left out
##   47001 => 'NC_003983', #  312185 | Equine rhinitis B virus          | Good refseq w/ 12 mat_peptides
 # There are problems with following refseq (NC_003077) for Equine rhinitis B virus, so the species is left out
##  168014 => 'NC_003077', #  312185 | Equine rhinitis B virus          | Good refseq w/ 12 mat_peptides
##  390157 => 'NC_011349', #  390157 | Seneca valley virus              | No mat_peptide in refseq
  442860 => 'NC_010415', #  442860 | Simian enterovirus SV6           | Good refseq w/ 11 mat_peptides, V1.1.6
  463676 => 'NC_009996', #  463676 | Human rhinovirus C               | Good refseq w/ 11 mat_peptides, V1.1.6
  471728 => 'NC_009891', #  471728 | Seal picornavirus type 1         | Good refseq w/ 11 mat_peptides, V1.1.6
  586419 => 'NC_012800', #  586419 | Human cosavirus A                | Good refseq w/ 11 mat_peptides, V1.1.6
  586420 => 'NC_012801', #  586420 | Human cosavirus B                | Good refseq w/ 11 mat_peptides, V1.1.6
  586422 => 'NC_012802', #  586422 | Human cosavirus D                | Good refseq w/ 11 mat_peptides, V1.1.6
  586423 => 'NC_012798', #  586423 | Human cosavirus E                | Good refseq w/ 11 mat_peptides, V1.1.6
##  655603 => 'NC_012986', #  655603 | Human klassevirus 1              | No mat_peptide in refseq
  686983 => 'NC_006553', #  686983 | Avian sapelovirus                | Good refseq w/ 12 mat_peptides, V1.1.6
 1002921 => 'NC_003987', #  686984 | Porcine sapelovirus              | Good refseq w/ 11 mat_peptides, V1.1.6
  686984 => 'NC_003987', #  686984 | Porcine sapelovirus              | Good refseq w/ 11 mat_peptides, V1.1.6
 # Species is Simian sapelovirus (686985)
 1002918 => 'NC_004451', #  686985 | Simian sapelovirus               | Good refseq w/ 12 mat_peptides, V1.1.6
  686985 => 'NC_004451', #  686985 | Simian sapelovirus               | Good refseq w/ 12 mat_peptides, V1.1.6
  651733 => 'NC_012957', #  688449 | Salivirus                        | Good refseq w/ 11 mat_peptides, V1.1.6
## 1006061 => 'NC_008250', #  691956 | Duck hepatitis A virus           | No mat_peptide in refseq
## 1006063 => 'NC_009750', #  691956 | Duck hepatitis A virus           | No mat_peptide in refseq
  693066 => 'NC_010354', #  693066 | Bovine rhinitis B virus          | Good refseq w/ 12 mat_peptides, V1.1.6
  871699 => 'NC_014411', #  871699 | Turdivirus 1                     | Good refseq w/ 11 mat_peptides, V1.1.6
  871700 => 'NC_014412', #  871700 | Turdivirus 2                     | Good refseq w/ 11 mat_peptides, V1.1.6
  871701 => 'NC_014413', #  871701 | Turdivirus 3                     | Good refseq w/ 11 mat_peptides, V1.1.6
  928289 => 'NC_015626', #  928289 | Pigeon picornavirus B            | Good refseq w/ 12 mat_peptides, V1.1.6
## 1074213 => 'NC_015936', # 1074213 | Mouse kobuvirus M-5/USA/2010     | No mat_peptide in refseq
 1074863 => 'NC_015940', # 1074863 | Bat picornavirus 1               | Good refseq w/ 14 mat_peptides, V1.1.6
 1074864 => 'NC_015941', # 1074864 | Bat picornavirus 2               | Good refseq w/ 14 mat_peptides, V1.1.6
 1074865 => 'NC_015934', # 1074865 | Bat picornavirus 3               | Good refseq w/ 12 mat_peptides, V1.1.6
## 1089138 => 'NC_016403', # 1089138 | Quail picornavirus QPV1/HUN/2010 | No mat_peptide in refseq
 1108810 => 'NC_016156', # 1108810 | Feline picornavirus              | Good refseq w/ 12 mat_peptides, V1.1.6
 # NC_011829 is not a high quality refseqs, therefore it's not used for the species Porcine kobuvirus (1156769)
  569195 => 'NC_011829', # 1156769 | Porcine kobuvirus                | Good refseq w/ 11 mat_peptides, V1.1.6
## 1136133 => 'NC_016769', # 1156769 | Porcine kobuvirus                | No mat_peptide in refseq
 1210913 => 'NC_018226', # 1210913 | Swine pasivirus 1                | Good refseq w/ 11 mat_peptides, V1.1.6
## 1214230 => 'NC_018400', # 1214230 | Turkey gallivirus                | No mat_peptide in refseq
 1233320 => 'NC_018668', # 1233320 | Bovine hungarovirus              | Good refseq w/ 12 mat_peptides, V1.1.6
#+---------+-----------+-----------+----------------------------------+
#73 rows in set (0.01 sec)
    },
           };

    # This variable is used to hold families in development
    my $tmp_ref = {
        };
    $debug && print STDERR "$subn: \$tmp_ref=\n".Dumper($tmp_ref)."end of \$tmp_ref\n\n";

    if ( 1 ) {
        for my $k (keys %$tmp_ref) {
            $REFSEQS->{refs}->{$k} = $tmp_ref->{$k};
        }
    }

    $debug && print STDERR "$subn: \$REFSEQS=\n".Dumper($REFSEQS)."end of \$REFSEQS\n\n";

    # Determine total numbers of families, species, strains
    my $nstrain = 0;
    for my $fam (keys %{$REFSEQS->{refs}}) {
        $nstrain += scalar keys %{$REFSEQS->{refs}->{$fam}};
        $debug && print STDERR "$subn: \$fam=$fam\t\$nstrain=$nstrain\n";
    }

    return $REFSEQS;
} # sub initRefseq


=head2 get_gene_symbol
Takes a mat_peptide feature from refseq object, returns the gene symbol. The gene symbols are
 defined as hash (accession) of hash (location of matpeptide) here, since the values in genbank
 file are inconsistent or unavailable for every mat_peptide
Returns the gene symbol in a string
=cut

sub get_gene_symbol {
    my ($reffeat,$exe_dir) = @_;

    my $debug = 0 && $debug_all;
    my $subn = 'get_gene_symbol';
#    print STDERR "get_gene_symbol: \$reffeat=\n".Dumper($reffeat)."end of \$reffeat\n\n";

    # These gene symbols are defined by CLarsen, after considering refseqs, e.g. NC_001477 & NC_009942
    my $gene_symbols;

#my $GENE_SYM = {
#               symbol_loaded => 0,
#               symbol_fn     => "Annotate_symbol_records.txt",
#            };
    # Read the definition of symbols from the text file
    $debug && print STDERR "$subn: \$exe_dir=$exe_dir\n";
    $debug && print STDERR "$subn: \$GENE_SYM->{symbol_fn}=$GENE_SYM->{symbol_fn}\n";
    if (!$GENE_SYM->{symbol_loaded} && -f "$exe_dir/$GENE_SYM->{symbol_fn}" ) {
        $debug && print STDERR "$subn: \$GENE_SYM->{symbol_loaded}=$GENE_SYM->{symbol_loaded}\n";
        $debug && print STDERR "$subn: found symbol_fn=$exe_dir/$GENE_SYM->{symbol_fn}\n";

        open my $symbol_file, '<', "$exe_dir/$GENE_SYM->{symbol_fn}"
           or croak("$subn: found '$exe_dir/$GENE_SYM->{symbol_fn}', but couldn't open: $OS_ERROR");
        while (<$symbol_file>) {
            chomp;
#            $debug && print STDERR "$subn: \$_='$_'\n";
            next if (m/^[#]/x);  # Skip any comment
            next if (m/^\s*$/x); # Skip empty lines
            s/'//g; # Remove all single quotes
            my $words = [split(/\s*;\s*/)];
#            $debug && print STDERR "$subn: \$_='$_'\n";
#            $debug && print STDERR "$subn: \$words($#{$words})='@$words'\n";
            if ($#{$words}<2) { # skip the lines without enough fields
                $debug && print STDERR "$subn: not enough data in gene_symbol: '@$words'\n";
                next;
            }

            my ($acc, $loc, $sym) = @{$words}[0..3];
#            $loc = '0..0' if (!$loc);
            next if (!$loc); # Skip if no location
            next if (!$sym); # Skip if no symbol
            $debug && print STDERR "$subn: \$acc=$acc \$loc=$loc \$sym=$sym\n";
            if (exists($GENE_SYM->{symbol}->{$acc}->{$loc}) && $sym ne $GENE_SYM->{symbol}->{$acc}->{$loc}) {
                print STDERR "$subn: conflicting data, \$GENE_SYM->{symbol}->{$acc}->{$loc}=$GENE_SYM->{symbol}->{$acc}->{$loc}\n";
                print STDERR "$subn: conflicting data, new \$sym=$sym\n";
            } else {
                $GENE_SYM->{symbol}->{$acc}->{$loc} = $sym;
            }
        }
        close $symbol_file or croak "$subn: Couldn't close $exe_dir/$GENE_SYM->{symbol_fn}: $OS_ERROR";
        $GENE_SYM->{symbol_loaded} = 1;

        $debug && print STDERR "$subn: finished reading list file: '$exe_dir/$GENE_SYM->{symbol_fn}'.\n";
        $debug && print STDERR "$subn: \$GENE_SYM->{symbol_loaded}=$GENE_SYM->{symbol_loaded}\n";
        $debug && print STDERR "$subn: \$GENE_SYM=\n".Dumper($GENE_SYM)."End of \$GENE_SYM\n\n";
    }

    # Search for the symbol based on the start/stop of the reference mat_peptide
    my $gene_symbol;
    $gene_symbol = ''; # Make the default value of gene_symbol as ''
    my $acc = $reffeat->seq->accession_number;
    my $loc = $reffeat->location->to_FTstring;
    $loc = $reffeat->location->start .'..'. $reffeat->location->end;
    $gene_symbol = $GENE_SYM->{symbol}->{$acc}->{$loc} ? $GENE_SYM->{symbol}->{$acc}->{$loc} : $gene_symbol;
    $debug && print "get_gene_symbol: primary_tag=".($reffeat->seq->accession_number);
    $debug && print "\tlocation=".($reffeat->location->to_FTstring);
    $debug && print "\tsymbol=$gene_symbol\n";

    return $gene_symbol;
} # sub get_gene_symbol


=head2 check_refseq
Takes 1 refseq and 1 target genome. Check if their taxids are same
=cut

sub check_refseq {
    my ($refseq, $inseq) = @_;

    my $debug = 0 && $debug_all;

    my $refseq_fit = 0;
    return $refseq_fit if (!$refseq || !$inseq);

    my $taxid = $inseq->species->ncbi_taxid;
    my $ref_taxid = $refseq->species->ncbi_taxid;
    if ($taxid == $ref_taxid) {
      $refseq_fit = 1;
    }

    return $refseq_fit;
} # sub check_refseq


=head2 similarAA
  Takes 2 AA, decides if they are similar. Returns 1 if similar
=cut

sub similarAA {
    my ($rc, $tc) = @_;

    my $debug = 0 && $debug_all;
    my $subn = 'similarAA';

    my $AA_groups = {
                'A' => 'AG',
                'R' => 'KR',
                'N' => 'NQ',
                'D' => 'DE',
                'C' => 'C',
                'Q' => 'NQ',
                'E' => 'DE',
                'G' => 'AG',
                'H' => 'H',
                'I' => 'ILV',
                'L' => 'ILV',
                'K' => 'KR',
                'M' => 'M',
                'F' => 'FY',
                'P' => 'P',
                'S' => 'ST',
                'T' => 'ST',
                'W' => 'W',
                'Y' => 'FY',
                'V' => 'ILV',
                };

    my $similar = 0;

    $rc = uc $rc;
    $tc = uc $tc;
    $similar = 1 if ($AA_groups->{$rc} eq $AA_groups->{$tc});

    return $similar;
} # sub similarAA


1;
