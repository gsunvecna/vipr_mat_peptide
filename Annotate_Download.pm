package Annotate_Download;

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

use version; our $VERSION = qv('1.1.6'); # Nov 07 2012

use Annotate_Def;

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


## //subroutines// ##

=head2
 sub downloadRefseq searches genbank for all refseqs for a family, then downloads all the refseqs
=cut

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


=head2
 sub downloadTaxon for a family downloads the the taxonnomy information for all strains/species
=cut

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


=head2
 sub loadXmlTaxon reads an XML file, searchs for each taxid, and its name, speciesid, species,
 genus, genusid, family, and familyid, to be compared with exsting data
=cut

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


=head2
 sub getTaxonString takes a hash of a taxon, and finds the relative info, and returns an array
=cut

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


=head2
 sub xml2hash takes an array of strings of XML, turns it to a hash
=cut

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


=head2
 sub checkAllTaxon looks over each family and see if there is any update in its taxonomy from NCBI
=cut

sub checkAllTaxon {
    my ( $exe_dir) = @_;

    my $debug = 0 && $debug_all;
    my $subn = 'checkAllTaxon';

    # Load the $REFSEQS if not already
    if ( !$TAXON->{taxon_loaded} ) {
        $TAXON = Annotate_Def::loadTaxonTable( $exe_dir);
        my $ctLines = scalar(keys %{$TAXON->{taxon}});
        $debug && print STDERR "$subn: loaded $ctLines lines from file '$exe_dir/$TAXON->{taxon_fn}'\n";
    }
#    $debug && print STDERR "$subn: \$TAXON=\n".Dumper($TAXON)."end of \$TAXON\n\n";
    for my $ti (sort {$a<=>$b} keys %{$TAXON->{taxon}}) {
#        $debug && print STDERR "$subn: \$ti=$ti: @{$TAXON->{taxon}->{$ti}}\n";
    }

    if ( !$REFSEQS->{refseq_loaded} ) {
        $REFSEQS = Annotate_Def::initRefseq();
        my $nstrain = 0;
        for my $fam (keys %{$REFSEQS->{refs}}) {
            $nstrain += scalar keys %{$REFSEQS->{refs}->{$fam}};
            $debug && print STDERR "$subn: \$fam=$fam\t\$nstrain=$nstrain\n";
        }
        my @fam = sort keys %{$REFSEQS->{refs}};
        printf STDERR ("$subn: loaded RefSeqs for $nstrain strains in %d families: \n%s\n", $#fam+1, join("\n", @fam));
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

        if ($DOWNLOAD_TAXON) {
            print STDERR "\n$subn: \$fam=$fam, performed live download for taxon\n";
        } else {
            print STDERR "\n$subn: \$fam=$fam, existing file used for taxon, no live download\n";
        }
        $outFileName = Annotate_Download::downloadTaxon( $fam, $exe_dir, $DOWNLOAD_TAXON);
         if (!$outFileName) {
            # Effectively try again to download if unsuccessful
            $outFileName = Annotate_Download::downloadTaxon( $fam, $exe_dir, $DOWNLOAD_TAXON);
        }
        print STDERR "$subn: \$outFileName=$outFileName\n";

        if (!$outFileName) {
            print STDERR "$subn: \$fam=$fam\n";
            print STDERR "$subn: ERROR: \$fam=$fam \$taxonFileName='$outFileName' is empty.\n";
            print STDERR "$subn: \$fam=$fam, perhaps turn on the download options in sub Annotate_Download::downloadTaxon\n\n";
            next;
        }

        my $newTax = loadXmlTaxon( $fam, $outFileName, $exe_dir);
#        $debug && print STDERR "$subn: \$newTax=".scalar(keys %$newTax)."\n".Dumper($newTax)."end of \$newTax\n\n";
        for my $ti (sort keys %{$newTax}) { $debug && print STDERR "$subn: \$ti=$ti: @{$newTax->{$ti}}\n"; }

        my $status = { total => 0, strain => 0, identical => 0, new => 0, };
        my $msg = '';
        my $str = sprintf("# Family=$fam taxon loaded from file=$outFileName \n");
        my $time = POSIX::strftime("%m/%d/%Y %H:%M:%S", localtime);
        $str .= sprintf("# Saved from Annotate_Dev.pm Ver$VERSION on %s \n", $time);
        my $str1 = '';
        my $ct = -1;
        for my $taxid (sort {$newTax->{$a}->[2]<=>$newTax->{$b}->[2] || $newTax->{$a}->[1]<=>$newTax->{$b}->[1] || $newTax->{$a}->[0]<=>$newTax->{$b}->[0]} keys %$newTax) {
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
#            $debug && print STDERR "$subn: \$ct=$ct \$taxid=$taxid \$newStr='$newStr'\n";

            my $oldStr = "";
            if (exists($TAXON->{taxon}->{$taxid})) {
                my $oldTax = $TAXON->{taxon}->{$taxid};
                $oldStr = sprintf("| %7d | %7d |", $oldTax->[0], $oldTax->[1]);
                $oldStr .= sprintf(" %7d | %7d |", $oldTax->[2], $oldTax->[3]);
                $oldStr .= sprintf(" %-28s |", $oldTax->[4]);
                $oldStr .= sprintf(" %-16s |", $oldTax->[5]);
                $oldStr .= sprintf(" %-15s |", $oldTax->[6]);
            }

            $status->{strain}++;
            if ($newStr eq $oldStr) {
                $status->{identical}++;
                $msg .= "#$ct\t---\t$newStr\n";
                $debug && print STDERR "$subn: \$ct=$ct \$taxid=$taxid identical Taxon info found. Skip...\n";
#                $debug && print STDERR "$subn: \$ct=$ct \$taxid=$taxid \$oldStr='$oldStr'\n";
#                $debug && print STDERR "$subn: \$ct=$ct \$taxid=$taxid \$newStr='$newStr'\n";
#                $debug && print STDERR "$subn: \$TAXON_UPDATED=$TAXON_UPDATED\n";
            } else {
                $TAXON_UPDATED = 1;
                $status->{new}++;
                $msg .= "#$ct\tnew\t$newStr\n";
                $debug && print STDERR "$subn: \$ct=$ct \$taxid=$taxid found new Taxon info from NCBI download for $fam:\n";
                $debug && print STDERR "$subn: \$ct=$ct \$taxid=$taxid \$oldStr='$oldStr'\n";
                $debug && print STDERR "$subn: \$ct=$ct \$taxid=$taxid \$newStr='$newStr'\n";
                $debug && print STDERR "$subn: \$TAXON_UPDATED=$TAXON_UPDATED\n";
            }
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
        print STDERR "$subn: Message: The new taxon file is saved in '$exe_dir/$taxonTempName'\n";
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
        my $cmd = "mv $exe_dir/$taxonTempName $exe_dir/$taxonFileName";
        $result = `$cmd`;
        print STDERR "$subn: Moved $taxonTempName to $taxonFileName\n";
        print STDERR "$subn: \$cmd='$cmd' \$result='$result'\n";
    }

    return $count;
} # sub checkAllTaxon


=head2
 sub checkAllRefseq look over each RefSeq (NC_??????) in $REFSEQS and see if there is any update in the genbank file
=cut

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

        if ($DOWNLOAD_REFSEQ) {
            print STDERR "\n$subn: \$fam=$fam, performed live download for RefSeq\n";
        } else {
            print STDERR "\n$subn: \$fam=$fam, existing file used for RefSeq, no live download\n";
        }
        $outFileName = Annotate_Download::downloadRefseq( $fam, $exe_dir, $DOWNLOAD_REFSEQ);
        if (!$outFileName) {
            print STDERR "$subn: \$fam=$fam\n";
            print STDERR "$subn: ERROR: \$fam=$fam \$outFileName='$outFileName' is empty.\n";
            print STDERR "$subn: \$fam=$fam Please turn on the download options in sub Annotate_Download::downloadRefseq\n\n";
            next;
        }
        print STDERR "\n$subn: \$fam=$fam \$outFileName=$outFileName\n";

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
            my $errcode = Annotate_Download::check1Refseq($new, $exe_dir);
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


=head2
 sub check1Refseq Bio::Seq object based on genbank file, checks if it has same feature annotation
 as anyfile in refseq directory. returs following values:
0: no difference, or existing gbk file
1. found difference in term of mat_peptide annotation, alert the user
=cut

sub check1Refseq {
    my ( $new, $exe_dir) = @_;

    my $debug = 1 && $debug_all;
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
                  # identical
                  $debug && print STDERR "$subn: \$newct=$newct \$newfeat='$newfeat' \$oldfeat='$oldfeat'\n";
                  $seen = 1;
                  last;
                } elsif ($newfeat1 eq $oldfeat1) {
                  # similar with only difference as '<' or '>'
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


1;
