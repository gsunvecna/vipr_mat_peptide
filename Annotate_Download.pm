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

#use version;
our $VERSION = qw('1.2.0'); # Apr 12 2013

use Annotate_Def;

my $debug_all = 0;

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

sub setDebugAll {
    my ($debug) = @_;
    $debug_all = $debug;
} # sub setDebugAll


=head2
 sub downloadRefseq searches genbank for all refseqs for a family, then downloads all the refseqs
=cut

sub downloadRefseq {
    my ($fam, $exe_dir, $download_REFSEQ) = @_;

    my $debug = 0 || $debug_all;
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
        my $bakFileName = "$outFileName.bak1";
        if (-e "$exe_dir/$outFileName") {
            my $copies = Annotate_Util::backupFiles( $exe_dir, $outFileName, 8);
            print STDERR "$subn: \$fam=$fam subroutine backupFiles made $copies backups for $outFileName\n";

        }
        # Download the NC_* as a set
        # retrieve data in batches of 100
        open (my $OUT, ">$exe_dir/$outFileName") || die "$subn: Can't open file'$outFileName'!\n";
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
              my $diff = `diff $exe_dir/$bakFileName $exe_dir/$outFileName 2>&1`;
              my $result = `ls -l $exe_dir/$outFileName*`;
              chomp $result;
              if ($diff) {
                print STDERR "$subn: ERROR: \$fam=$fam difference between $bakFileName and $outFileName\n";
                print STDERR "$subn: \$result=\n$result\n";
              } else {
                print STDERR "$subn: \$fam=$fam identical $bakFileName and $outFileName\n";
                print STDERR "$subn: \$result=\n$result\n";
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

    my $debug = 0 || $debug_all;
    my $subn = 'downloadTaxon';

    my $count = 0;
    $exe_dir = './' if (!$exe_dir);
    $debug && print STDERR "$subn: \$exe_dir='$exe_dir'\n";

    if ($download_TAXON) {
        print STDERR "\n$subn: \$fam=$fam, performing live download for taxon\n";
    } else {
        print STDERR "\n$subn: \$fam=$fam, existing file to be used, no live download requested\n";
    }
    # Search for NC_* for each family from genbank
    my $web = '';
    my $key = '';
    my $base = 'http://eutils.ncbi.nlm.nih.gov/entrez/eutils/';
    my $outFileName = '';
#    $outFileName = "RefSeq_${fam}.gbk";
#    print STDERR "\n$subn: \$fam=$fam \$outFileName=$outFileName\n";

    my $taxonFileName = "taxon_${fam}.xml";
    $download_TAXON = 0 if (!$download_TAXON);
    if ( !$download_TAXON ) {
        $debug && print STDERR "$subn: \$fam=$fam, \$download_TAXON=$download_TAXON file='$exe_dir/$taxonFileName'\n";
        if (!-e "$exe_dir/$taxonFileName") {
            print STDERR "$subn: ERROR: \$fam=$fam, file '$exe_dir/$taxonFileName' non-existant, while no live download is requested\n";
            $taxonFileName = '';
        } else {
            print STDERR "$subn: \$fam=$fam, file '$exe_dir/$taxonFileName' located\n";
        }
    } else {
        #assemble the esearch URL
        my $query = "${fam}[orgn]";
        print STDERR "$subn: \$fam=$fam \$outFileName=$outFileName \$query='$query'\n";
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
        my $tmpfile = $taxonFileName . ".0";
for (my $i=0; $i<3; $i++) {
        open (my $OUT, ">$exe_dir/$tmpfile") || die "$subn: Can't open file='$tmpfile'!\n";
        my $retmax = 5000;
        # download the taxon
        for (my $retstart = 0; $retstart < $count; $retstart += $retmax) {
            $debug && print STDERR "$subn: \$fam=$fam \$i=$i \$count=$count Downloading \$retstart=$retstart\n";
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
            for my $ntry (0 .. 4) {
                ($ntry>0) && sleep( $delay * $ntry );
                $efetch_out = get($efetch_url);
                $debug && print STDERR "$subn: \$ntry=$ntry \$efetch_out=".substr($efetch_out, 0, 300)."\n";
                last if ($efetch_out && $efetch_out !~ /Unable to obtain query/);
            }
            if (!$efetch_out || $efetch_out =~ /Unable to obtain query/) {
                $debug && print STDERR "$subn: ERROR: \$count=$count \$retstart=$retstart \$efetch_url='$efetch_url'\n";
                $debug && print STDERR "$subn: \$efetch_out=$efetch_out\n";
#                return undef ;
            } else {
                print $OUT "$efetch_out";
            }
        }
        close $OUT;
        $debug && print STDERR "$subn: \$i=$i Downloaded $count genomes for $fam, saved to $exe_dir/$tmpfile\n";

        last if (!-e "$exe_dir/$taxonFileName");
        my $s1 = -s "$exe_dir/$taxonFileName";
        my $s0 = -s "$exe_dir/$tmpfile";
        $debug && print STDERR "$subn: \$i=$i \$s0=$s0 for new file:$exe_dir/$tmpfile\n";
        $debug && print STDERR "$subn: \$i=$i \$s1=$s1 for old file:$exe_dir/$taxonFileName\n";
        last if ($s0>=$s1);
}

        # Move $tmpfile to $taxonFileName is there is update
if (!-z "$exe_dir/$tmpfile") {
        my $cmd = "diff $exe_dir/$tmpfile $exe_dir/$taxonFileName 2>&1";
        my $result = ` $cmd `;
        chomp $result;
        $debug && print STDERR "$subn: \$cmd=$cmd\n";
        $debug && print STDERR "$subn: \$result='\n$result'\n" if ($result);
        if ($result) {
            if (-e "$exe_dir/$taxonFileName") {
                my $copies = Annotate_Util::backupFiles( $exe_dir, $taxonFileName, 8);
                print STDERR "$subn: for taxon file, subroutine backupFiles made $copies backups\n";
            }
            $cmd = "mv $exe_dir/$tmpfile $exe_dir/$taxonFileName";
            $result = ` $cmd `;
            $debug && print STDERR "$subn: \$cmd=$cmd\n";
            $debug && print STDERR "$subn: \$result='\n$result'\n" if ($result);
            print STDERR "$subn: Downloaded $count genomes for $fam, saved to $exe_dir/$taxonFileName\n";
        } else {
            print STDERR "$subn: Downloaded $count genomes for $fam, no difference with stored taxon.\n";
            print STDERR "$subn: Downloaded $count genomes for $fam, new download saved in $tmpfile.\n";
        }
    }
}

    return $taxonFileName;
} # sub downloadTaxon


=head2
 sub loadXmlTaxon reads an XML file, searchs for each taxid, and its name, speciesid, species,
 genus, genusid, family, and familyid, to be compared with exsting data
=cut

sub loadXmlTaxon {
    my ($fam, $taxonFileName, $exe_dir) = @_;

    my $debug = 0 || $debug_all;
    my $subn = 'loadXmlTaxon';

    my $newTax = {};
    $exe_dir = './' if (!$exe_dir);
    $debug && print STDERR "$subn: \$exe_dir='$exe_dir'\n";
    if (!-e "$exe_dir/$taxonFileName") {
        print STDERR "$subn: ERROR: \$fam=$fam file '$exe_dir/$taxonFileName' non-existant, abort\n";
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
#            $debug && print STDERR "$subn: \$_='$_'\n";
            last if ($_ =~ /^<\/Taxon>/);
        }
#        $debug && print STDERR "$subn: \@set=".Dumper(@set)."\n";

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
#        $debug && print STDERR "$subn: \$tax=\n".Dumper($tax)."\n";

    }
    $debug && print STDERR "$subn: \$newTax=\n".Dumper($newTax)."\n";

    return $newTax;
} # sub loadXmlTaxon


=head2
 sub getTaxonString takes a hash of a taxon, and finds the relative info, and returns an array
=cut

sub getTaxonString {
    my ($tax) = @_;

    my $debug = 0 || $debug_all;
    my $subn = 'getTaxonString';

#    $debug && print STDERR "$subn: Just got in the sub: \$tax=\n".Dumper($tax)."\n";
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
#        $debug && print STDERR "$subn: \$i=$i \$line='$line'\n";
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
 sub taxToString takes a taxon in the Array form, and turn it into formated string
=cut

sub taxToString {
    my ( $tax) = @_;

    my $debug = 0 || $debug_all;
    my $subn = 'taxToString';

    my $str = '';
    $str .= sprintf("| %7d |", $tax->[0]);
    $str .= sprintf(" %7d |", $tax->[1]);
    $str .= sprintf(" %7d |", $tax->[2]);
    $str .= sprintf(" %7d |", $tax->[3]);
    $str .= sprintf(" %-28s |", $tax->[4]);
    $str .= sprintf(" %-16s |", $tax->[5]);
    $str .= sprintf(" %-15s |", $tax->[6]);

    return $str;
} # sub taxToString

=head2
 sub checkFamilyTaxon looks over one family and see if there is any update in its taxonomy from NCBI
=cut

sub checkFamilyTaxon {
    my ( $fam, $exe_dir, $outFileName) = @_;

    my $debug = 0 || $debug_all;
    my $subn = 'checkFamilyTaxon';

    my $count = { total => 0 };
    $exe_dir = './' if (!$exe_dir);
    $debug && print STDERR "$subn: \$exe_dir='$exe_dir'\n";

    # Check strains for each family from genbank
    my $status = {
                   TAXON_UPDATED => 0,
                   file  => $outFileName,
                   old   => 0,
                   total => 0,
                   strain => 0,
                   identical => 0,
                   new   => 0,
                 };
    if (!$outFileName) {
        $debug && print STDERR "$subn: \$fam=$fam\n";
        print STDERR "$subn: ERROR: \$fam=$fam \$taxonFileName='$outFileName' can't be found\n";
        print STDERR "$subn: \$fam=$fam, please check such taxon file\n\n";
        return $status;
    }

    my $newTax = loadXmlTaxon( $fam, $outFileName, $exe_dir);
#    $debug && print STDERR "$subn: \$newTax=".scalar(keys %$newTax)."\n".Dumper($newTax)."end of \$newTax\n\n";
    print STDERR "$subn: loaded ".scalar(keys %$newTax)." taxons from $outFileName\n";
    for my $ti (sort keys %{$newTax}) { $debug && print STDERR "$subn: \$ti=$ti: @{$newTax->{$ti}}\n"; }

    for my $k (sort keys %{$TAXON->{taxon}}) {
            ++$status->{old} if ($fam eq $TAXON->{taxon}->{$k}->[6]);
    }
    my $time = POSIX::strftime("%m/%d/%Y %H:%M:%S", localtime);
    my $str = sprintf("# Saved from Annotate_Download.pm Ver$VERSION on %s \n", $time);
    $status->{comment} = $str;
    my $msg = '';
    my $str1 = '';
    my $ct = -1;
    my $cskipped = 0;
    for my $taxid (keys %$newTax) {
        if ($newTax->{$taxid}->[1]<0) {
            $cskipped++;
            print STDERR "$subn: skip #$cskipped:'@{$newTax->{$taxid}}'.\n";
            next;
        }
        $ct++;
        $status->{total}++;
        my $newStr = '';
        $newStr = Annotate_Download::taxToString($newTax->{$taxid});
        $str .= $newStr . "\n";
#        $debug && print STDERR "$subn: \$ct=$ct \$taxid=$taxid \$newStr='$newStr'\n";

        my $oldStr = "";
        if (exists($TAXON->{taxon}->{$taxid})) {
            $oldStr = Annotate_Download::taxToString( $TAXON->{taxon}->{$taxid} );
        }

        $status->{strain}++;
        if ($oldStr && ($newStr eq $oldStr)) {
            $status->{identical}++;
            $msg .= "#$ct\t---\t$newStr\n";
            $TAXON->{newtax}->{$taxid} = $newTax->{$taxid};
            $debug && print STDERR "$subn: \$ct=$ct \$taxid=$taxid identical Taxon info found. Skip...\n";
#            $debug && print STDERR "$subn: \$ct=$ct \$taxid=$taxid \$oldStr='$oldStr'\n";
#            $debug && print STDERR "$subn: \$ct=$ct \$taxid=$taxid \$newStr='$newStr'\n";
#            $debug && print STDERR "$subn: \$TAXON_UPDATED=$TAXON_UPDATED\n";
        } else {
            $status->{TAXON_UPDATED} = 1;
            $status->{new}++;
            $msg .= "#$ct\tnew\t$newStr\n";
#            $debug && print STDERR "$subn: \$newTax=".Dumper($newTax)."\n";
            # add new taxon to the list of all taxons
            $TAXON->{newtax}->{$taxid} = $newTax->{$taxid};
            $debug && print STDERR "$subn: \$ct=$ct \$taxid=$taxid found new Taxon info from NCBI download for $fam:\n";
            print STDERR "$subn: \$ct=$ct \$oldStr='$oldStr'\n";
            print STDERR "$subn: \$ct=$ct \$newStr='$newStr'\n";
            $debug && print STDERR "$subn: \$status->{TAXON_UPDATED}=$status->{TAXON_UPDATED}\n";
        }
    } # for my $taxid (sort {}) {

    my $msg1 = "# $fam found taxon from NCBI: $status->{total} total:";
    $msg1 .= " $status->{strain} strains,";
    $msg1 .= " $status->{identical} identical,";
    $msg1 .= " $status->{new} new";
    $status->{comment} .= "# Family=$fam taxon loaded from file=$outFileName, total=$status->{total}\n";
    $status->{comment} .= $msg1;
    # Save to Anntate_taxon_records.txt
    print STDERR "$subn: \$status=\n$status->{comment}\n";
    print STDERR "$subn: \$msg=\n$msg\n";

    return $status;
} # sub checkFamilyTaxon

=head2
 sub checkAllTaxon looks over each family and see if there is any update in its taxonomy from NCBI
=cut

sub checkAllTaxon {
    my ( $exe_dir, $DOWNLOAD_TAXON, $update) = @_;

    my $debug = 0 || $debug_all;
    my $subn = 'checkAllTaxon';

    # Load the $TAXON if not already
    if ( !$TAXON->{taxon_loaded} ) {
        $TAXON = Annotate_Def::loadTaxonTable( $exe_dir);
        my $ctLines = scalar(keys %{$TAXON->{taxon}});
        $debug && print STDERR "$subn: loaded $ctLines lines from file '$exe_dir/$TAXON->{taxon_fn}'\n";
    }
#    $debug && print STDERR "$subn: \$TAXON=\n".Dumper($TAXON)."end of \$TAXON\n\n";
    for my $ti (sort {$a<=>$b} keys %{$TAXON->{taxon}}) {
#        $debug && print STDERR "$subn: \$ti=$ti: @{$TAXON->{taxon}->{$ti}}\n";
    }

    # Load the $REFSEQS if not already
    if ( !$REFSEQS->{refseq_loaded} ) {
        $REFSEQS = Annotate_Def::initRefseq();
        my $nrefseq = 0;
        for my $fam (sort keys %{$REFSEQS->{refs}}) {
            my $nstrain = scalar keys %{$REFSEQS->{refs}->{$fam}};
            $nrefseq += $nstrain;
            $debug && printf STDERR ("$subn: \$fam=%15s\tRefSeq=$nstrain\n", $fam);
        }
        my @fam = sort keys %{$REFSEQS->{refs}};
#        printf STDERR ("$subn: loaded RefSeqs for $nstrain strains in %d families: \n%s\n", $#fam+1, join("\n", @fam));
        printf STDERR ("$subn: loaded RefSeqs for $nrefseq strains in %d families: \n", $#fam+1);
        for my $i(0..$#fam) { printf STDERR ("%2d\t%s\n", $i+1, $fam[$i]); }
    }
    $debug && print STDERR "$subn: \$REFSEQS=\n".Dumper($REFSEQS)."end of \$REFSEQS\n\n";

    my $count = { total => 0 };
    $exe_dir = './' if (!$exe_dir);
    $debug && print STDERR "$subn: \$exe_dir='$exe_dir'\n";

    # Check strains for each family from genbank
#    my $DOWNLOAD_TAXON = 0;
    my $genera = {};
    my $TAXON_UPDATED = 0;
    for my $fam (sort keys %{$REFSEQS->{refs}}){
#        next if ($fam ne 'Caliciviridae');

        my $outFileName = '';
        $outFileName = Annotate_Download::downloadTaxon( $fam, $exe_dir, $DOWNLOAD_TAXON);
        if (!$outFileName && $DOWNLOAD_TAXON) {
            # Effectively try again to download if unsuccessful
            $outFileName = Annotate_Download::downloadTaxon( $fam, $exe_dir, $DOWNLOAD_TAXON);
        }
        $debug && print STDERR "$subn: \$outFileName=$outFileName after calling downloadTaxon\n";

        # Take the downloaded file, read the taxon in it and compare with existing taxon data
        my $status = Annotate_Download::checkFamilyTaxon( $fam, $exe_dir, $outFileName);
        $TAXON_UPDATED = $TAXON_UPDATED || $status->{TAXON_UPDATED};
        $debug && print STDERR "$subn: \$fam=$fam \$status=\n".Dumper($status)."End of \$status\n";

        $count->{$fam} = $status;
        $count->{total} += $status->{total};

    } # for my $fam (sort keys %{$REFSEQS->{refs}}){

    # If there is any update, write all data to file, backup older file
    my @fams = ( );
    for my $f (sort keys %{$count}) { next if ($f !~ m/viridae$/i); push @fams, $f; }
    my $taxonFileName = "Annotate_taxon_records.txt"; # plain text file storing the taxon info for all families
    my $taxonTempName = "$taxonFileName.0";

    my $taxonFile = undef;
    open($taxonFile, ">$exe_dir/$taxonTempName") || die "$subn: Can't open file='$taxonTempName'!\n";

    for my $fam (@fams ) {
        my $newTax = {};
        # gather all taxon for a family
        for my $taxid (keys %{$TAXON->{newtax}}) {
            $newTax->{$taxid} = $TAXON->{newtax}->{$taxid} if ($fam eq $TAXON->{newtax}->{$taxid}->[6]);
        }
        my $oldTax = {};
        # gather all taxon for a family
        for my $taxid (keys %{$TAXON->{taxon}}) {
            $oldTax->{$taxid} = $TAXON->{taxon}->{$taxid} if ($fam eq $TAXON->{taxon}->{$taxid}->[6]);
        }
        my $oldLen = scalar(keys %$oldTax);
        my $newLen = scalar(keys %$newTax);
        $debug && print STDERR "$subn: \$fam=$fam \$oldTax=$oldLen \$newTax=$newLen\n";
        my $diff = {};
        for my $k (keys %$oldTax) {
            $diff->{$k} = $oldTax->{$k} if (!exists($newTax->{$k}));
        }
        print STDERR "$subn: \$fam=$fam \$oldTax=$oldLen \$newTax=$newLen \$diff=".scalar(keys%$diff)."\n";
        $debug && print STDERR "$subn: \$fam=$fam \$diff=\n".Dumper($diff)."End of \$diff\n";
        my @ti = (keys %$diff);
        if ($#ti>=0) {
            print STDERR "$subn: \$fam=$fam following are the $#ti missing taxons in new update:\n";
            for my $i (0 .. $#ti) {
                print STDERR "$subn: #$i: ".Annotate_Download::taxToString( $diff->{$ti[$i]})."\n";
            }
        }
        my $status = $count->{$fam};
        $debug && print STDERR "$subn: \$fam=$fam \$status=\n".Dumper($status)."End of \$status\n";
        my $msg = '';
        my $outFileName = $status->{file};
        my $str = sprintf("# Family=$fam taxon loaded from file=$outFileName \n");
        my $time = POSIX::strftime("%m/%d/%Y %H:%M:%S", localtime);
        $str .= sprintf("# Saved from Annotate_Download.pm Ver$VERSION on %s \n", $time);
        my $str1 = '';
        my $ct = -1;
        for my $taxid (sort {$newTax->{$a}->[2]<=>$newTax->{$b}->[2] || $newTax->{$a}->[1]<=>$newTax->{$b}->[1] || $newTax->{$a}->[0]<=>$newTax->{$b}->[0]} keys %$newTax) {
            next if ($newTax->{$taxid}->[1]<0);
            $ct++;
            my $newStr = Annotate_Download::taxToString( $newTax->{$taxid});
            $str .= $newStr . "\n";
            $msg .= "#$ct\tnew\t$newStr\n";
        } # for my $taxid (sort {}) {
        $debug && print STDERR "$subn: \$fam=$fam \$genera=\n".Dumper($genera)."End of \$genera\n";

        my $msg1 = "$fam found taxon from NCBI: $status->{total} total:";
        $msg1 .= " $status->{strain} strains,";
        $msg1 .= " $status->{identical} identical,";
        $msg1 .= " $status->{new} new\n";
        $msg = $msg1 . $msg;
        # Save to Anntate_taxon_records.txt
        $str .= sprintf("# Family=$fam taxon loaded from file=$outFileName, total=%d\n", $status->{total});
        print $taxonFile "$str\n" || die "$subn: Can't write to file='$taxonFileName'!\n";
        print STDERR "$subn: \$str='\n$str'\n";

    } # for my $fam (sort keys %{$REFSEQS->{refs}}){
    printf $taxonFile ("# Updated $count->{total} strains in %d families: @fams\n", scalar(@fams));
    close $taxonFile || die "$subn: Can't close file='$taxonFileName'!\n";

    # save the result to temp file, and report the result
#    printf STDERR ("$subn: # Checked $count strains in %d families: \n%s\n\n", scalar(@fams), join("\n", @fams));
    if ($DOWNLOAD_TAXON) {
        print STDERR "$subn: performed live download\n";
    } else {
        print STDERR "$subn: existing file used, no live download\n";
    }
    printf STDERR ("$subn: # Checked $count->{total} strains in %d families: \n", $#fams+1);
    for my $i (0 .. $#fams) {
        my $fam = $fams[$i];
        printf STDERR ("%d\t", $i+1);
        printf STDERR "$fam\t";
        printf STDERR ("old=%5d\t", $count->{$fam}->{old});
        printf STDERR ("total=%5d\t", $count->{$fam}->{total});
        printf STDERR ("strain=%5d\t", $count->{$fam}->{strain});
        printf STDERR ("identical=%5d\t", $count->{$fam}->{identical});
        printf STDERR ("new=%3d\n", $count->{$fam}->{new});
    }
    printf STDERR ("\n");

    # if needed, backup files, move new data to designated file
    my $cmd = '';
    my $bak1 = "$taxonFileName.bak1";
    my $TAXON_SAVED = 0;
    if ( !$TAXON_UPDATED ) {
        print STDERR "$subn: Message: No update is found between taxon from NCBI and stored in '$exe_dir/$taxonFileName'\n";
        print STDERR "$subn: Message: No change is made to the current file '$exe_dir/$taxonFileName'\n";
        print STDERR "$subn: Message: The new taxon file is saved in '$exe_dir/$taxonTempName'\n";
    } else {
        print STDERR "$subn: Found update in taxon from NCBI, as compared to file='$exe_dir/$taxonFileName'\n";
        my $s1 = -s "$exe_dir/$taxonFileName";
        my $s0 = -s "$exe_dir/$taxonTempName";
        if ($s0<$s1) {
            print STDERR "$subn: Found update in taxon from NCBI, however the new file is smaller than existing file\n";
            print STDERR "$subn: This probably indicates problems during download. Please try again.\n";
        } elsif ($update) {
            my $result = '';
            print STDERR "$subn: Found update in taxon from NCBI, full info will be saved to file='$exe_dir/$taxonFileName'\n";
            my $copies = Annotate_Util::backupFiles( $exe_dir, $taxonFileName, 8);
            print STDERR "$subn: for taxon file='$taxonFileName', subroutine backupFiles made $copies backups\n";

            $TAXON_SAVED = 1;
            $cmd = "mv $exe_dir/$taxonTempName $exe_dir/$taxonFileName";
            $result = `$cmd`;
            print STDERR "$subn: Moved $taxonTempName to $taxonFileName\n";
            print STDERR "$subn: \$cmd='$cmd' \$result='\n$result'\n";
        } else {
            print STDERR "$subn: Found update in taxon from NCBI\n";
            print STDERR "$subn: However, file='$exe_dir/$taxonFileName' won't be updated because \$update='$update'\n";
            print STDERR "$subn: In order to update file='$exe_dir/$taxonFileName', set the -update option in commandline\n";
        }
    }
    $cmd = "ls -l $exe_dir/$taxonFileName $exe_dir/$bak1";
    $cmd .= " $exe_dir/$taxonTempName " if (!$TAXON_SAVED);
    print STDERR ` $cmd `;
    $cmd = "ls -l $exe_dir/taxon_*xml";
    print STDERR ` $cmd `;

    return $count;
} # sub checkAllTaxon


=head2
 sub checkAllRefseq look over each RefSeq (NC_??????) in $REFSEQS and see if there is any update in the genbank file
 $DOWNLOAD_REFSEQ controls if live download will be performed
 $update controls if the updated gene symbols will be writen to the text file.
=cut

sub checkAllRefseq {
    my ( $exe_dir, $DOWNLOAD_REFSEQ, $update) = @_;

    my $debug = 0 || $debug_all;
    my $subn = 'checkAllRefseq';

    # Load the $REFSEQS if not already
    if ( !$REFSEQS->{refseq_loaded} ) {
        $REFSEQS = Annotate_Def::initRefseq();
        my $nstrain = 0;
        my @fam = sort keys %{$REFSEQS->{refs}};
        for my $fam (@fam) { $nstrain += scalar keys %{$REFSEQS->{refs}->{$fam}}; }
        printf STDERR ("$subn: loaded RefSeqs for $nstrain strains in %d families: \n", $#fam+1);
        printf STDERR ("%s\t%s\t%s\n", '#', 'allowed', 'Family');
        for my $i(0..$#fam) { printf STDERR ("%2d\t%3d\t%s\n", $i+1, scalar keys %{$REFSEQS->{refs}->{$fam[$i]}}, $fam[$i]); }
    }
    $debug && print STDERR "$subn: \$REFSEQS=\n".Dumper($REFSEQS)."end of \$REFSEQS\n\n";

    my $count = {total=>0};
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
#    my $DOWNLOAD_REFSEQ = 0;
    # Search for NC_* for each family from genbank
    my $meta = 1; # To load the content of Annotate_symbol_records.txt
    Annotate_Def::load_gene_symbol( $exe_dir, $meta);
    $debug && print STDERR "$subn: \$meta='$meta'\n\n";

    my $refseqs_excluded = {
                             'Bunyaviridae' => {
                                                   'NC_004158' => 'Having only 1 mat_peptide',
                                                 },
                             'Caliciviridae' => {
                                                   'NC_007916' => 'Having only 2 mat_peptides',
                                                   'NC_006554' => 'Having only 1 mat_peptide',
                                                   'NC_010624' => 'Having only 1 mat_peptide',
                                                 },
                             'Flaviviridae' => {
                                                   'NC_009942' => 'Some minor problems as of 5/14/2014',
                                                   'NC_009824' => 'E2/NS1 combined at 1489..2544',
                                                   'NC_012671' => 'Having only 1 mat_peptide',
                                                 },
                             'Hepeviridae' => {
                                                   'NC_015521' => 'huge gape between the 4 mat_peptides',
                                                 },
                             'Picornaviridae' => {
                                                   'NC_003983' => 'double coverage between 3A-3B, 3B-3C',
                                                   'NC_013695' => 'species=1330521 has 2 refseqs, using NC_010415 instead',
                                                 },
                             'Reoviridae' => {
                                                   'NC_003758' => 'Ask Richard for opinion on Reoviridae',
                                                   'NC_004278' => 'Ask Richard for opinion on Reoviridae',
                                                 },
                             'Rhabdoviridae' => {
                                                   'NC_001560' => 'Rhabdoviridae doesn\'t have mat_peptide',
                                                 },
                             'Togaviridae' => {
                                                   'NC_016959' => 'Having only 1 mat_peptide',
                                                   'NC_016960' => 'Having only 1 mat_peptide',
                                                   'NC_016961' => 'Having only 1 mat_peptide',
                                                 },
                           };

    my @fams = sort keys %{$REFSEQS->{refs}};
    for my $fam (@fams){
        my $web = '';
        my $key = '';
        my $ct = 0;
        my $base = 'http://eutils.ncbi.nlm.nih.gov/entrez/eutils/';
        my $outFileName = '';
        my $status = {
                       total => 0,
                       old   => 0,
                       no_file_on_disk => 0,
                       has_matpeptide  => 0,
                       change_in_matpeptide => 0,
                       allowed  => 0,
                       newSymbols  => 0,
                       updated     => 0,
                     };
        $status->{old} = scalar(keys %{$REFSEQS->{refs}->{$fam}});
        $count->{$fam} = $status;
#        $outFileName = "RefSeq_${fam}.gbk";

        my $fam_tested;
#        $fam_tested = 'Flaviviridae';
        if ($fam_tested) {
            next if ($fam ne $fam_tested);
        }

        if ($DOWNLOAD_REFSEQ) {
            print STDERR "\n$subn: \$fam=$fam, performed live download for RefSeq\n";
        } else {
            print STDERR "\n$subn: \$fam=$fam, existing file used for RefSeq, no live download\n";
        }
        $outFileName = Annotate_Download::downloadRefseq( $fam, $exe_dir, $DOWNLOAD_REFSEQ);
        if ($DOWNLOAD_REFSEQ && !$outFileName) {
            # Effectively try again to download if unsuccessful
            $outFileName = Annotate_Download::downloadRefseq( $fam, $exe_dir, $DOWNLOAD_REFSEQ);
        }
        if (!$outFileName) {
            $debug && print STDERR "$subn: \$fam=$fam\n";
            print STDERR "$subn: ERROR: \$fam=$fam \$outFileName='$outFileName' is empty, while \$DOWNLOAD_REFSEQ=$DOWNLOAD_REFSEQ\n";
            print STDERR "$subn: \$fam=$fam Please turn on the download options in the calling program\n\n";
            next;
        }
        print STDERR "\n$subn: \$fam=$fam \$outFileName=$outFileName after calling Annotate_Download::downloadRefseq\n";

        # Check each genbank file, if it has mat_peptide, compare with stored file, alert if different
        my $seqio = Bio::SeqIO->new( -file => "$exe_dir/$outFileName");
        my $nseq = 0;
        my $result = {}; # one line for each RefSeq
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
#            $debug && print STDERR "$subn: \$nseq=$nseq \$acc=$acc\n";
            my $errcode = Annotate_Download::check1Refseq( $new, $exe_dir, $update );
            $debug && print STDERR "$subn: \$nseq=$nseq \$acc=$acc \$errcode=\n".Dumper($errcode)."\n";

#            print STDERR "$subn: \$nseq=$nseq \$acc=$acc\n";
            if (exists($refseqs_excluded->{$fam}) && 
                  exists($refseqs_excluded->{$fam}->{$acc}) && 
                  $refseqs_excluded->{$fam}->{$acc}) {
                print STDERR "$subn: \$nseq=$nseq \$acc=$acc of family=$fam has been excluded: '$refseqs_excluded->{$fam}->{$acc}'\n";
            }
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
            $status->{allowed}++ if ($errcode->{allowed});
            $status->{updated}++ if ($errcode->{updated});
            $status->{newSymbols} += scalar(@{$errcode->{newSymbols}}) if (scalar(@{$errcode->{newSymbols}})>=0);
            my $msg1 = "#$nseq\t";
            $msg1 .= "$fam\t";
            $msg1 .= "$acc\t";
            $msg1 .= ($errcode->{has_matpeptide}) ? "$errcode->{has_matpeptide}\t" : "-\t";
            $msg1 .= ($errcode->{no_file_on_disk}) ? "$errcode->{no_file_on_disk}\t" : "-\t";
            $msg1 .= ($errcode->{change_in_matpeptide}) ? "$errcode->{change_in_matpeptide}\t" : "-\t";
            $msg1 .= ($errcode->{allowed}) ? "$errcode->{allowed}\t" : "-\t";
            $msg1 .= ($errcode->{updated}) ? "$errcode->{updated}\t" : "-\t";
            $msg1 .= (scalar(@{$errcode->{newSymbols}})>=0) ? sprintf("%d\t",scalar(@{$errcode->{newSymbols}})) : "-\t";
            $msg1 .= "$taxid\t";
            $msg1 .= "$speciesid\t";
            $msg1 .= "'$species'\n";

            $result->{$speciesid}->{$acc} = $msg1;
        }

        # collect all summaries
        my $msg = "#\tfamily  \taccession\t#mat\tno_file\tchanged\tallowed\tupdated\tnewSyms\ttaxid\tspeciesid\tspecies\n";
        for my $speciesid (sort {$a<=>$b} keys %{$result}) {
          for my $acc (sort keys %{$result->{$speciesid}}) {
            $msg .= $result->{$speciesid}->{$acc};
          }
        }
        $msg .= "$subn: ERROR: total RefSeq $status->{total} doesn't match \$ct=$ct" if ($ct && $ct!=$status->{total});
        $count->{total} += $status->{total};

        $msg .= "family: $fam\t";
        $msg .= sprintf("total=%3d\t", $status->{total});
        $msg .= ":$status->{has_matpeptide}\t";
        $msg .= ":$status->{no_file_on_disk}\t";
        $msg .= ":$status->{change_in_matpeptide}\t";
        $msg .= ":$status->{allowed}\t";
        $msg .= ":$status->{updated}\t";
        $msg .= ":$status->{newSymbols}\n";

        printf STDERR ("$msg\n");
        printf STDERR ("$subn: family=$fam total RefSeqs checked:  $status->{total}\n");
        printf STDERR ("$subn: family=$fam Refseqs w/ mat_peptide: $status->{has_matpeptide}/$status->{total}\n");
        printf STDERR ("$subn: family=$fam Refseqs w/o gbk file:   $status->{no_file_on_disk}/$status->{total}\n");
        printf STDERR ("$subn: family=$fam Refseqs w/ changes:     $status->{change_in_matpeptide}/$status->{total}\n");
        printf STDERR ("$subn: family=$fam Refseqs allowed:        $status->{allowed}/$status->{total}\n");
        printf STDERR ("$subn: family=$fam Refseqs updated:        $status->{updated}/$status->{total}\n\n");

    }

    $debug && print STDERR "$subn: \$count=\n".Dumper($count)."\n";
    my $newSymbols = 0;
    if ($DOWNLOAD_REFSEQ) {
        print STDERR "$subn: performed live download\n";
    } else {
        print STDERR "$subn: existing file used, no live download\n";
    }
    printf STDERR ("$subn: # Checked $count->{total} strains in %d families: \n", scalar(@fams));
    for my $i (0 .. $#fams) {
        my $fam = $fams[$i];
        printf STDERR ("%2d  ", $i+1);
        printf STDERR ("%16s\t", $fam);
        printf STDERR ("old=%3d\t", $count->{$fam}->{old});
        printf STDERR ("total=%3d  ", $count->{$fam}->{total});
        printf STDERR ("has_mat=%3d  ", $count->{$fam}->{has_matpeptide});
        printf STDERR ("no_file=%3d  ", $count->{$fam}->{no_file_on_disk});
        printf STDERR ("changed=%3d  ", $count->{$fam}->{change_in_matpeptide});
        printf STDERR ("newSymbols=%3d  ", $count->{$fam}->{newSymbols});
        printf STDERR ("allowed=%3d  ", $count->{$fam}->{allowed});
        printf STDERR ("updated=%3d\n", $count->{$fam}->{updated});
        $newSymbols += $count->{$fam}->{newSymbols};
    }
    printf STDERR ("\n");

    # Need to save any update
    if ($newSymbols) {
        $debug && printf STDERR ("$subn: There are $newSymbols new symbols to save for ".scalar(@fams)." families\n");
        my $t = Annotate_Def::saveSymbolText();
        $debug && printf STDERR ("$subn: \$t=\n$t\n");
        if ($update) {
            my $symbolFile = "Annotate_symbol_records.txt";
            my $tempFile = "$symbolFile.0";
            open (my $OUT, ">$exe_dir/$tempFile") || die "$subn: Can't open file '$tempFile'!\n";
            print $OUT "$t";
            close $OUT || die "$subn: Can't close file '$tempFile'!\n";
            print STDERR "$subn: new gene symbol file saved to '$exe_dir/$tempFile'\n";
            my $copies = Annotate_Util::backupFiles( $exe_dir, $symbolFile, 8);
            print STDERR "$subn: subroutine backupFiles made $copies backups of '$exe_dir/$symbolFile'\n";

            $copies = `mv $exe_dir/$tempFile $exe_dir/$symbolFile`;
            $debug && print STDERR "$subn: moved new gene symbol file to '$exe_dir/$symbolFile' from '$tempFile'\n";

            # open output file for writing
        } else {
            $debug && printf STDERR ("$subn: \$update=$update, no saving of gene symbol is done\n");
        }
    } else {
        $debug && printf STDERR ("$subn: There is no new symbols to save for ".scalar(@fams)." families\n");
    }

    return $count->{total};
} # sub checkAllRefseq


=head2
 sub check1Refseq Bio::Seq object based on genbank file, checks if it has same feature annotation
 as anyfile in refseq directory. returs following values:
0: no difference, or existing gbk file
1. found difference in term of mat_peptide annotation, alert the user
=cut

sub check1Refseq {
    my ( $new, $exe_dir, $update ) = @_;

    my $debug = 0 || $debug_all;
    my $subn = 'check1Refseq';

#    $debug = 1 if ($acc eq 'NC_013695');
    $exe_dir = './' if (!$exe_dir);
    my $errcode = {
                    no_file_on_disk => 0, # Set to 1 if there is no such genbank file
                    refseq_name     => '', # NC_001786
                    has_matpeptide  => 0,
                    manual_curation => 0,
                    change_in_matpeptide => 0,
                    change_in_file  => 0,
                    allowed         => 0, # is listed in the approved list
                    mat_peptides    => [], # 'NC_007916|CDS=GI:90403549|Loc=1063..1920',
                    new_mat_peptides=> [], # 'NC_007916|CDS=GI:90403549|Loc=1063..1920',
                    newSymbols      => [],
                    updated         => 0,
                  };
    my $acc = $new->accession_number;

    #check if the new RefSeq is among those approved to be used as standard
    my $taxid = $new->species->ncbi_taxid;
    my $speciesid;
    $speciesid = Annotate_Def::getTaxonInfo( $taxid);
    $debug && print STDERR "$subn: $acc \$speciesid='@$speciesid'\n";
    my $fam = $speciesid->[6];
    my $strain = $speciesid->[0];
    $debug && print STDERR "$subn: $acc \$fam='$fam' \$strain='$strain'\n";
    if (exists($REFSEQS->{refs}) 
        && exists($REFSEQS->{refs}->{$fam}) 
        && exists($REFSEQS->{refs}->{$fam}->{$strain}) 
        && ($acc eq $REFSEQS->{refs}->{$fam}->{$strain})
       ) {
        $errcode->{allowed} = 1;
        $debug && print STDERR "$subn: $acc \$errcode=".Dumper($errcode)."\n";
    }

    # Gather all mat_peptides
    my $newfeats = [ $new->get_SeqFeatures ];
    $debug && print STDERR "$subn: $acc has $#{$newfeats} features in total\n";

    my $has_matpeptide = 0;
    for (my $newct = 0; $newct<=$#{$newfeats}; $newct++) {
        my $newfeat = $newfeats->[$newct];
        my $newfeat_tag = $newfeat->primary_tag;
        $debug && print STDERR "$subn: $acc $newct: $newfeat_tag: ".$newfeat->location->to_FTstring."\n";
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
            my $newfeat_tag = $newfeat->primary_tag;
            $debug && print STDERR "$subn: $acc $newct: $newfeat_tag: ".$newfeat->location->to_FTstring."\n";
            next if ($newfeat->primary_tag ne 'mat_peptide');
            $errcode->{has_matpeptide} += 1;
            my $str = $acc."|$cdsid|Loc=".$newfeat->location->to_FTstring;
            push @{$errcode->{mat_peptides}}, $str;
        }
    }

    $debug && print STDERR "$subn: $acc \$errcode=".Dumper($errcode)."\n";
    ($errcode->{has_matpeptide}==0) && return $errcode;

    # Read in existing genbank file
    my $old = undef;
    my $oldfeats = [ ];
    # See if genbank file exists on disk
    my ($refacc, $specid, $spec) = Annotate_Def::get_refseq_acc( $new, $exe_dir);
    my $refseq_fn = '';
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
                $refseq_fn = $refacc.'.gb' if ($refacc);
            }
    $errcode->{refseq_name} = $refseq_fn;
    $debug && print STDERR "$subn: \$acc=$acc \$errcode->{refseq_name}='$errcode->{refseq_name}'\n";
#    if ($refseq_fn && $errcode->{refseq_name} !~ m/$acc/i) {
#        print STDERR "$subn: \$acc=$acc is of specie=$specid \$errcode->{refseq_name}='$errcode->{refseq_name}' don't match. Problem.\n";
#    }

    if ($refseq_fn !~ m/^$acc[.]gb$/i) {
        $errcode->{manual_curation} = 1;
#        return $errcode;
    }

    if (!-e "$exe_dir/refseq/$errcode->{refseq_name}") { $errcode->{no_file_on_disk} = 1;
    } elsif (!$errcode->{refseq_name}) {
        # empty strain
        $debug && print STDERR "$subn: \$acc=$acc \$errcode->{refseq_name}='$errcode->{refseq_name}'\n";
    } else {
        # Check between new and old RefSeq
        $old = Bio::SeqIO->new( -file => "$exe_dir/refseq/$errcode->{refseq_name}")->next_seq();
        $oldfeats = [ $old->get_SeqFeatures ];
    }
    # Check against stored RefSeq
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
                  $debug && print STDERR "$subn: \$newct=$newct identical \$newfeat='$newfeat' \$oldfeat='$oldfeat'\n";
                  $seen = 1;
                  last;
                } elsif ($newfeat1 eq $oldfeat1) {
                  # similar with only difference as '<' or '>'
                  $debug && print STDERR "$subn: \$newct=$newct similar \$newfeat='$newfeat' \$oldfeat='$oldfeat'\n";
                }
            }
        }
        $debug && print STDERR "$subn: \$newct=$newct \$seen=$seen \$newfeat='$newfeat'\n";
        if (!$seen) {
            $errcode->{change_in_matpeptide} += 1;
            push @{$errcode->{new_mat_peptides}}, $newfeat;
        }
    }
    $debug && print STDERR "$subn: \$errcode=\n".Dumper($errcode)."End of \$errcode\n";

    # Get the gene symbol for all the mat_peptides if there is any update in term of mat_peptide
    # Need to update the gene symbol list, and put the genbank file in the refseq directory
    if ($errcode->{change_in_matpeptide}>0) {
      if ($errcode->{allowed}) {
        Annotate_Download::updateRefseq( $exe_dir, $acc, $update, $errcode );
      } else {
        $debug && print STDERR "$subn: Found change for $acc, but it's not approved\n";
      }
    }
    return $errcode;
} # sub check1Refseq


=head2
 sub updateRefseq does following for the accession:
1. Get the genbank file from NCBI, save to the refseq directory and backup existing file
2. Compare the mat_peptides in the genbank file with eixsting list
3. For any new or changed mat_peptide, update the gene symbol, by changing the meta and meta1 fields in $GENE_SYM
=cut

sub updateRefseq {
    my ( $exe_dir, $acc, $update, $errcode) = @_;

    my $debug = 0 || $debug_all;
    my $subn = 'updateRefseq';

    my $base = 'http://eutils.ncbi.nlm.nih.gov/entrez/eutils/';
#    my $outFileName = undef;
    my $outFileName = '';

    # Check if the refseq is standard, ie. without any special treatment, such as in NC_004162 for Chikungunya virus, and NC_006558
    $debug && print STDERR "$subn: Refseq for \$acc=$acc: refseq_name=$errcode->{refseq_name}\n";
    if ($errcode->{manual_curation}==1 ) {
        print STDERR "$subn: manual_curation needed for \$acc=$acc: refseq_name=$errcode->{refseq_name}. Skip update.\n";
        return;
    }

    # First, download the genbank file
    # modeled after http://www.ncbi.nlm.nih.gov/books/NBK25498/#chapter3.Application_3_Retrieving_large -11/28/2012
    # Have to get the summary before downloading gbk files
    $outFileName = "refseq/$errcode->{refseq_name}";
    my $outTemp = "$outFileName.0";
    if ( $update ) {
        #append [accn] field to each accession
        #join the accessions with OR
    #    my $query = join('+OR+',@acc_array);
        my $query = $acc."[accn]";

        #assemble the esearch URL
        my $url = $base . "esearch.fcgi?db=nucleotide&term=$query&usehistory=y";
        my $output = get($url); #post the esearch URL
        $debug && print STDERR "$subn: \$acc=$acc \$output='$output'\n";

        #parse WebEnv and QueryKey
        my $web = $1 if ($output =~ /<WebEnv>(\S+)<\/WebEnv>/);
        my $key = $1 if ($output =~ /<QueryKey>(\d+)<\/QueryKey>/);
        $debug && print STDERR "$subn: \$acc=$acc \$web='$web' \$key='$key'\n";

        #assemble the efetch URL
        # Download the NC_* as a set
        my $efetch_url = $base . "efetch.fcgi?db=nucleotide&query_key=$key&WebEnv=$web";
        $efetch_url .= "&rettype=gb&retmode=text";
    #    $debug && print STDERR "$subn: \$acc=$acc \$efetch_url='$efetch_url'\n";
        my $efetch_out = '';
        for my $ntry (0 .. 2) {
            $efetch_out = get($efetch_url);
    #        $debug && print STDERR "$subn: \$acc=$acc \$efetch_out='$efetch_out'\n";
            last if ($efetch_out !~ /Unable to obtain query/);
            sleep(10);
        }
        return undef if ($efetch_out =~ /Unable to obtain query/);

        # open output file for writing
        open (my $OUT, ">$exe_dir/$outTemp") || die "$subn: Can't open file '$outTemp'!\n";
        print $OUT "$efetch_out";
        close $OUT || die "$subn: Can't close file '$outTemp'!\n";
        print STDERR "$subn: \$acc=$acc downloaded gb file and saved to $exe_dir/$outTemp\n";
    }
    # Get the gene symbol for all mat_peptides if there is any update in any mat_peptide
    # Need to update the gene symbol list, and put the genbank file in the refseq directory
    # Assumes the new file is in refseq/$outTemp, and old file in refseq/$outFileName
    $errcode->{noNewGenbank} = 0 if (!exists($errcode->{noNewGenbank})); # Set to 1 if there is no such genbank file
    $errcode->{noOldGenbank} = 0 if (!exists($errcode->{noOldGenbank})); # Set to 1 if there is no such genbank file
    $errcode->{newGenbank}   = $outTemp if (!exists($errcode->{newGenbank}));
    $errcode->{oldGenbank}   = $outFileName if (!exists($errcode->{oldGenbank}));
    $errcode->{newSymbols}   = [] if (!exists($errcode->{newSymbols})); # 'NC_007916|CDS=GI:90403549|Loc=1063..1920',

#    if ($errcode->{change_in_matpeptide}>0) {
        Annotate_Def::updateGeneSymbol( $exe_dir, $acc, $errcode);
#    }


#    my $bakFileName = "$outFileName.bak1";
    my $result = '';
    if (-e "$exe_dir/$outFileName") {
        $result = ` diff $exe_dir/$outFileName $exe_dir/$outTemp 2>&1 `;
        $debug && print STDERR "$subn: \$acc=$acc \$result='$result'\n";
    } else {
        $result = "$exe_dir/$outFileName: File not found";
    }
    if (!$result) {
        print STDERR "$subn: no difference between file $outFileName vs. $outTemp. not updated\n";
    } else {
      if ( $update ) {
        my $copies = Annotate_Util::backupFiles( $exe_dir, $outFileName, 8);
        print STDERR "$subn: \$acc=$acc subroutine backupFiles made $copies backups\n";
        ` mv $exe_dir/$outTemp $exe_dir/$outFileName `;
        print STDERR "$subn: $acc updated to $outFileName\n";
        $errcode->{updated} = 1;
      } else {
        print STDERR "$subn: $acc: No update is performed as \$update=$update. Please use '-update' to actually save RefSeq.\n";
      }
    }

    return ;
} # sub updateRefseq


1;
