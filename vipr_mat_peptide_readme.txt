
##//README//##
#
# vipr_mat_peptide.pl
#
# This script takes a viral genome in genbank format, uses a refseq to annotate any polyproptein within.
# It outputs a file named as <accession>_matpept_msa.faa with the annotated mat_peptides in fasta format
# if the result comes from alignment. Otherwise, outputs to a file <accession>_matpept_gbk.faa when the
# result comes from genbank.
#
# INPUT: directory of input file, genome file name
#
# OUTPUT: fasta file containing the annotated mat_peptide sequence
#
# DEPENDENCIES:
# This script calls perl and uses BioPerl modules.
#
# USAGE:
# For single genome
# ./vipr_mat_peptide.pl -d [dir_path] -i [inputFile.gb]
# For multiple input genomes within a directory
# ./vipr_mat_peptide.pl -d [dir_path] -l [directory]
# e.g.
# ./vipr_mat_peptide.pl -d ./ -i NC_001477_test.gb >> out.txt 2>> err.txt
# ./vipr_mat_peptide.pl -d ./ -l test >> test/out.txt 2>> test/err.txt
#
#    Authors: Guangyu Sun, gsun@vecna.com;
#    Authors: Chris Larsen, clarsen@vecna.com;
#    September 2011
#
#################

To-do:
1. What to do with the 2 RefSeqs of the Reoviridae family?

Changes

V1.3.1, Jan 2016
1. Added Zika virus (64320), Bagaza virus (64290), and Kedougou virus (64311). The reference mat_peptide
annotation was obtain from Arch Virol 2007, 152, 687-696.
2. Changed name to vipr_mat_peptide.

V1.3.0
1. Re-structured the program.

Changes in V1.2.1, May, 2013
1. Add the family to the output, so that it's clear what family the genome
belongs to.
2. Add function to include more reference CDS to the RefSeq CDS. The additional 
ref CDS were chosen as the centroid sequences in the groups clustered by 
uclust by RC Edgar, based on 90% identity search. This is done to make the
alignment align better in term of aligning the cleavage sites. Without the
additional CDS, alignment are more prone to error, such as following:
==========
ref=NC_001959|CDS=5..5374      ETEEEESEDQ IQMVPSDAVP E-G-KNKGKTKK GRGRKNNYNAFSRRGLSDEE
ACC=KC175410|CDS=3..5102       .......DDE EFVISSDDIK T-E-GKKGKNKT GRGKK..HTAFSSKGLSDEE
==========
2581 aagatgatga ggagttcgtc atttcatctg acgacattaa aact-gag-ggt aagaaaggga
===>                                    DDDIIIKK KTTT-EEE-GGG KKKKKKGGG
==========
  After adding 4 more polyproteins, the alignment form most of the genomes are
correct. Exceptions include: AF097917 and AY126474, where the correct cleavage
site can't be reliably determined based on alignment.

3. Added MERS to the list of approved species. The refseq NC_019843 is
annotated based on SARS refseq NC_004718 and the annotation from UniProt.
4. Debugged the MSA=*..* tag in the note of the features annotated from
different CDS based on same refseq mat_peptide, in the annotation of
NC_019843 based on NC_004718. However, this turned out to be because of the
different cleavage site. Since in the real run of the script, the refseq is
closer to the target, this issue is not further pursued.

Changes in V1.2.0, April, 2013
1. Add a third nucleotide if the last codon has 2 nucleotides that can
determine the amino acid already, in order to overcome the change from BioPerl
1.5.8 to 1.6.1 where the default behavior to deal with partial codon changed.
2. Added 13 species in the Flaviviridae to the apporved RefSeq list.
3. For family Reoviridae, there are 2 RefSeqs that have mat_peptides:
NC_004278 and NC_003758. What to do with these?

Changes in V1.1.8, March 19, 2012
1. Tweaked the sub get_new_translation, in order to handle cases like
AF379201, and EF424619.
=====================
 Test #122, Flaviviridae,   AF379201.gb,  $result=diff
./test/output/AF379201_matpept_msa.faa ./test/AF379201_matpept_msa.faa
4c4
< .KTYTTGGVVGRSTSGLASFLSPGPQQKIQ
---
> XKTYTTGGVVGRSTSGLASFLSPGPQQKIQ
=====================
2. Changed the way duplicate mat_peptides are detected, so that in cases like
following where 2 mat_peptide with identical locations are created but are from 
different reference mat_peptides, only the one without "Partial=Y" label is
listed. The sub Annotate_Align::combineFeatures was created to handle this.
=====================
> src=MSA|ACC=AJ299464|Ver=AJ299464.3|CDS=GI:74381881|ref=NC_001489.1|RM=741..803|Loc=644..706||AA=1..21|symbol=1A(VP4b)|Partial=N|product=1A VP4b mature peptide (alt.)|*new*
MSRQGIFQTVGSGLDHILSLA
> src=MSA|ACC=AJ299464|Ver=AJ299464.3|CDS=GI:74381881|ref=NC_001489.1|RM=735..803|Loc=644..706||AA=1..21|symbol=1A(VP4a)|Partial=Y|product=1A VP4a mature peptide (alt.)|*new*
MSRQGIFQTVGSGLDHILSLA
=====================
3. Improved the way to update the taxon and RefSeqs from NCBI. See following
options:
   -update: save the result to file. Default is only to output to screen.
   -checkrefseq=n: to check the refseqs. n=1: use saved genbank files for each
family; 2: live download from genbank server.
   -checktaxon=n: to check the taxon info. n=1: use saved xml files for each
family; 2: live download from genbank server.
4. There are noticed "product" name change in some Refseqs in Togaviridae.
Examples include NC_001449.

Changes in V1.1.7, February 05, 2012
1. Added NC_018138 of Bunyaviridae to approved list of refseqs.
2. Separated the 3 versions of 1C and 1D in NC_001430 by symbols 1Cv1, 1Cv2 
and 1CV3, so that all version of the mat_peptides are retained.
3. For a split refseq mat_peptide, requires the resulting mat_peptide to be
also split, unless the resulting mat_peptide is shorter than or equal to the
longest fragment of the refseq mat_peptide. This results in potential changes 
in all species/families with split mat_peptide.
4. Added commandline option -checktaxon to check any update for the taxon 
structure, as this has occurred for Coronaviridae recently. Sub checkAllTaxon 
is added for this.  To actually download the taxon from NCBI, the variable
$DOWNLOAD_TAXON=1 needs to be set in sub checkAllTaxon. Otherwise the stored
xml files would be used (this is mainly because of Caliciviridae that has
14000+ strains that causes problem during download). Then the new taxon is
read from the xml file for each family and compared with those read from
Annotate_taxon_records.txt. Any difference will be noted.
  With this feature and the -checkrefseq option, any update from NCBI regarding
the taxon and refseq could be detected. Notice, these checking operations need to 
be managed by a competent administrator. Therefore, the final process isn't 
automated.

Changes in V1.1.6, November 16, 2012
1. Added following new family/species/strain.
1a. Added the family Picornaviridae, which has 69 strains/species. 
1b. Also added 8 strains/species to family Caliciviridae in addition to 
species Norwalk virus, 11983.
2. Added the location of reference mat_peptide to the defline in output fasta
file. Only the start/stop are included as this is for identification purpose
only. Example is below:
=====================
> src=MSA|ACC=AM710685|Ver=AM710685.1|CDS=GI:146349919|ref=NC_004102.1|reffeat=915..1490|Loc=1..169|cstart=2|AA=1..56|symbol=E1|Partial=Y|product=E1 protein
PTTALVVSQVLRIPQAIVDMVAGAHWGVLMGLAYYSMVGNWAKVLVVLLLFAGVDA
=====================
3. Modified the handling of codon_start with complement location, such as in
AF211032. Added capability to handle complement strand with codon_start in RefSeq for
Simple and Split reffeat. Need review for Fuzzy location. Example:
=====================
> src=MSA|ACC=AF211032|Ver=AF211032.1|CDS=GI:6856791|ref=NC_004102.1|reffeat=915..1490|Loc=complement(100..155)|cstart=3|AA=1..18|symbol=E1|Partial=Y|product=E1 protein|*new*
GELAKVLIVMLLFAGVDG
=====================
4. Added -checkrefseq option, to check any update in RefSeq from NCBI, then
exit. With this option, each family is searched, then the RefSeqs are
downloaded in genbank format. The mat_peptides are then compared with from
files stored with the script. Any change is printed out. Sample output:
====================
downloadRefseq: $fam=Hepeviridae identical RefSeq_Hepeviridae.gbk.bak1 and RefSeq_Hepeviridae.gbk $result=''
checkAllRefseq: $nseq=1 $acc=NC_015521	newf='NC_015521|CDS=GI:335352416|Loc=284..835'
checkAllRefseq: $nseq=1 $acc=NC_015521	newf='NC_015521|CDS=GI:335352416|Loc=1388..1897'
checkAllRefseq: $nseq=1 $acc=NC_015521	newf='NC_015521|CDS=GI:335352416|Loc=2999..3763'
checkAllRefseq: $nseq=1 $acc=NC_015521	newf='NC_015521|CDS=GI:335352416|Loc=3770..5221'
checkAllRefseq:   family        accession       has_mat no_file changed taxid	speciesid       species
checkAllRefseq: Hepeviridae     NC_018382       -       -       -	1216472 -1      ''
checkAllRefseq: Hepeviridae     NC_001434       -       -       -       12461	12461   'Hepatitis E virus'
checkAllRefseq: Hepeviridae     NC_015521       4       1       4	1016879 1016879 'Cutthroat trout virus'
checkAllRefseq: Hepeviridae     total   3       1       1       1
====================
5. When no input is given, print out the list of viral families, species/strains,
 and their refseq in the usage in following format:
====================
Family=Bunyaviridae
Bunyaviridae     #1      species=11591   strain=-----	ref=NC_005220        "Uukuniemi virus"
Bunyaviridae     #2      species=11599   strain=-----	ref=NC_005219        "Hantaan virus"
====================
##### Change #6 has been reverted in V1.1.7 #####
6. The criteria for a good mat_peptide in term of internal gaps is: 
   1). if the run of good alignment is longer than 10AA, and
   2). the total of gapped alignment is less than 33% of total alignment.
   Examples include AM710685.
   With this change, it's necessary to re-compute all mat_peptides as any
partial mat_peptide could be affected.
   One problem with this change is Genome AJ507016, where the 16AA at N end
doesn't look like good fit for "core protein". And there are many like this.
##### Change #6 has been reverted in V1.1.7 #####
7. Internally, changed the way the sequences around the gap is extracted by 
taking the values while the feature location is being calculated, to get the 
actual sequences from the alignment.
8. Internally, created module, Annotate_Def.pm, for functions dealing with 
standard data such as RefSeq, taxonomy, and gene_symbol.


Changes in V1.1.5, October 09, 2012
1. To eliminate any extremely short sequence to be presented as mat_peptide,
included a minimum requirement (10AA or 50% of reference mat_peptide if
shorter) to exclude such. This change would prevent those short sequences to
be identified as mat_peptides, while allowing sequences >10AA or 50% of
reference mat_peptide to be included as mat_peptide.
2. Switched the alignment program to CLUSTALW. This program has been shown to
create less gap, thus less segments, when the ends of the sequences have large
number of mutations. This is compared with MUSCLE, which tends to create more
gaps.
3. To handle the incomplete CDS where codon_start is not 1, it is decided
during 10/9/2012 Weekly Science Meeting that any such mat_peptide should
follow the CDS and use codon_start. Example output is:
=====================
> src=MSA|ACC=AM272362|Ver=AM272362.1|CDS=GI:119668529|ref=NC_004102.1|Loc=1..169|cstart=2|AA=1..56|symbol=E1|Partial=Y|product=E1 protein
PTTALVVSQLLRIPQAVVDMVAGAHWGVLAGLAYYSMVGNWAKVLIVMLLFAGVDG
=====================
4. BioPerl version 1.006001 is required. It has been observed that BioPerl 
v1.006901 doesn't translate trailing codon with only 2 nucleotides. Since 
such codons are translated in genbank CDS, it's better to follow this 
practice. It is shown BioPerl v1.006001 does this translation. 


Update of V1.1.4, June 27, 2012
1. Added Coronaviridae family to the list of covered viral families. This
means that the refseqs have been reviewed and all suitable (with mat_peptide
annotation) have been added the the script. This is aimed mainly at SARS, but
also covers other species within the Coronaviridae family.
2. Code change mostly aimed at finding out and handling the duplicate mat_peptides
resulted from >1 CDS from input genome.

Features of vipr_mat_peptide.pl V1.1.3
Perl script: vipr_mat_peptide.pl, 
Modular design: functions grouped into 6 modules for easy upkeep
  Easy input/output
  Takes genbank file as input, output mat_peptides in fasta file with name <accession>_matpept_msagbk.faa. Also writes messages to stdout and stderr for monitoring/debugging. 
  Validated species/refseq via opt-in process
  MSA-based annotation is only available for genomes in the species that have been validated. For other species, the annotation in genbank is taken if present, or left blank otherwise
  Easy addition of more validated species
  Simply add the genbank file to refseq/ directory, add a line in sub Annotate_Util::get_refseq_acc (format: speciesid=>accession, this should probably be changed to a file, instead of a subroutine), and define the gene symbols in the Annotate_symbol_records.txt.


###################
Why annotate the mat_peptide in viral genomes

Some viruses code all propteins in one long coding sequence (CDS), which produces a polyprotein after translation
The polyprotein is then processed to generate final products, the mature peptides, by proteinase. These mat_peptides are the real functional proteins
Many viral genomes in genbank lack the mat_peptide annotation, which is essential to study the sequence, structure, and function of a protein, as mat_peptide is the functional unit
Other genomes might have mat_peptide annotation that doesn't confirm with latest development

Using refseq to annotate mat_peptides

Many viral species have refseqs in genbank w/ mat_peptide annotation, which can be used as standard to annotate other genomes of the same species
The refseqs in genbank have been created by expert biologist as well as NCBI staff. 
Some have undergone review, while others are provisional
Most have complete sequence in full length

Annotation Process

With a normal genome, one would do following:
1. find an appropriate refseq,
2. find the corresponding CDS,
3. align the CDS, then
4. find the start/end of each mat_peptide in the target,
5. define the mat_peptides in the target, and finally
6. save the annotated mat_peptide.

Workflow of the Script

1. find an appropriate refseq:
First determin the taxon id of the target genome, then look for a defined refseq for that strain. If no refseq for that strain, go up to the species level and look for any defined refseq for the species
2. find the corresponding CDS:
There might be >1 CDSs in either target or refseq. To find the corresponding CDSs, run blast with a target CDS and a refseq CDS. Anything with <50% conserved residues is ignored
3. align the CDS:
Ran global alignment with CLUSTALW (MUSCLE was used before Oct 2012).
Blast produces local alignment, which might leave out any overhang, or mis-matched terminals
Also considered hidden Markov model, but this is more of a scoring method than a alignment method
4. find the start/end of each mat_peptide in the target:
Based on the alignment of target CDS and refseq CDS, find the start/end of each mat_peptide in the target genome, taking into consideration following factors: different start of the polyprotein, terminal gaps in MSA, and internal gaps in MSA, any skip/buckle in either refseq/target CDS. Then convert the AA range back to DNA range
5. define the mat_peptides in the target:
For each mat_peptide, identify with genome accession (with version), CDS id, refseq accession, DNA location in genome, AA location in polyprotein, gene symbol, and product name
6. save the annotated mat_peptide:
Save the annotation in fasta for database loading. May also save updated genbank to file for other purpose


Validate more species/families: 

Add capability to cover other families of viprbrc.org:
SS+ RNA
 (x) Caliciviridae
 (x) Coronaviridae
 (x) Flaviviridae
 (*) Hepeviridae (problems with the 2 available refseqs)
 (x) Picornaviridae
 (x) Togaviridae
SS- RNA
 (x) Arenaviridae
 (x) Bunyaviridae
 ( ) Filoviridae
 ( ) Paramyxoviridae
 ( ) Rhabdoviridae
DS RNA
 ( ) Reoviridae
DS DNA
 ( ) Herpesviridae
 ( ) Poxviridae

Validated species (Most up-to-date numbers are available from the script)
So far, validated species include:
Flaviviridae: 
  Flavivirus, 27 strain/species
  Hepacivirus, 2 strain/species
  other species.

Caliciviridae: 
  Norovirus, 2 strain/species

Togaviridae: 
  Alphavirus, 15 strain/species
  Rubivirus, 1 strain/species

Coronaviridae:
  10 strain/species
  SARS is the focal species, but other 7 species are also covered.

Arenaviridae:
  1 strain/species
  Only Junin virus has a good refseq with mat_peptide for the segment S. It
has 2 mat_peptides for glycoprotein precursor.

Bunyaviridae:
  2 strain/species
  The genomes in this family is divided into 3 segments: L, M, and S. Segment
M has a polyprotein, at least for 3 species: 
NC_005220: Uukuniemi virus. Has 2 mat_peptides, and sig_peptide. Included in
V1.1.4
NC_005219: Hantaan virus. Has 2 mat_peptides, and sig_peptide. Included in
V1.1.4
NC_004158: Dugbe virus. Has only 1 mat_peptide. Excluded from V1.1.4

Hepeviridae:
  Found 2 refseqs:
=======================================================================
|   12461 | NC_001434 |              181 | Hepatitis E virus     |   12461 |
No mat_peptide in refseq
| 1016879 | NC_015521 |                2 | Cutthroat trout virus | 1016879 |
Has 4 mat_peptide, but has gaps at every cleavage site
=======================================================================
  Because of the problems with the refseqs, don't annotate this family.

Picornaviridae:
  Found 70 refseqs (NC_0*) for this family, of which 56 have mat_peptides in
their genbank file.
  This family has many common viruses, such as polio, foot-and-mouth diease,
enterovirus, etc.
  There are quite some mat_peptide annotation for Foot-and-mouth diease,
however, they are of low quality, some might be from 20 years ago. The newly
available refseqs should improve the quality of mat_peptide annotation.


Further direction

Question:
  With the current design/implementation of vipr_mat_peptide.pl script, it's relatively easy to add new species. However, much time has been spent in checking if the refseq is valid, and if the refseq works with all genomes in the given species. Though some minor problem has been observed, such as GT3 of hepacivirus as compared with GT1, most refseqs are used as-is without much change. Is it necessary to do extensive checking going forward?
Answer:
  There were quite a few problems found duing update for Coronaviridae (Version 1.1.4). So, at least
for now, careful review is needed.

Change the list of validated species to a file, instead of subroutine 

Add capability to automatically check refseq against genbank 
Important for long-term use.

Some refseqs have sig_peptide in addition to mat_peptide, such as in
NC_005220:
=======================================================================
     sig_peptide     18..68
                     /locus_tag="UUKVsMgp1"
                     /note="membrane glycoprotein G1 signal peptide"
     mat_peptide     69..1505
                     /locus_tag="UUKVsMgp1"
                     /product="membrane glycoprotein G1"
                     /protein_id="NP_941985.1"
                     /db_xref="GI:38371705"
     sig_peptide     1506..1556
                     /locus_tag="UUKVsMgp1"
                     /note="membrane glycoprotein G2 signal peptide"
     mat_peptide     1557..3041
                     /locus_tag="UUKVsMgp1"
                     /product="membrane glycoprotein G2"
                     /protein_id="NP_941986.1"
                     /db_xref="GI:38371706"
=======================================================================
Should this be also included in the output?

Other request?


