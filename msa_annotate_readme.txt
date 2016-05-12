
##//README//##
#
# msa_annotate.pl
#
# This script takes a viral genome in genbank format, uses a refseq to annotate any polyproptein within.
# It outputs a file named as <accession>_matpept_msagbk.faa with the annotated mat_peptides in fasta format.
#
# INPUT: directory of input file, genome file name
#
# OUTPUT: fasta file containing the annotated mat_peptide sequence
#
# DEPENDENCIES:
# This script calls perl and uses BioPerl modules.
## Path to the MUSCLE binaries. You need to configure this, if muscle is not already in the path
#       BEGIN {$ENV{MUSCLEDIR} = '/net/home/gsun/prog/muscle/mus37'}
# Specify the MUSCLE executable location in your environment as directed above!
#
# USAGE:
# For single genome
# ./msa_annotate.pl -d [dir_path] -i [inputFile.gb]
# For multiple input genomes within a directory
# ./msa_annotate.pl -d [dir_path] -l [directory]
# e.g.
# ./msa_annotate.pl -d ./ -i NC_001477_test.gb >> out.txt 2>> err.txt
# ./msa_annotate.pl -d ./ -l test >> test/out.txt 2>> test/err.txt
#
#    Authors: Chris Larsen, clarsen@vecna.com;
#             Guangyu Sun, gsun@vecna.com;
#    September 2011
#
#################

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
First determin the taxon id of the target genome, then look for a defined refseq for that strain. If no refseq for that strain, go up one level and look for any defined refseq for the species
2. find the corresponding CDS:
There might be >1 CDSs in either target or refseq. To find the corresponding CDSs, run blast with a target CDS and a refseq CDS. Anything with <50% conserved residues is ignored
3. align the CDS:
Ran global alignment with MUSCLE. 
Blast produces local alignment, which might leave out any overhang, or mis-matched terminals
Also considered hidden Markov model, but this is more of a scoring method than a alignment method
4. find the start/end of each mat_peptide in the target:
Based on the alignment of target CDS and refseq CDS, find the start/end of each mat_peptide in the target genome, taking into consideration following factors: different start of the polyprotein, terminal gaps in MSA, and internal gaps in MSA, any skip/buckle in either refseq/target CDS. Then convert the AA range back to DNA range
5. define the mat_peptides in the target:
For each mat_peptide, identify with genome accession (with version), CDS id, refseq accession, DNA location in genome, AA location in polyprotein, gene symbol, and product name
6. save the annotated mat_peptide:
Save the annotation in fasta for database loading. May also save updated genbank to file for other purpose

Features of msa_annovate.pl V1.1.3

Perl script: msa_annotate.pl, 
Modular design: functions grouped into 6 modules for easy upkeep
  Easy input/output
  Takes genbank file as input, output mat_peptides in fasta file with name <accession>_matpept_msagbk.faa. Also writes messages to stdout and stderr for monitoring/debugging. 
  Validated species/refseq via opt-in process
  MSA-based annotation is only available for genomes in the species that have been validated. For other species, the annotation in genbank is take if present, or left blank otherwise
  Easy addition of more validated species
  Simply add the genbank file to refseq/ directory, add a line in sub Annotate_Util::get_refseq_acc (format: speciesid=>accession, this should probably be changed to a file, instead of a subroutine), and define the gene symbols in the Annotate_symbol_records.txt.

Update of V1.1.4, June 27, 2012
1. Added Coronaviridae family to the list of covered viral families. This
means that the refseqs have been reviewed and all suitable (with mat_peptide
annotation) have been added the the script. This is aimed mainly at SARS, but
also covers other species within the Coronaviridae family.
2. Code change mostly aimed at finding out and handling the duplicate mat_peptides
resulted from >1 CDS from input genome.


Validated species
So far, validated species include:
Flaviviridae: 
  Flavivirus, 27 strain/species
  Hepacivirus, 2 strain/species

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


Further direction
Validate more species/families: 

Add capability to cover other families of viprbrc.org:
SS+ RNA
 (x) Caliciviridae
 (x) Coronaviridae
 (x) Flaviviridae
 ( ) Hepeviridae
 ( ) Picornaviridae
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

Question:
  With the current design/implementation of msa_annotate.pl script, it's relatively easy to add new species. However, much time has been spent in checking if the refseq is valid, and if the refseq works with all genomes in the given species. Though some minor problem has been observed, such as GT3 of hepacivirus as compared with GT1, most refseqs are used as-is without much change. Is it necessary to do extensive checking going forward?
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


