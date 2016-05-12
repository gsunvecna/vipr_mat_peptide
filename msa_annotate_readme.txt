
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
#    Authors: Chris Larsen, clarsen@vecna.com; Guangyu Sun, gsun@vecna.com;
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

Further direction
Validate more species/families: 
Queation: with the current design/implementation of msa_annotate.pl script, it's relatively easy to add new species. However, much time has been spend in checking if the refseq is valid, and if the refseq works with all genomes in the given species. Though some minor problem has been observed, such as GT3 of hepacivirus as compared with GT1, most refseqs are used as-is without much change. Is it necessary to do extensive checking going forward?
Change the list of validated species to a file, instead of subroutine 
Add capability to automatically check refseq against genbank 
Important for long-term use.
Other request?
