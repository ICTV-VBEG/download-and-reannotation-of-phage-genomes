# Snakemake workflow: `download-and-reannotation-of-phage-genomes`

[![Snakemake](https://img.shields.io/badge/snakemake-≥6.3.0-brightgreen.svg)](https://snakemake.github.io)
[![GitHub actions status](https://github.com/ICTV-VBEG/download-and-reannotation-of-phage-genomes/workflows/Tests/badge.svg?branch=main)](https://github.com/ICTV-VBEG/download-and-reannotation-of-phage-genomes/actions?query=branch%3Amain+workflow%3ATests)


A Snakemake workflow for `A snakemake workflow for downloading and reannotating phage genomes. The workflow is based on prior work by Dann Turner.`


## Usage

The usage of this workflow is described in the [Snakemake Workflow Catalog](https://snakemake.github.io/snakemake-workflow-catalog/?usage=ICTV-VBEG%2Fdownload-and-reannotation-of-phage-genomes).

If you use this workflow in a paper, don't forget to give credits to the authors by citing the URL of this (original) download-and-reannotation-of-phage-genomessitory and its DOI (see above).

# TODO for workflow implementation

This workflow uses a combination of bash scripts, some sourced online and others written by DT.

## 1. Download genomes using Efetch
This script uses Efetch from the NCBI EUtilities package to download fasta records for a list specified in a txt file - one accession per line. Ryan Cook has probably got a much better implementation in his Inphared pipeline

```python
#!/usr/bin/env python3
import sys
from Bio import Entrez

#define email for entrez login
db           = "nuccore"
Entrez.email = "some_email@somedomain.com"

#load accessions from arguments
if len(sys.argv[1:]) > 1:
  accs = sys.argv[1:]
else: #load accesions from stdin
  accs = [ l.strip() for l in sys.stdin if l.strip() ]
#fetch
sys.stderr.write( "Fetching %s entries from GenBank: %s\n" % (len(accs), ", ".join(accs[:10])))
for i,acc in enumerate(accs):
  try:
    sys.stderr.write( " %9i %s          \r" % (i+1,acc))
    handle = Entrez.efetch(db=db, rettype="fasta", id=acc)
    #print output to stdout
    sys.stdout.write(handle.read())
  except:
    sys.stderr.write( "Error! Cannot fetch: %s        \n" % acc)
 ```
 
 Usage:
 ```bash
 ./efetch_fasta.py < accessions.txt > accessions.fasta
 ```
 
 ## 2. Split multi-fasta file into separate records
 I rename the fasta records in the file as some annotation tools have issues with fasta ids over a certain number of characters or which contain whitespace
 ```bash
 sed -i accessions.bak 's/\.[0-9].*//g' accessions.fasta
 ```
 
 The following script splits a multi-fasta file into individual files, named by the sequence id:
 
 ```python
#!/usr/bin/env python3
import os
from Bio import SeqIO

def split(fastafile = "test_fasta.fasta",
          outfastadir = "splitoutput"):
    """Extract multiple sequence fasta file and write each sequence in separate file"""
    os.system("mkdir -p %s"% (outfastadir))
    with open (fastafile) as FH:
        record = SeqIO.parse(FH, "fasta")
        file_count = 0
        for seq_rec in record:
            file_count = file_count + 1
            header = seq_rec.id
            with open("%s/%s.fasta" % (outfastadir,str(header)), "w") as FHO:
                SeqIO.write(seq_rec, FHO, "fasta")
    if file_count == 0:
        raise Exception("No valid sequence in fasta file")
    return "Done"

if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser()

    parser.add_argument('-f','--fastafile',
                        action  ="store",
                        default ="test_fasta.fasta",
                        help="Fasta File for parsing")
    parser.add_argument('-d','--outfastadir',
                        action  ="store",
                        default ="splitoutput",
                        help    ="Fasta File output directory")

    args = parser.parse_args()
    split(fastafile     =   args.fastafile,
          outfastadir   =   args.outfastadir)
 ```
 
 Usage: 
 ```bash
 ./split_multifasta.py -f accessions.fasta -d split_output
 ```
 
## 3. Reannotate genomes
This is simply a step to ensure that there is consistency across the gene calls. There are both pros and cons to this step.
 
Annotation can be performed using either [Prokka](https://github.com/tseemann/prokka) or [Phannotate](https://github.com/gbouras13/pharokka) - whichever is the user's preference. For the former, a custom database is set up using a [PHROGs HMM db](https://s3.climb.ac.uk/ADM_share/all_phrogs.hmm.gz).

The bash script is simply a call to the annotation pipeline, but takes the base name of the fasta file and passes this as a variable for naming the locus_tag as well as the output directory and files.

**Option 1: Prokka** 
```bash
#!/bin/bash
for file in *.fasta; do
	tag=${file%.fasta};
	string=${tag##*_};
	echo $string	
	prokka --hmms /bioinformatics/prokka/db/hmm/all_phrogs_ann.hmm --locustag $string --outdir $tag --prefix $tag --gcode 11 --cpus 8 --rfam $file
done
```

**Option 2: Pharokka**
```bash
#!/bin/bash
for file in *.fasta; do
	tag=${file%.fasta};
	string=${tag##*_};
	echo $string	
	pharokka -i $file -l $string -o $tag -p $tag -t 8
done
```

Usage: script is in PATH or copied to directory containing the split FASTA files 
```
./batch_pharokka.sh
```

## 4. Collate reannotated GFF or FAA files for pan-genome analysis
Recursively copy all .gff and .faa files to an appropriate directory for downstream analysis
```
cp **/*.gff genomes_gff/
cp **/*.faa genomes_faa/
```

## 5. [OPTION] Pangenome analysis using PIRATE
Again, a variety of tools are available here. Note that none are specifically designed for phage genomes, rather they are focussed on bacterial core genomes.
Some examples are [PIRATE], [Panaroo], [Roary]

For ease, this example uses PIRATE with multiple % id thresholds for pangenome construction
```
PIRATE -i genomes_gff/ -o genomes_pirate --min-len 10 -s 35,40,45 -t 24
```

In case PIRATE changed any of the locus_tag values, and to ensure that it matches back to the sequence identifiers in our .faa files, this utility script is run

```
subsample_outputs.pl -i genomes_pirate/PIRATE.gene_families.tsv -g genomes_pirate/modified_gffs/ -o genomes_pirate/genomes_pirate_gene_families.prev_locus.tsv --field "prev_locus"
```

At this point, some manual processing is required. I create a copy of the genomes_pirate_gene_families.prev_locus.tsv file, sort by number_genomes column and transform the data to a binary presence/absence matrix

## 6. [OPTION] Using MMSeqs2 to cluster proteins

Using a concatenated file of all of the proteins, create an MMSeqs database. There is code within PanACoTA which also performs these steps.
```
cat *.faa > proteins.faa
```

```
mmseqs createdb proteins.faa protein_db --createdb-mode 1
```

Cluster the proteins
```
mmseqs cluster proteins_db proteins_clust_db --cluster-mode 1 -c 0.7 --alignment-mode 3 --min-seq-id 0.7
```

Convert the cluster database to a tsv file
```
mmseqs createtsv proteins_db proteins_db proteins_clust_db clusters.tsv
```

I started to write a script to summarise the clusters from long to short format, sorted by genome accession as columns. Could add to this to convert the dataframe to a count or binary.

```python
#!/usr/bin/env python3
import os
import sys
import numpy as np
import pandas as pd

def genomes(df):
         #Create a dataframe consisting of genome accession numbers
         df['genomes'] = df['members'].str.replace('_.*', '', regex=True)
         df2 = df.pivot_table(index='rep_seq', columns='genomes', values='members', aggfunc=lambda x: ', '.join(str(v) for v in x))
         df2.to_csv('clusters_by_accession.tsv', sep='\t', encoding='utf-8', index=False)
         print(df.columns.values)
         #print(df.to_markdown())
         return df


mmseq_df = pd.read_csv('results.tsv', header=None, sep='\t')
mmseq_df.columns = ['rep_seq', 'members']
genomes(mmseq_df)

```

## 7. Extract protein sequences correlating to a "signature" gene product
We can concatenate all of the proteins called by the annotation pipeline
```
cat *.faa > proteins.faa
```

Then, using a list of accessions representing a conserved protein group these can be pulled from the proteins.faa file into a separate file for alignment and phylogeny estimation.

The following script was developed by John Nash at the NRCC. I've not yet found the time to write a python solution as this works.

```perl
#! /usr/bin/perl
use warnings;
use strict;
use Text::Wrap;

## Program Info:
#
# Name:  subset_fasta
#
# Function:  Takes a multiple fasta file and removes a set of 
#    sequences to makes a second fasta file.  Useful for pulling
#    subsets of sequences from entire genomes.
#
# Author: John Nash
#  Copyright (c) National Research Council of Canada, 2000-2003,
#  all rights reserved.
#
# Licence: This script may be used freely as long as no fee is charged
#    for use, and as long as the author/copyright attributions
#    are not removed.
#
# History:
#   Version 1.0 (June 04, 2001): first non-beta release.
#   Version 1.1 (June 05, 2001): Fixed up header of subset file to 
#     give more info
#   Version 1.2 (June 26, 2001): used stdin instead of a file from cmd line
#     Can define sequences to be retrieved on an optional command line list
#   Version 1.3 (Oct 4, 2002): Preserves gene order from incoming list
#   Version 1.4 (Oct 7, 2002); better behaved with ">" characters!
#   Version 1.5 (July 22, 2005); Take name from first word of 
#     header of input (-i) file.
#     No more dependence on Tie:IxHash 
#
##

my $title = "subset_fasta";
my $version = "1.5";
my $date = "22 July, 2005";
$Text::Wrap::columns = 65;

# Error message:
my $error_msg = "Type \"$title -h\" for help.";

# Get and process the command line params:
# Returns array of $fasta_file and $orf_file;
my (@cmd_line, $orf_file, $opt_orfs);

@cmd_line = ();
@cmd_line = process_command_line();
$orf_file = $cmd_line[0];
$opt_orfs = $cmd_line[1];

my @orf_list;
# handle the possiblilty of command-line ORFS:

if ($opt_orfs eq '')  {
# open the ORF list:
  open ORF_FILE, $orf_file or 
    die "Cannot open $orf_file: $!\n$error_msg\n";
	
# read each value into an array:
  while (<ORF_FILE>) {
    s/\r\n/\n/g;
    chomp;
# We don't like blank lines!!!
    next if (/^\s*$/);
    my @temp = split / /, $_;
    $_ = $temp[0];

# if there is a fasta header char ( > ), remove it:
    $_ = substr ($_, 1, length $_) if (/^>/);
    push @orf_list, $_;
  }
}
else  {
  @orf_list = split /:/, $opt_orfs;
}

## Read in the sequence from a multiple FASTA file:
my %sequence;
my $seq_name;

while (<>) {
  s/\r\n/\n/g;
  chomp;
# read the whole header as we display it later:
  if (/^>/)  {
# Take the header as the sequence name:
# remove the header character ( > ):
    $seq_name = substr($_, 1, length $_);
  }
  else {
    $sequence{$seq_name} .= $_;
  }
}

# if the header is in the ORF_list
# take the sequence
# print the header and sequence:
my $orf_name;
my $hit;

foreach $hit (sort keys %sequence) {

# grab the first word which is the sequence name:
  my @temp = split / /, $hit;
  my $temp_name = $temp[0];

# Compare with list
  foreach $orf_name (@orf_list) {
    if ($orf_name eq $temp_name) {

# Add the full header, not just the sequence name:
      print ">$hit\n";
      print wrap('', '', "$sequence{$hit}\n");
    }
  }
}

### end of main:

##### SUBROUTINES:
sub process_command_line {
#
# Expects: 
# Returns: @my_cmd_line = ($orf_file, $opt_orfs)
# Uses:
	
# Variables:
  my %opts = ();    # command line params, as entered by user
  my @my_cmd_line;  # returned value
  my @list;	    # %opts as an array for handling
  my $cmd_args;     # return value for getopts()
	
# Holders for command line's files:
  my $myorf_file = '';
  my $myopt_orfs = '';
	
# Scratch:
  my $item;
	
# Get the command=line parameters:
  use vars qw($opt_f $opt_i $opt_o $opt_h);
  use Getopt::Std;
  $cmd_args = getopts('i:o:h', \%opts);
	
# Die on illegal argument list:
  if ($cmd_args == 0) {
    die ("Error: Missing or incorrect command line parameter(s)!\n",
         $error_msg, "\n");
  }
	
# Check and enact each command-line argument:
  if (!%opts)  {
    die ($error_msg, "\n");
  }
	
# Make the hashes into an array:
  @list = keys %opts;
	
# Do a quick check for "help" and the compulsory parameters:
#   If the compulsory files are not there, squawk and die:
  foreach $item (@list)  {
# Help:
    if ($item eq "h")  {
      help();
    }
# ORF file:
    elsif ($item eq "i") {
      $myorf_file = $opts{$item};
    }
# optional command-line ORFs:
    elsif ($item eq "o") {
			$myopt_orfs = $opts{$item};
    }
  }
	
# Put it in an array:
	@my_cmd_line = ($myorf_file, $myopt_orfs);
	return @my_cmd_line;
	
} #end of sub process_command_line()

sub help {
	
print <<EOHelp;
$title, version $version, $date

Function:  Takes a multiple fasta file and removes a set of 
      sequences to makes a second fasta file.  Useful for pulling
      subsets of sequences from entire genomes. The input fasta file
      comes from stdin, and output is written to stdout, so $title can 
      be a filter. 

Syntax:  $title -i list_file < fasta_file > subset_file
or:   $title -o header1:header2:header3 < fasta_file > subset_file
or:   $title -h for help

Compulsory arguments:

Either:
 -i   list_file:  a list of sequences, one per line, of headers of
      sequences to be retrieved from the larger file.
or:
 -o   header1:header2:headerN which is a colon-delimited set of headers
      of sequences to be retrieved from the larger file. This parameter
      takes priority over \"-i\".  If both are supplied, \"-o\" is used
      preferentially.


Example:
 $title -i myorfs < genome.fasta > subset.fasta
      where genome.fasta is the multiple fasta file containing lots of 
      genes, and subset.fasta is the new output.

**Warning**  Make sure the headers in your list are the same as the
      headers in the fasta file !!!  The characters from the \> up to 
      the first space in each entry are used as the \"sequence name\".

EOHelp
die ("\n");
} # end of sub help
```

## 8. Signature gene alingment and calculating ML phylogeny using IQTree2
This step does require manual input if the user wishes to add an outgroup sequence to use later as a root for the tree.  

There are a variety of different packages that could be used for the alignment step e.g. [MAFFT](https://mafft.cbrc.jp/alignment/software/), [Clustal Omega](http://www.clustal.org/omega/) etc.

For a rooted ML phylogeny, I generally use the following settings in IQTree2 to construct an ML tree with 1000 UFBoot and SH-alrt replicates.

```
 iqtree2 -s alignment.file -o outgroup_sequence_id -st AA -m MFP -nt AUTO -mem 128G -bb 1000 -alrt 1000 -cptime 120
 ```
 
The output .treefile can then be visualised and annotated using the iTOL webservice
