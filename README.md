# samToPolyA
A utility to detect poly-adenylated sequencing reads, call on-genome polyA sites and infer the reads' strand based on reads-to-genome alignments in SAM format

## Input
Read-to-genome alignments in SAM format, and the corresponding genome sequence in multifasta.

## Options
This script maps polyA sites on the genome based on read mappings in SAM format, and according to following provided parameters:
  - `minClipped` (integer)
            = minimum length of A or T tail required to call a PolyA site. 
            **Default**: '10'.
  - `minAcontent` (float)
            = required A (or T, if minus strand) content of the tail. 
            **Default**: '0.8'.
              Note: minAcontent affects both the A tail and the upstream A stretch.
  - `minUpMisPrimeAlength` (integer)
            = minimum length of genomic A stretch immediately upstream a putative site required to call a false positive (presumably due to internal RT priming), and hence not report the corresponding site.
            **Default**: '10'.
  - `genomeFasta` (string)
            = path to multifasta of genome (+ spike-in sequences if applicable), used to extract upstream genomic sequence.

## Output
The script will output [BED6] (https://genome.ucsc.edu/FAQ/FAQformat#format1) with the following columns:

* column 1: chromosome
* column 2: start of polyA site (0-based)
* column 3: end of polyA site
* column 4: ID of the read containing a polyA tail
* column 5: length of the polyA tail on read
* column 6: genomic strand of the read (inferred from the mapping of the read, i.e. reads where a polyA tail was detected at their 3’ end are assigned a ‘+’ genomic strand, whereas reads with a polyT tail at their 5’ end are deduced to originate from the ‘-’ strand.)

## Dependencies (CPAN)
Getopt::Long

Bio::DB::Fasta

## Example run (on a BAM file)

`$ samtools view $file.bam |samToPolyA.pl --minClipped=20 --minAcontent=0.9 --minUpMisPrimeAlength=10 --genomeFasta=hg38.fa - > ${file}_polyAsites.bed`


## Author
Julien Lagarde, CRG, Barcelona, contact julienlag@gmail.com
