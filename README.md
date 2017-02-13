# NAME

samToPolyA

# SYNOPSIS

A utility to detect poly-adenylated sequencing reads, call on-genome polyA sites and infer the reads' strand based on reads-to-genome alignments in SAM format.

**Usage example** (on a BAM file):

`samtools view $file.bam |samToPolyA.pl --minClipped=20 --minAcontent=0.9 --minUpMisPrimeAlength=10 --genomeFasta=hg38.fa - > ${file}_polyAsites.bed`

## INPUT

Read-to-genome alignments in SAM format, and the corresponding genome sequence in multifasta.

The script looks for terminal soft-clipped A/T sequences (marked as "S" in the CIGAR string).

## OPTIONS

This script maps polyA sites on the genome based on read mappings in SAM format, and according to the following provided parameters:

- **minClipped** (integer) = minimum length of A or T tail required to call a PolyA site.

    Default: '10'.

- **minAcontent** (float) = required A (or T, if minus strand) content of the tail.

    Default: '0.8'.

    Note: minAcontent affects both the A tail and the upstream A stretch.

- **minUpMisPrimeAlength** (integer) = minimum length of genomic A stretch immediately upstream a putative site required to call a false positive (presumably due to internal RT priming), and hence not report the corresponding site in the output.

    Default: '10'.

- **genomeFasta** (string) = path to multifasta of genome (+ spike-in sequences if applicable), used to extract upstream genomic sequence.

## OUTPUT

The script will output BED6 with the following columns:

- column 1: chromosome
- column 2: start of polyA site (0-based)
- column 3: end of polyA site
- column 4: ID of the read containing a polyA tail
- column 5: length of the polyA tail on read
- column 6: genomic strand of the read (see DESCRIPTION below)

# DESCRIPTION

The script will search for read alignment patterns such as:

`XXXXXXXXXXXAAAAAAAAAAAAAAA(YYYY) [read]`

`|||||||||||..................... [match]`

`XXXXXXXXXXXZ-------------------- [reference sequence]`

or

`(YYYY)TTTTTTTTTTTTTTTTXXXXXXXXXX [read]`

`......................|||||||||| [match]`

`---------------------ZXXXXXXXXXX [reference sequence]`

Where:

- `|` / `.` = a position mapped / unmapped to the reference, respectively
- `X` = the mapped portion of the read or reference sequence
- `(Y)` = an optional soft-clipped, non-(A|T)-rich sequence (possibly a sequencing adapter)
- `Z` = the position on the reference sequence where the alignment breaks
- The `A` / `T` streches are soft-clipped ('S' in CIGAR nomenclature) in the alignment
- `-` = the portion of the reference sequence unaligned to the read

The genomic strand of the read + polyA site is inferred from the mapping of the read, _i.e._, reads where a polyA tail was detected at their 3' end are assigned a '+' genomic strand, whereas reads with a polyT tail at their 5' end are deduced to originate from the '-' strand. In that example, the first / second alignment would lead to a called polyA site at position Z on the '+' / '-' strand of the reference sequence, respectively.

# DEPENDENCIES

CPAN: Bio::DB::Fasta

# AUTHOR

Julien Lagarde, CRG, Barcelona, contact julienlag@gmail.com
