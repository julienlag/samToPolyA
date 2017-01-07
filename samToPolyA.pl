#!/usr/bin/env perl

use strict;
use warnings;
use Getopt::Long;
use Bio::DB::Fasta;

$|=1;

###### this script maps polyA sites on the genome based on read mappings in SAM format, and according to following provided parameters:
#
#
##  minClipped: integer
#              = minimum length of A or T tail required to call a PolyA site
##  minAcontent: float
#              = required A (or T, if minus strand) content of the tail
#              Note: minAcontent affects both the A tail and the upstream A stretch.
##  minUpMisPrimeAlength: integer
#              = minimum length of genomic A stretch immediately upstream a putative site required to call a false positive (presumably due to internal RT priming), and hence not report the corresponding site
##  genomeFasta: string
#              = path to multifasta of genome (+ spike-in sequences if applicable), used to extract upstream genomic sequence.
#
#
## The script will output BED6 with the following columns:
# col1 -> chr
# col2 -> start of polyA site (0-based)
# col3 -> end of polyA site
# col4 -> ID of the read containing a polyA tail
# col5 -> length of the polyA tail on read
# col6 -> strand of the read (inferred from the mapping of the read, i.e. '-' if the detected A tail is at the beginning of the read (polyT tail), and '+' if it's at the end of it (polyA tail).)
## author: Julien Lagarde, CRG, Barcelona, contact julienlag@gmail.com

my $minSoftClippedSeqLength=10;
my $minAcontent=0.8;
my $minUpMisPrimeAlength=10;
my $genomeFa;
GetOptions ('minClipped=i' => \$minSoftClippedSeqLength,
            'minAcontent:f' => \$minAcontent,
            'minUpMisPrimeAlength=i' => \$minUpMisPrimeAlength,
            'genomeFasta=s' => \$genomeFa)
  or die("Error in command line arguments\n");

my $chrdb = Bio::DB::Fasta->new("$genomeFa", -reindex => 0);



while (<STDIN>){
	my $line=$_;
	chomp;
	next if($_=~/^\@SQ/); #skip sequence header
	my @line=split "\t";
	next if($line[5] eq '*'); #skip unmapped reads
	my @cigarNumbers=split (/[A-Z]/,$line[5]);
	my @cigarLetters=split(/\d+/,$line[5]);
	shift(@cigarLetters); #the first element is empty, remove it
	#check if there's any soft-clipped sequence at the end(s) of the read (i.e. unmapped tail)
	if($cigarLetters[0] eq 'S' && $cigarNumbers[0] > $minSoftClippedSeqLength){ #tail is likely at the beginning of the read
	#it means strand is -
		my $tailSeq=substr($line[9], 0, $cigarNumbers[0]);
		#reverse string so that we start scanning the sequence where it drops from the genome
		$tailSeq=reverse($tailSeq);
		$tailSeq=~ tr/ACGTacgt/TGCAtgca/;
		my $strand='-';
		#there might be an exogenous sequence adapter at the end of the sequence.
		my $aTailLength=countAs($tailSeq, $minSoftClippedSeqLength);
		if ($aTailLength>0){
			my $upstreamStop=$line[3]+$minUpMisPrimeAlength;
			my $upstreamGenomeSeq=$chrdb->seq($line[2], $line[3], $upstreamStop);
			$upstreamGenomeSeq=reverse($upstreamGenomeSeq);
			$upstreamGenomeSeq=~ tr/ACGTacgt/TGCAtgca/;
			my $upstreamGenomeAsLength=countAs($upstreamGenomeSeq, $minUpMisPrimeAlength);
			unless ($upstreamGenomeAsLength >= $minUpMisPrimeAlength){
				my $start=$line[3]-1;
				print "$line[2]\t$start\t$line[3]\t$line[0]\t$aTailLength\t$strand\n";
				next; # found, no need to look at the other end of the read
			}
			}
	}

	if($cigarLetters[$#cigarLetters] eq 'S' && $cigarNumbers[$#cigarNumbers] > $minSoftClippedSeqLength){ #tail is likely at the end of the read
	# this is no ELSIF!! need to check second end if first fails
		my $tailSeq=substr($line[9], length($line[9]) - $cigarNumbers[$#cigarNumbers], $cigarNumbers[$#cigarNumbers]);
		my $strand='+';
		#there might be an exogenous sequence adapter at the end of the sequence.
		my $aTailLength=countAs($tailSeq, $minSoftClippedSeqLength);
		if($aTailLength>0){
			#calculate where the start of the tail is on the genome
			my $genomicLength=0;
			for (my $i=0; $i<=$#cigarLetters;$i++){
				if($cigarLetters[$i] =~ /[MDNXP]/){
					$genomicLength+=$cigarNumbers[$i];
				}
			}
			my $start=$line[3]+$genomicLength-1;
			my $end=$start+1;
			my $upstreamStop=$start;
			my $upstreamStart=$start-$minUpMisPrimeAlength;
			my $upstreamGenomeSeq=$chrdb->seq($line[2], $upstreamStart, $upstreamStop);
			my $upstreamGenomeAsLength=countAs($upstreamGenomeSeq, $minUpMisPrimeAlength);
			unless ($upstreamGenomeAsLength >= $minUpMisPrimeAlength){
				print "$line[2]\t$start\t$end\t$line[0]\t$aTailLength\t$strand\n";
			}
		}
	}
}

sub countAs{
	my $minSeqLength=$_[1];
	my $seq=lc($_[0]);
	my @seq=split("", $seq);
	my $countAs=0;
	my $maxMismatches=$minSeqLength-($minAcontent*$minSeqLength);
	my $SeqLength=0;
	my $countMismatches=0;
	for(my $i=0; $i<=$#seq;$i++){
	 my $nt=$seq[$i];
	 $SeqLength=$i+1;
	 if($countMismatches>$maxMismatches && $SeqLength>$minSeqLength){
		 last;
	 }
	 if($nt eq "a"){
	 	$countAs++;
	 }
	 else{
	 	$countMismatches++;
	 }
	}
	if($SeqLength >= $minSeqLength && $countAs > $minSeqLength-$maxMismatches){# && $countAs/$SeqLength >= $minAcontent){

		return $SeqLength;
	}
	else{
		return 0;
	}

}
