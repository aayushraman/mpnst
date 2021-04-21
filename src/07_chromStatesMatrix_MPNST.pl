################################################################################
# @Author: Ayush T. Raman
# Date: May 09, 2017
#
# Program is used for:
# 1. Chrom States Matrix for 2 MB/100 Kb/10 Kb region
################################################################################

#!/tools/bin/perl
use strict;
use warnings;

print "The arguments should be given in this order: Number_of_States_used_in_calling_ChromHMM State_Number_user_is_interested_in Bin_Size Sample_Name_File ChromHMM_Folder \n";

## Input from user
my $LearnStateNumber = $ARGV[0];
my $STATE_NUMBER = "E".$ARGV[1];
my $bin = $ARGV[2];
my $sampleNameFile = $ARGV[3];
my $ChromHMM_folder = $ARGV[4];

## Variables
my @chr = ("chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7", "chr8", "chr9", "chr10", "chr11", "chr12", "chr13", "chr14", "chr15", "chr16", "chr17", "chr18", "chr19", "chr20", "chr21", "chr22", "chrX");
my @cellLines;
my %chromSizeHash;
my %chromGC; ## chromosome genomic coordinates
my %chromStateMatrix;
my $count_gc = 0;
my $start = time;

## Reading sampleName Files
open my $FILE, $sampleNameFile or die "Cannot open file with sample names $! \n"; ## $ChromHMM_folder."/".
while(my $line = <$FILE>){
	chomp $line;
	my @cols = split("\t", $line);
	push @cellLines, $cols[0];
	print $cols[0],"\n";
}
close($FILE);

## Working directory, Output File for Chrom State Matrix and Chromosome Length File
open my $OUT, ">", $ChromHMM_folder."/CombinedMatrix-".$STATE_NUMBER."-".$bin."bps.txt" or die "Output Folder does not exists $! \n";

## Chromosome Length File
open my $IN, "../dat/ChromSize_hg19.txt" or die "Chromosome Size File does not exists $! \n";
while(my $line = <$IN>){
	chomp $line;
	my @cols = split("\t", $line);

	## Defining Genomic Coordinates of region
	$chromSizeHash{$cols[0]} = $cols[1];
	for(my $i = 0; $i <= $cols[1]; $i+=$bin){
		my $start = $i;
		my $end = $i+$bin-1;

	# if the $end is bigger than the size of the chromosome
		$end = $cols[1] if($end >= $cols[1]);
		my $gc = $cols[0].":".$start."-".$end;
		$chromGC{$cols[0]}{$gc} = 1;
		$count_gc++;
	}
}
close($IN);
print "Number of genomic bins in hg19 for binsize $bin basepairs are = ",$count_gc,"\n";

## Reading the segment bed files
foreach my $cell (@cellLines){
	print "Running: ",$cell."_".$LearnStateNumber."_segments.bed file \n";
	my $bedFile = $ChromHMM_folder."/".$cell."_".$LearnStateNumber."_segments.bed";
	open my $segment_bed, $bedFile or die "Cannot open the file $!" ;

	## Reading the Segment Bed file
	while(my $line = <$segment_bed>){
		chomp $line;
		my @cols = split("\t", $line);
		my $chrBed = $cols[0];
		my $startBed = $cols[1];
		my $endBed = $cols[2];
		my $chromStateBed = $cols[3];
		my $chromSize = $chromSizeHash{$chrBed};
		my $lengthBed = $startBed - $endBed;

		## Defining the cols of the bed file
		if(exists $chromGC{$chrBed} && $chromStateBed eq $STATE_NUMBER){

			## genomic cord as possible bin
			my $genomic_StartBin = int($startBed/$bin);
			my $genomic_EndBin = int($endBed/$bin);

			## check if the $binNumber_start and $binNumber_end are in the same bin or not
			if($genomic_StartBin == $genomic_EndBin){
				my $start =  $genomic_StartBin * $bin;
				my $end = $start + $bin - 1;
				$end = $chromSize if($end > $chromSize);
				my $posBin = $chrBed.":".$start."-".$end;
				$chromStateMatrix{$chrBed}{$posBin}{$chromStateBed}{$cell} += 1;
			}
			elsif($genomic_EndBin > $genomic_StartBin){
				for(my $i = $genomic_StartBin; $i<= $genomic_EndBin; $i++){
					my $start =  $i * $bin;
					my $end = $start + $bin - 1;
					$end = $chromSize if($end > $chromSize);
					my $posBin = $chrBed.":".$start."-".$end;
					$chromStateMatrix{$chrBed}{$posBin}{$chromStateBed}{$cell} += 1;
				}
			}
		}
	}
	close($segment_bed);
	print "Calculation of the number of states in each bin of $bin for ".$cell."_".$LearnStateNumber."_segments.bed"." is done \n";
	my $duration = time - $start;
	print "Execution time to complete the run for $cell: $duration s\n"; sleep(2);
  print "\n\n";
}

## Output header
print $OUT "Genomic_Loci";
foreach my $cell (@cellLines){
  print $OUT "\t",$cell,"-",$STATE_NUMBER;
}
print $OUT "\n";

## Forming the Chrom State Matrix
foreach my $chr (@chr){
  foreach my $posBin (keys $chromStateMatrix{$chr}){
		print $OUT $posBin;
		foreach my $cell (@cellLines){
			print $OUT "\t",$chromStateMatrix{$chr}{$posBin}{$STATE_NUMBER}{$cell}
            if(exists $chromStateMatrix{$chr}{$posBin}{$STATE_NUMBER}{$cell});
			print $OUT "\t0"
            if(!exists $chromStateMatrix{$chr}{$posBin}{$STATE_NUMBER}{$cell});
		}
		print $OUT "\n";
	}
}
close($OUT);

## Total Time taken to run the script
print "Execution time to complete the entire run: ",time - $start,"s\n";
