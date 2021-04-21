#############################################################################################################
# @Author: Ayush T. Raman
# Rai Lab, MDACC
# Date: August 30th, 2016
#
# Program is used for:
# 1. RNA-Seq Alignment
#
# Comments:
# 1. MPNST project
#############################################################################################################

#!/tools/bin/perl
use strict;
use warnings;

## Arguments
my $file = $ARGV[0];
my $result = $ARGV[1];
my $logs_dir = $ARGV[2];

## dat folder
open my $FILE, $file or die $!;
while(my $line = <$FILE>){
  chomp $line;
  my $read1 = "araman/2016-08-27_MPNST/".$line;
  my $read2 = "araman/2016-08-27_MPNST/".$line;
  $read2 =~ s/R1/R2/;
  my $gtf = "araman/2016-08-27_MPNST/STAR_index/genes.gtf";
  my $genome_dir = "araman/2016-08-27_MPNST/STAR_index/index";

  ## dir mkdir
  my @dir = split("\/",$line);
  my $dir1 = $result."/".$dir[2];
  my $dir2 = $dir1."/".$dir[3];
  system("mkdir $dir1") if(!-d $dir1);
  system("mkdir $dir2") if(!-d $dir2);
  my $out_filename = $dir2."/".$dir[3];

  ## read1 and read2 files for RNA-Seq
  my $cmd = "/nfs/nfs_ayush/aayushraman/software/STAR_2.4.2a/bin/Linux_x86_64/STAR --runThreadN 8 --genomeDir $genome_dir --readFilesCommand zcat --readFilesIn $read1 $read2 --sjdbGTFfile $gtf --sjdbOverhang 75 --outFileNamePrefix $out_filename --outSAMtype BAM SortedByCoordinate --quantMode GeneCounts";
  my $runName = $dir[3];
  print $runName,"\n\n",$cmd,"\n\n";
  parallelRun($runName, $cmd);
}
close($FILE);

## Run parallel
sub parallelRun{
	## Arguments for PBS file
	my @args = @_;
	my $runName = $args[0];
	my $emailID = "araman\@mdanderson.org";
	my $pbsFile = $logs_dir."/".$runName.".pbs";

	## PBS File
	open my $OUT, ">", $pbsFile or dies $!;
	print $OUT "#PBS -N $runName\n";
	print $OUT "#PBS -l nodes=1:ppn=8,walltime=8:00:00\n";
	print $OUT "#PBS -l mem=64gb\n";
	print $OUT "#PBS -M $emailID\n";
	print $OUT "#PBS -m n\n";
  print $OUT "#PBS –d	/scratch/genomic_med/araman/2016-08-27_MPNST/ \n";
	print $OUT "#PBS -o $logs_dir\n";
	print $OUT "#PBS -e $logs_dir\n";
	print $OUT "##\n\n";
	print $OUT "echo \"My job ran on: \$(date)\"\n";
	print $OUT "cat \$PBS_NODEFILE\n";
	print $OUT "$args[1]\n";
  print $OUT "##\n\n";
	close($OUT);

	## running the pbs file
	my $jobID = `msub $pbsFile`;
	$jobID =~ s/^\s+|\s+$//g;
	print "Job ID = ",$jobID," for sample name:",$pbsFile,"\n";
}
