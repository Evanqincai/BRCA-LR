#! /use/bin/perl

use strict;
use warnings;
use File::Basename qw(basename dirname);
#use Statistics::Descriptive;
####################################### Set parameter
if (@ARGV != 6){
	print "perl calculate_zscore.pl bamfile panel12.bed  baselinesfile  sampleID  zscore panel2_BRCA_2.bed \n";
	die;
}

my $bam = $ARGV[0];
my $panel = $ARGV[1];
my $baselines = $ARGV[2];
my $sampleID = $ARGV[3];
my $outdir = $ARGV[4];
my $panel2 = $ARGV[5];
`mkdir $outdir` unless (-d $outdir);
`mkdir $outdir/work_sh` unless (-d "$outdir/work_sh");
######################################## calculate depth

open (CMD,">$outdir/work_sh/$sampleID.stat_depth.sh") or die $!;
#print CMD  "/usr/bin/bedtools  multicov -bams  $bam -bed $panel >$outdir/$sampleID.depth.txt";
print CMD "/usr/bin/bedtools map -a $panel -b $bam -c 10,10 -o count,concat|awk -v OFS=\"\t\" '{n=length(\$6); gc=gsub(\"[gcGC]\", \"\", \$6); print \$1\"\\t\"\$2\"\\t\"\$3\"\\t\"\$4\"\\t\"\$5\"\\t\"gc/n}' >$outdir/$sampleID.depth.txt\n";
close CMD;
`perl /lustre/rdi/user/licq/tools/qsub/qsub-sge.pl -maxjob 50 --queue rd.q --resource  mf=5G $outdir/work_sh/$sampleID.stat_depth.sh`;
###################################### baselines
my (%base_mean,%base_std,%base_up,%base_down);
open (BASE,$baselines) or die $!;
while (<BASE>){
	chomp;next if (/^#/ || /^$/);
	next if ($. == 1);
	my @lines = split(/\s+/,$_);
	my $index = "$lines[0]_$lines[1]_$lines[2]_$lines[3]";
	$base_mean{$index} = $lines[-6];
	$base_std{$index} = $lines[-5];
	$base_up{$index} = $lines[-2];
	$base_down{$index} = $lines[-1];
} 
close BASE;


####################################### normalization
my (%depth,@total_depth,%GC_sample,%GC_value,%GC_aver,%depth_correct,@total_depth_correct);
open (DEPTH,"$outdir/$sampleID.depth.txt") or die $!;
while (<DEPTH>){
	chomp;next if (/^$/ || /^#/);
	my @lines = split(/\s+/,$_);
	push @total_depth,$lines[4];
	my $gc = sprintf "%0.2f",$lines[-1];
	push @{$GC_sample{$gc}},$lines[4];
	my $index = "$lines[0]_$lines[1]_$lines[2]_$lines[3]";
	$depth{$index}{$sampleID} = $lines[4];
	$GC_value{$index} = $gc;
}
close DEPTH;
#my $stat = Statistics::Descriptive::Full->new();
#$stat->add_data(@total_depth);
#my $sample_mean = $stat -> mean();
my $sample_mean = &aver(\@total_depth);
foreach my $gc_value (keys %GC_sample){
	$GC_aver{$gc_value} = &aver_gc(\@{$GC_sample{$gc_value}},$sample_mean);
}
 

foreach my $pos_info_tmp (keys %depth){
        $depth_correct{$pos_info_tmp}{$sampleID} = $depth{$pos_info_tmp}{$sampleID}*$GC_aver{$GC_value{$pos_info_tmp}};
	push @total_depth_correct,$depth_correct{$pos_info_tmp}{$sampleID};
}
my $sample_mean_gc_correct = &aver(\@total_depth_correct);

###################################### z_score
my ($status,$z_score);
open (OUT,">$outdir/$sampleID.z_score.tmp.txt") or die $!;
foreach my $pos_info (keys %depth){
	my @info = split(/_/,$pos_info);
	my $info_tmp = join("\t",@info);
	my $ND = 0;
	if ($base_mean{$pos_info} && $base_std{$pos_info} && $base_up{$pos_info} && $base_down{$pos_info}){
		if ($base_std{$pos_info} != 0){
			$ND = $depth_correct{$pos_info}{$sampleID}/$sample_mean_gc_correct;;
			$z_score = ($ND - $base_mean{$pos_info})/$base_std{$pos_info};
		}else{
			$z_score = 0;
		}
		if ($z_score < $base_down{$pos_info}){
			if ($z_score == 0){
				$status = "NO Read";
			}else{
				$status = "Del";
			}
		}elsif ($z_score > $base_up{$pos_info}){
			$status = "Dup";
		}else{
			$status = "no cnv";
		}
		print OUT "$info_tmp\t$sample_mean_gc_correct\t$depth_correct{$pos_info}{$sampleID}\t$ND\t$base_up{$pos_info}\t$base_down{$pos_info}\t$z_score\t$status\n";
	}
	
}
close OUT;
`sort -k 1,1 -V -k 2n $outdir/$sampleID.z_score.tmp.txt | awk 'BEGIN{print \"chr\tstart\tend\texon\tmean_depth\tdepth\tND_depth\tmean2std_up\tmean2std_down\tZ_Score\tStatus\"}{print \$0}' >$outdir/$sampleID.z_score.txt && rm $outdir/$sampleID.z_score.tmp.txt`;

##################################### PNG
`awk 'NR==FNR{a[\$1\"\\t\"\$2\"\\t\"\$3]=1}NR!=FNR{if(a[\$1\"\\t\"\$2\"\\t\"\$3]){print \$0}}' $panel2 $outdir/$sampleID.z_score.txt  >$outdir/$sampleID.Z_score.braca.txt`;
`/lustre/rde/user/rde_admin/bin/Rscript /lustre/rdi/user/licq/project/BRCA1_BRCA2/Indel/CNV_script/z_score.dis.R  -s $sampleID -z $outdir/$sampleID.Z_score.braca.txt -o $outdir`;
##################################### Sub

sub aver {
        my $arr = shift;
        my $s = 0;
        grep {$s +=$_}@$arr;
        return ($s/@$arr);
}
sub aver_gc {
        my $arr = shift;
        my $sample_aver = shift;
        my $s = 0;
        grep {$s +=$_}@$arr;
        my $gc_aver = $s/@$arr;
        if ($gc_aver == 0){
		return 0;
	}else{
		return ($sample_aver/$gc_aver);
	}
}

