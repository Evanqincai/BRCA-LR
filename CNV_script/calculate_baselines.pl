#! /use/bin/perl

use strict;
use warnings;
use File::Basename qw(basename dirname);
#use Statistics::Descriptive;
use Data::Dumper;
####################################### Set parameter
if (@ARGV != 3){
	print "perl calculate_baselines.pl  <config>  <panel18.bed> <baselines>\n";
	die;
}
my $bamconfig = $ARGV[0];
my $panel = $ARGV[1];
my $outdir = $ARGV[2];
`mkdir $outdir` unless (-d $outdir);
`mkdir $outdir/work_sh` unless (-d "$outdir/work_sh");
######################################## calculate depth

open (CONFIG,$bamconfig) or die $!;
while (<CONFIG>){
	chomp;
	next if (/^#/ || /^$/);
	my @lines = split (/\s+/,$_);
	open (CMD,">$outdir/work_sh/$lines[1].sh") or die $!;
	#print CMD  "/usr/bin/bedtools  multicov -bams  $lines[0] -bed $panel >$outdir/$lines[1].depth.txt\n";
	print CMD "/usr/bin/bedtools map  -a $panel -b $lines[0] -c 10,10 -o count,concat|awk -v OFS=\"\t\" '{n=length(\$6); gc=gsub(\"[gcGC]\", \"\", \$6); print \$1\"\\t\"\$2\"\\t\"\$3\"\\t\"\$4\"\\t\"\$5\"\\t\"gc/n}' >$outdir/$lines[1].depth.txt\n";
	close CMD;
	`perl /lustre/rdi/user/licq/tools/qsub/qsub-sge.pl -maxjob 50 --queue rd.q --resource  mf=5G  $outdir/work_sh/$lines[1].sh`;
}
close CONFIG;

#`sh $outdir/work_sh/stat_depth.sh`;
####################################### normalization and GC Correct
my @depth = glob("$outdir/*depth.txt");
my (%depth_ND,%depth_ND_Sample,%depth,@sample_list,%depth_correct);
for (my $num = 0 ;$num <@depth ;$num++){
	my (@total_depth,%GC_sample,%GC_aver,%GC_value,@total_depth_correct);
	my $filename = (split(/\./,basename($depth[$num])))[0];
	push @sample_list,$filename;
	open (DEPTH,$depth[$num]) or die $!;
	while (<DEPTH>){
		chomp;next if (/^$/ || /^#/);
		my @lines = split(/\s+/,$_);
		push @total_depth,$lines[4];
		my $gc = sprintf "%0.2f",$lines[-1];
		push @{$GC_sample{$gc}},$lines[4];
		my $index = "$lines[0]_$lines[1]_$lines[2]_$lines[3]";
		$depth{$index}{$filename} = $lines[4];
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
		$depth_correct{$pos_info_tmp}{$filename} = $depth{$pos_info_tmp}{$filename}*$GC_aver{$GC_value{$pos_info_tmp}};
		push @total_depth_correct,$depth_correct{$pos_info_tmp}{$filename};
	}
	my $sample_mean_gc_correct = &aver(\@total_depth_correct);

	foreach my $pos_info (keys %depth_correct){
		my @info = split(/_/,$pos_info);
		my $info_tmp = join("\t",@info);
		if ($depth_correct{$pos_info}{$filename} == 0){
			$depth_ND_Sample{$pos_info}{$filename} = 0;
		}else{
			$depth_ND_Sample{$pos_info}{$filename} = $depth_correct{$pos_info}{$filename}/$sample_mean_gc_correct;
		}
	}
	
}
my $sample_list_depth = join("\t",(map {"$_-Depth"}@sample_list));
my $sample_list_ND = join("\t",(map{"$_-ND"}@sample_list));
my $sample_list_Z_Score = join("\t",(map{"$_-Z_Score"}@sample_list));
####################################### mean and std
open (BASE,">$outdir/baselines.tmp.txt") or die $!;
foreach my $pos_info (keys %depth_ND_Sample){
	my @info = split(/_/,$pos_info);
	my $info_tmp = join("\t",@info);
	#my $stat = Statistics::Descriptive::Full->new();
	#$stat->add_data(@{$depth_ND{$pos_info}});
	#my $ba_mean = $stat -> mean();
	#my $ba_std = $stat -> variance();
	my @ND_value = values %{$depth_ND_Sample{$pos_info}};
	my $ND_mean = &aver(\@ND_value);
	my $ND_std =  &var(\@ND_value);
	my ($up,$down,$mean,$std,%z_score) = &quantile(\%{$depth_ND_Sample{$pos_info}},$ND_mean,$ND_std);
	print BASE "$info_tmp\t";
	for (my $i = 0 ;$i < @sample_list; $i++){ 
		print BASE "$depth{$pos_info}{$sample_list[$i]}\t";
	}
	for (my $j = 0;$j < @sample_list; $j++){
		print BASE "$depth_ND_Sample{$pos_info}{$sample_list[$j]}\t"
	}

	for (my $k = 0;$k < @sample_list; $k++){
                print BASE "$z_score{$sample_list[$k]}\t"
        }
	my $mean2std_up = $mean + 2 *$std;
	my $mean2std_down = $mean - 2*$std;
	print BASE "$ND_mean\t$ND_std\t$up\t$down\t$mean2std_up\t$mean2std_down\n";
}
close BASE;

`sort -k 1,1 -V -k 2n $outdir/baselines.tmp.txt | awk 'BEGIN{print \"chr\tstart\tend\texon\t$sample_list_depth\t$sample_list_ND\t$sample_list_Z_Score\tMean_ND\tStd_ND\tquantile95\tquantile05\tmean2std_up\tmean2std_down\"}{print \$0}' >$outdir/baselines.txt && rm $outdir/baselines.tmp.txt`;

#####################################  Sub
sub quantile {
	my ($count,$list,$ba_mean,$ba_std,@z_score,@z_score_sort,%z_score);
	($list,$ba_mean,$ba_std) = @_;
	foreach my $sampleID (keys %{$list}){
		if ($ba_std == 0){
			$z_score{$sampleID} = 0;
		}else{
			$z_score{$sampleID} = (${$list}{$sampleID} - $ba_mean)/$ba_std;
		}

	}
	@z_score = values %z_score;
	@z_score_sort = sort {$a <=> $b} @z_score;
	if (scalar(@z_score_sort) == 0){
		return 0;die;
	}
	my $up = $z_score_sort[int((scalar(@z_score_sort)-1)*0.95)];
	my $down = $z_score_sort[int((scalar(@z_score_sort) -1)*0.05)];
	my $mean = &aver(\@z_score_sort);
	my $std = &var(\@z_score_sort);
	#my $up = $list[round(($count-1)*0.95)];
	#my $down = $list[round(($count -1)*0.05)];
	#return ($up,$down);
	return ($up,$down,$mean,$std,%z_score);
}


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

sub var {
	my $arr = shift;
	my $v = &aver($arr);
	my $d = 0;
	grep {$d += (($_ - $v)**2);}@$arr;
	my $std = sqrt($d/((@$arr)-1));
	return $std;
	
}
