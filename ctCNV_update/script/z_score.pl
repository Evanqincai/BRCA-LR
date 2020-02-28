
use strict;
use warnings;
use File::Basename qw(basename dirname);

if (@ARGV != 4){
	print "perl z_score.pl ref_cov analysis_dir outdir flag\n";
	die();
}

#####################################################################
`mkdir $ARGV[2]` unless (-d $ARGV[2]);
my %hash_stat;
if ($ARGV[3] == 1){
	(%hash_stat) = &save_hash($ARGV[0]);
}else{
	my $update_file = &score($ARGV[0],$ARGV[2]);
	(%hash_stat) = &hash($update_file);
}
my @RZfile = glob("$ARGV[1]/*RZ.txt");
for (my $i = 0 ;$i <@RZfile ; $i++) {
	my $filename = (split(/.RZ/,basename($RZfile[$i])))[0];
	open (IN,"paste $ARGV[1]/${filename}_COV.txt $RZfile[$i]|cut -f 1-4,8-|");
	open (OUT,">$ARGV[2]/$filename.RZ.update.txt") or die $!;
	while (<IN>){
		chomp;
		my @lines = split(/\s+/,$_);
		if ($lines[0] eq "Chr"){
			print OUT "$_\tup\tdown\tStatus\n";
		}else{
			my $index = "$lines[0]_$lines[1]_$lines[2]_$lines[3]";
			my ($up,$down) = split("_",$hash_stat{$index}) if ($hash_stat{$index});
			if ($up && $down){
				if ($lines[8] > $up){
					print OUT "$_\t$up\t$down\tDup\n";
				}elsif ($lines[8] < $down){
					print OUT "$_\t$up\t$down\tDel\n";
				}else{
					print OUT "$_\t$up\t$down\tno cnv\n";
				}
			}
		}
	}
	close IN;
	close OUT;
	`less $ARGV[2]/$filename.RZ.update.txt |grep -w -E 'BRCA1|BRCA2|Chr' >$ARGV[2]/$filename.RZ.BRCA.update.txt`;
	`/lustre/rde/user/rde_admin/bin/Rscript /lustre/rdi/user/licq/project/BRCA1_BRCA2/Indel/ctCNV_update/z_score.dis.R  -s $filename -z $ARGV[2]/$filename.RZ.BRCA.update.txt -o $ARGV[2]`;
}
######################################################################
sub hash {
	my $infile = shift;
	my (%hash_stat);
	open (IN,$infile) or die $!;
	while (<IN>){
		chomp;
		next if (/^$/ || /^#/ || /^Chr/);
		my @lines = split(/\s+/,$_);
		my $index = "$lines[0]_$lines[1]_$lines[2]_$lines[3]";
		$hash_stat{$index} = "$lines[-2]_$lines[-1]";
	}
	close IN;
	return (%hash_stat);
}

sub score {
	my $file = shift;
	my $outdir = shift;
	my $POS;
	my $filename = (split(/.txt/,basename($file)))[0];
	open (IN,$file) or die $!;
	open (OUT,">$outdir/$filename.update.txt") or die $!;
	while (<IN>){
		chomp;
		my ($score_value,$aver,$std,@SCORE,@title,$mean_plus_2std,$mean_sub_2std);
		my @lines = split (/\s+/,$_);
		if ($lines[0] =~/Chr/){
			for (my $i= 0 ;$i <@lines-3; $i++){
				if ($lines[$i] =~ /log2_ND/){
					if (!$POS){
						$POS = $i;
					}
					push @title,"$lines[$i]_SCORE";
				}
			}
			my $title = join("\t",@title);
			print OUT "$_\t$title\tmean_score\tstd_score\tmean_plus_2std\tmean_sub_2std\n";	
		}else{
			for (my $i=$POS ;$i <@lines-3;$i++){
				push @SCORE,($lines[$i]-$lines[-3])/$lines[-2];
			}
			$aver = &aver(\@SCORE);
			$std = &var(\@SCORE);
			$score_value = join ("\t",@SCORE);
			$mean_plus_2std = $aver + 2*$std;
			$mean_sub_2std = $aver - 2*$std;
			print OUT "$_\t$score_value\t$aver\t$std\t$mean_plus_2std\t$mean_sub_2std\n";
		}
	}
	close IN;
	close OUT;
	return ("$outdir/$filename.update.txt");

}


sub aver {
        my $arr = shift;
        my $s = 0;
        grep {$s +=$_}@$arr;
        return ($s/@$arr);
}


sub var {
        my $arr = shift;
        my $v = &aver($arr);
        my $d = 0;
        grep {$d += (($_ - $v)**2);}@$arr;
        my $std = sqrt($d/((@$arr)-1));
        return $std;

}
