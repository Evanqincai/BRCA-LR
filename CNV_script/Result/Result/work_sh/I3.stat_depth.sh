/usr/bin/bedtools map -a /lustre/rdi/user/licq/project/BRCA1_BRCA2/Indel/CNV_script/panel12.bed -b /beegfs/work/commercial_test/pipeline_V2/181114_A00168_0230_AH7NN5DSXX.20181116123002/R-181110-535748-BLD-521005_BC_E18111112-3-DNA_F1811112948-16_L181111-00002-A5_20181113-C12_panel12_P181114419-3_sample20181114-A3-A00168.20181116123002/basic_analysis/align/R-181110-535748-BLD-521005_BC_E18111112-3-DNA_F1811112948-16_L181111-00002-A5_20181113-C12_panel12_P181114419-3_sample20181114-A3-A00168.sorted.rmdup.realign.bam -c 10,10 -o count,concat|awk -v OFS="	" '{n=length($6); gc=gsub("[gcGC]", "", $6); print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"gc/n}' >Result/I3.depth.txt
