/usr/bin/bedtools map  -a /lustre/rdi/user/licq/project/BRCA1_BRCA2/Indel/CNV_script/panel12.bed -b /lustre/rdi/user/licq/project/panel12_second/BC_Data/sample20180810-A2_capture7_I23.sorted.rmdup.realign.bam -c 10,10 -o count,concat|awk -v OFS="	" '{n=length($6); gc=gsub("[gcGC]", "", $6); print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"gc/n}' >/lustre/rdi/user/licq/project/BRCA1_BRCA2/Indel/CNV_script/script/baselines//I29.depth.txt
