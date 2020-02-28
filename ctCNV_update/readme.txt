BRCA-LR:
	1)用ctCNV流程对样本进行cnv的分析；
	2)基于ctCNV的结果文件，算z_score值。

cmd:

perl /lustre/rdi/user/licq/work/BRCA-LR/ctCNV_update/script/z_score.pl /lustre/rdi/user/licq/work/BRCA-LR/ctCNV_update/demo/panel12_Baseline_BC_mkdup_ref_COV.txt /lustre/rdi/user/licq/work/BRCA-LR/ctCNV_update/demo/analysis_dir/  /lustre/rdi/user/licq/work/BRCA-LR/ctCNV_update/demo/Result 

perl /lustre/rdi/user/licq/work/BRCA-LR/ctCNV_update/script/z_score.pl + 基线文件(ctCNV baselines) + ctCNV的结果文件 + 输出目录 + 基线文件是否有log2ND的z_score值(有:值为1,没有:值为0)
