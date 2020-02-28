校正方法:

    先统计每个GC含量（0, 1, 2, 3,…, 100%）下的捕获区间的平均覆盖度，再计算所有捕获区间的平均覆盖度，用来校正测序得到的覆盖度

    特定捕获区间校正后的覆盖度 = 该捕获区间的原始覆盖度 *（所有捕获区间的平均覆盖度/与该捕获区间的有相同GC含量的所有捕获区间的平均覆盖度）

命令行:
1)建立基线
   perl /lustre/rdi/user/licq/work/BRCA-LR//CNV_script/calculate_baselines.pl /lustre/rdi/user/licq/work/BRCA-LR//CNV_script/Result/baselines.contig  /lustre/rdi/user/licq/work/BRCA-LR//CNV_script/panel12.bed /lustre/rdi/user/licq/work/BRCA-LR//CNV_script/Result/baselines/

  第一个参数：配置文件,即建立基线样本的Bam文件
  第二个参数：样本对应的bed文件
  第三个参数：输出目录

2)Z-score模型

  perl /lustre/rdi/user/licq/work/BRCA-LR//CNV_script/calculate_zscore.pl /beegfs/work/commercial_test/pipeline_V2/181114_A00168_0230_AH7NN5DSXX.20181116123002/R-181108-995555-BLD-101100_BC_E18111102-16-DNA_F1811112948-6_L181111-00002-A4_20181113-C11_panel12_P181114419-3_sample20181114-A3-A00168.20181116123002/basic_analysis/align/R-181108-995555-BLD-101100_BC_E18111102-16-DNA_F1811112948-6_L181111-00002-A4_20181113-C11_panel12_P181114419-3_sample20181114-A3-A00168.sorted.rmdup.realign.bam /lustre/rdi/user/licq/work/BRCA-LR//CNV_script/panel12.bed baselines/baselines.txt I1 Result /lustre/rdi/user/licq/work/BRCA-LR//CNV_script/panel12.BRCA.bed

   第一个参数：待测样本的bam文件
   第二个参数：待测样本对应的panel的bed文件
   第三个参数：基线文件
   第四个参数：待测样本ID
   第五个参数：输出目录
   第六个参数：BRCA1/2基因对应的bed文件
