library(getopt)
opt_spec <- matrix(c('help', 'h', 0, 'logical', 'help manual',
                     'sampleid', 's', 1, 'character', 'id of the sample',
                     'z_score', 'z', 1, 'character', 'tomor sample z_score.txt',
                     'output', 'o', 1, 'character', 'output directory'),
                   byrow=TRUE, ncol=5)
opt = getopt(opt_spec, commandArgs(TRUE))
if (!is.null(opt$help)){
  cat(getopt(opt_spec, usage=TRUE))
  q(save='no', status=1)
}
if ( is.null(opt$z_score) )  { cat(getopt(opt_spec,usage=TRUE));q(status=1)}
if ( is.null(opt$sampleid))  { cat(getopt(opt_spec,usage=TRUE));q(status=1)}
if ( is.null(opt$output))    { cat(getopt(opt_spec,usage=TRUE));q(status=1)}
#stopifnot(file.exists(opt$snv))
#----------------------------------------------------------------  z_score
library(ggplot2)
data <- read.table(opt$z_score,header = T,sep = "\t")
p<- ggplot(data,aes(start/1000000,Z_Score)) + geom_point(aes(shape = Status,colour = Status)) + facet_wrap(~chr,scales = "free_x") + geom_line(aes(x = data$start/1000000 , y = data$mean2std_up),color = "red",linetype = 1) + geom_path(aes(x = data$start/1000000 , y = data$mean2std_down),color = "red",linetype = 1) + xlab("Pos(M)")  
outfile = paste(opt$output,"/",opt$sample,".z_score.pdf",sep = "")
#----------------------------------------------------------------- png
pdf(file=outfile)
print(p)
dev.off()

#--------------------------------------------------------------- ND_depth
p1<- ggplot(data,aes(factor(start/1000000),ND_depth)) + geom_bar(stat = "identity") + facet_wrap(~chr,scales = "free_x") + xlab("Pos(M)") + theme(axis.text.x = element_text(angle = 90, hjust = 1))
#---------------------------------------------------------------- png
outfile = paste(opt$output,"/",opt$sample,".DN_depth.pdf",sep = "")
pdf(file=outfile)
print(p1)
dev.off()
