for (package in c('optparse')) {
  if (!require(package, character.only=T, quietly=T)) {
    install.packages(package,repos = "http://cran.us.r-project.org")
    library(package, character.only=T)
  }
}

option_list = list(
  make_option(c("-i", "--input"), type="character", default=NULL, 
              help="mgmt bed file", metavar="character"),
  make_option(c("-p", "--probes"), type="character", default=NULL, 
              help="top probes", metavar="character"),
  make_option(c("-m", "--model"), type="character", default=NULL, 
              help="mgmt prediction model", metavar="character"),
  make_option(c("-o", "--out_dir"), type="character", default=NULL, 
              help="output directory", metavar="character"),
  make_option(c("-s", "--sample"), type="character", default=NULL, 
              help="sample", metavar="character"),
  make_option(c("-t", "--threads"), type="character", default=NULL,
              help="threads", metavar="character")
)

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

load(opt$probes)
load(opt$model)
mgmt_meth <- read.delim(opt$input,header = F)
mgmt_meth <- subset(mgmt_meth, V2 %in% pred_pos)
## temp file added to check effect of na.rm added below
#write.csv(mgmt_meth,file=paste0(opt$out_dir,"/",opt$sample,"_mgmt_meth_to_check.csv"),row.names = T)
## na.rm=T added by GF
#mgmt_average <- mean(mgmt_meth$V11, na.rm=T)
mgmt_average <- mean(t(as.numeric(as.data.frame(strsplit(mgmt_meth$V10, ' ', fixed=TRUE))[2,])), na.rm=T)
#mgmt_average <- mean(mgmt_meth$V11)
mgmt <- data.frame(average=mgmt_average)
#pred <- predict(log.model,newdata = mgmt,type = "response")
## GF added a cores argument
pred <- predict(log.model,newdata = mgmt,type = "response", cores=opt$threads)
mgmt$pred <- pred[[1]]
if (pred[[1]] <0.5){mgmt$status <- "Unmethylated"} else {mgmt$status <- "Methylated"}
write.csv(mgmt,file=paste0(opt$out_dir,"/",opt$sample,"_mgmt_status.csv"),row.names = F)
