for (package in c('optparse')) {
  if (!require(package, character.only=T, quietly=T)) {
    install.packages(package,repos = "http://cran.us.r-project.org")
    library(package, character.only=T)
  }
}

option_list = list(
  make_option(c("-i", "--input"), type="character", default=NULL, 
              help="input_file", metavar="character"),
  make_option(c("-s", "--sample"), type="character", default=NULL, 
              help="sample", metavar="character"),
  make_option(c("-o", "--output"), type="character", default=NULL, 
              help="output_file", metavar="character")
)

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

var_file <- read.csv(opt$input)

#no_na <- var_file[which(var_file$cosmic68 != "." | var_file$X1000g2015aug_eur != "." | var_file$Func.refGene == "upstream"),]

## swapped cosmic68 for cosmic70
no_na <- var_file[which(var_file$cosmic70 != "." | var_file$X1000g2015aug_eur != "." | var_file$Func.refGene == "upstream"),]

#cosmic <- no_na[which(no_na$cosmic68 != "." | no_na$X1000g2015aug_eur < 0.01 | no_na$Func.refGene == "upstream"),]
## swapped cosmic68 for cosmic70
cosmic <- no_na[which(no_na$cosmic70 != "." | no_na$X1000g2015aug_eur < 0.01 | no_na$Func.refGene == "upstream"),]

no_syn <- cosmic[which(cosmic$ExonicFunc.refGene != "synonymous SNV"),]

no_syn <- no_syn[order(no_syn$X1000g2015aug_eur),]


table <- no_syn[which(no_syn$X1000g2015aug_eur < 0.01),c("Chr","Start","End","Ref","Alt","Func.refGene","Gene.refGene","ExonicFunc.refGene","AAChange.refGene","cytoBand","avsnp147","X1000g2015aug_eur","cosmic70" )]

var_list <- strsplit(table$AAChange.refGene[1],split = ":")
alteration <- lapply(table$AAChange.refGene, function(x){
  var_list <- strsplit(x,split = ":")
  alter <- paste(var_list[[1]][3],var_list[[1]][4],strsplit(var_list[[1]][5],split = ",")[[1]][1],sep = ";")
  return(alter)
})

cosmic <- lapply(table$cosmic70, function(x){
  var_list <- strsplit(x,split = ";")
  alter <- strsplit(var_list[[1]][1],split = ",")[[1]][1]
  return(alter)
})

filtered <- cbind("Chr"=table$Chr,"Start"=table$Start,"End"=table$End,"Func"=table$Func.refGene,"Gene"=table$Gene.refGene,"ExonicFunc"=table$ExonicFunc.refGene,"AAChange"=reshape2::melt(alteration)["value"],"cytoBand"=table$cytoBand,"1000g_EUR"=table$X1000g2015aug_eur,"COSMIC"=reshape2::melt(cosmic)["value"])

colnames(filtered) <- c("Chr","Start","End","Func","Gene","ExonicFunc","AAChange","cytoBand","1000g_EUR","COSMIC")


write.csv(filtered,file = opt$output,row.names = F)
