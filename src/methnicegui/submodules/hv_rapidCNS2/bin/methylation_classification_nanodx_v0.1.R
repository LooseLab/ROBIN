#Check and install packages

for (package in c('optparse', 'GenomicRanges','ranger','matrixStats','data.table', 'glmnet' )) {
  if (!require(package, character.only=T, quietly=T)) {
    install.packages(package,repos = "http://cran.us.r-project.org")
    library(package, character.only=T)
  }
}

#Parse arguments
option_list = list(
  make_option(c("-s", "--sample"), type="character", default=NULL, 
              help="sample", metavar="character"),
  make_option(c("-o", "--out_dir"), type="character", default=NULL, 
              help="output directory", metavar="character"),
  make_option(c("-i", "--in_file"), type="character", default=NULL,
              help="path to modified base called file", metavar="character"),
  make_option(c("-p", "--probes_file"), type="character", default="top_probes_hm450.Rdata"),
  make_option(c("-a", "--array_file"), type="character", default="top_probes_hm450.Rdata"),
  make_option(c("-d", "--training_data"), type="character", default="capper_betas.RData",
              help="capper dataset", metavar="character"),
  make_option(c("-t", "--threads"), type="numeric", default=16,
              help="number of threads", metavar="character")
)

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

#Create output directory
dir.create(file.path(opt$out_dir), showWarnings = FALSE)

#Read megalodon methylation values
#filein <- paste0(opt$megalodon_dir,"/modified_bases.5mC.bed")
meth <- read.delim(opt$in_file,header=FALSE)

#Keep relevant columns
#meth_filter <- meth[,c(1:3,6,10:11)]
meth_filter <- as.data.frame(cbind(meth[,c(1:3,6)], t(as.data.frame(strsplit(meth$V10, " ")))[,1:2]))

rm(meth)

#Set colnames
colnames(meth_filter) <- c("chr","start","end","strand","cov","methylation_percent")

#Remove X and Y chromosomes
meth_filter <- subset(meth_filter, !(chr %in% c("chrX","chrY")))
meth_filter <- meth_filter[complete.cases(meth_filter), ]

## convert str to numeric
meth_filter$cov <- as.numeric(meth_filter$cov)
meth_filter$methylation_percent <- as.numeric(meth_filter$methylation_percent)

#Convert to decimal for comparison to array beta values
meth_filter$methylation <- meth_filter$methylation_percent/100
print("GRanges")
#Make GRanges object
cur38 <- makeGRangesFromDataFrame(meth_filter,keep.extra.columns=TRUE)
print("Load array")
#Load hg38 450K array sites 
genome(cur38) <- "hg38"
load(opt$array_file)
load(opt$probes_file)
#Identify overlaps with 450K array CpG probes
print("Identify overlaps")
meth_38_450k <- subsetByOverlaps(hm450,cur38)
index <- findOverlaps(hm450,cur38)
hm450.matched <- hm450[queryHits(index)]
mcols(hm450.matched) <- cbind.data.frame(mcols(hm450.matched),mcols(cur38[subjectHits(index)]))
methylation_sample <- hm450.matched[,c("cov","methylation")]
print("Write file")
#Write GRanges object as RDS file
fileout <- paste0(opt$out_dir,"/",opt$sample,"_methylation_hg38_HM450.RDS")

saveRDS(methylation_sample,fileout)
print("Start methylation classification")
#####################################################################
# Methylation classification#
load(opt$probes_file)
case <- methylation_sample
case$cpg <- names(case)
probes <- intersect(unique(names(case)), top_probes)
message(paste(length(probes)," overlapping CpG sites between sample and reference set.",sep=""))
write.csv(probes,file = paste0(opt$out_dir,"/",opt$sample,"_probes_for_training.csv"),row.names = FALSE, col.names = FALSE,quote = FALSE)

#Read training data
#betas_train <- fread(opt$training_data,select = c(probes,"Dx"),stringsAsFactors = FALSE)
load(opt$training_data)

probes <- intersect(unique(names(case)), colnames(betas))
probes <- c(probes,"Dx")
betas_train <- betas[,match(probes, colnames(betas))]
#betas_train <- betas[,c(probes,"Dx")]#Subset of overlapping beta values
#betas_train <- as.data.frame(betas_train)
rm(betas)
#betas_train <- t(betas_train)
#betas_train <- betas[,c(probes,"Dx") %in% colnames(betas)]

#Select most variable probes
max_CpG=10000
sds <- colSds(as.matrix(betas_train[,-ncol(betas_train)]), na.rm=F)
maxSDs <- head(order(sds,decreasing=T),n = min(ncol(betas_train)-1, max_CpG))
ts <- betas_train[,c(maxSDs,ncol(betas_train))]
#ts <- (ts >= .5) +0 #Binarise probes
rm(betas_train)
cols <- colnames(ts[,-ncol(ts)])

#Get fractions
ts$Dx <- as.factor(ts$Dx)
Dx_fractions <- min(summary(ts$Dx,maxsum=10000)) / summary(ts$Dx,maxsum=10000)

#rf <- ranger(dependent.variable.name = "Dx", data = ts[,cols], num.trees=20000, probability = T, sample.fraction = Dx_fractions,verbose=TRUE,num.threads=16)

rf <- ranger(dependent.variable.name = "Dx", data = ts[,c(cols,"Dx")], num.trees=20000, probability = T, sample.fraction = Dx_fractions,verbose=TRUE,num.threads=opt$threads)

probs <- predict(rf, ts, predict.all = F)$predictions
scores <- unlist(lapply(1:dim(probs)[1], function(i){probs[i,which(probs[i,] == max(probs[i,]))]}))
pred <- attr(scores, "name")
classes <- levels(ts$Dx)

# A GLM is trained for each class individually (1 vs all)
glm_models <- lapply(classes, function(type){
  scores <- probs[, type]
  scale_this <- data.frame(class = ifelse(pred == type & pred == ts$Dx, 1, 0), score = scores)
  scaled_scores <- glm(class ~ score, scale_this, family = binomial)
  return(scaled_scores)
})

###### For paper
details <- paste(paste0("450k_overlap: ", length(unique(names(case)))), 
                           paste0("Top_100k_overlap: ", length(probes)),
                           paste0("OOB_error: ", rf$prediction.error),
                           paste0("Training_probes: ", length(cols)),
                           sep = "\n")
write.table(details, file = paste0(opt$out_dir,"/",opt$sample,"_rf_details.tsv"), row.names = F, col.names = F, quote = F)


### predict case

case <- as.data.frame(unique(case))
case <- case[match(cols,case$cpg),]
case <- case[,c("methylation","cpg")]
case$methylation <- ((case$methylation >= .5) +0)
df <- subset(case,select=methylation)
df <- t(df)
x_probs <- predict(rf, rbind(df[,cols],df[,cols]), predict.all=F)$predictions[1,]
x_score <- x_probs[which(x_probs == max(x_probs))]
x_pred <- attr(x_score, "name")


votes <- data.frame(x_probs)
colnames(votes) <- c("Freq")

votes$Freq <- votes$Freq / sum(votes$Freq) * 100

# Apply recalibration to case score
x_scaled <- lapply(levels(ts$Dx), function(type){
  x_scores <- x_probs[type]
  scaled_scores <- glm_models[[which(classes == type)]]
  x_scaled <- predict(scaled_scores, newdata = data.frame(score = x_scores), type = "response")
  return(x_scaled)
})
x_mean_this <- unlist(x_scaled)
x_calibrated_scores <- x_mean_this/sum(x_mean_this)
x_calibrated_score <- x_calibrated_scores[x_pred]

votes$cal_Freq <- x_calibrated_scores
votes$cal_Freq <- votes$cal_Freq / sum(votes$cal_Freq) * 100

votes <- votes[order(votes$Freq),, drop = FALSE] 

### Save calibration report

report <- paste(paste0("Number of features: ", rf$num.independent.variables), 
                paste0("Predicted Class: ", x_pred),
                paste0("Initial Score: ", x_score),
                paste0("Calibrated Score: ", x_calibrated_score),
                sep = "\n")
write.table(report, file = paste0(opt$out_dir,"/",opt$sample,"_calibrated_classification.tsv"), row.names = F, col.names = F, quote = F)

write.table(votes, file = paste0(opt$out_dir,"/",opt$sample,"_votes.tsv"), quote = F)
