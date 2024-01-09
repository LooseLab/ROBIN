
for (package in c('optparse', 'rmarkdown','kableExtra','knitr')) {
  if (!require(package, character.only=T, quietly=T)) {
    install.packages(package,repos = "http://cran.us.r-project.org")
    library(package, character.only=T)
  }
}

#Parse arguments
option_list = list(
  make_option(c("-p", "--prefix"), type="character", default=NULL, 
              help="prefix", metavar="character"),
  make_option(c("-m", "--mutations"), type="character", default=NULL, 
             help="mutations table", metavar="character"),
  make_option(c("-c", "--cnv_plot"), type="character", default=NULL,
              help="cnv plot", metavar="character"),
  make_option(c("-r", "--rf_details"), type="character", default=NULL,
              help="RF details tsv", metavar="character"),
  make_option(c("-v", "--votes"), type="character", default=NULL,
              help="votes file", metavar="character"),
  make_option(c("-o", "--output_dir"), type="character", default=NULL,
              help="output directory", metavar="character"),
  make_option(c("-n", "--patient"), type="character", default=NULL, 
              help="patient", metavar="character"),
  make_option(c("-e", "--coverage"), type="character", default=NULL,
              help="coverage summary", metavar="character"),
  make_option(c("-s", "--sample"), type="character", default=NULL,
              help="sample", metavar="character"),
  make_option(c("-t", "--mgmt"), type="character", default="false", 
             help="mgmt prediction", metavar="character"),
  ## following options added by GF
  make_option(c("-a", "--methylartist"), type="character", default="false",
              help="methylartist mgmt plot", metavar="character"),
  make_option(c("-u", "--report_UKHD"), type="character", default=NULL,
              help="report_UKHD R markdown doc", metavar="character"),
  make_option(c("-b", "--promoter_mgmt_coverage"), type="integer", default=NULL,
              help="average coverage at mgmt promoter", metavar="character"),
  make_option(c("-k", "--sturgeon"), type="character", default="false",
              help="check the sturgeon flag", metavar="character"),
  make_option(c("-d", "--sturgeon_pdf"), type="character", default=NULL,
              help="sturgeon output pdf", metavar="character"),
  make_option(c("-f", "--sturgeon_csv"), type="character", default=NULL,
              help="sturgeon output csv", metavar="character"),
  make_option(c("-g", "--igv_report"), type="character", default="false",
              help="IGV-report html output", metavar="character"),
  make_option(c("-i", "--seq"), type="character", default="Unknown",
             help="Platform used to sequencing; F=MinION/GridION, P=PromethION", metavar="character"),
  make_option(c("-j", "--nextflow_ver"), type="character", default=NULL,
              help="Include the version of the Nextflow pipeline used to generate the report", metavar="character")
)

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);
prefix <- opt$prefix
mutations <- opt$mutations
rf_details <- opt$rf_details
cnv_plot <- opt$cnv_plot
votes <- opt$votes
coverage <- opt$coverage
patient <- opt$patient
seq <- opt$seq
sample <- opt$sample
#mgmt <- opt$mgmt
# add option for report_UKHD above
report_UKHD <- opt$report_UKHD
#methylartist_plot <- opt$methylartist
cov <- opt$promoter_mgmt_coverage
igv_report <- opt$igv_report
#report_full <- opt$report_full
nextflow_ver <- opt$nextflow_ver

mgmt = "false"
if (file.exists(opt$mgmt)) {
    mgmt="true"
}

methylartist_plot = "false"
if (file.exists(opt$methylartist)) {
    methylartist_plot = "true"
}


## check for existience of sturgeon files (it might not have run)
#if (file.exists(opt$sturgeon_pdf)){
if (opt$sturgeon == "true"){
	sturgeon_pdf <- opt$sturgeon_pdf
	sturgeon_csv <- opt$sturgeon_csv
} else {
	sturgeon_pdf <- "NULL"
	sturgeon_csv <- "NULL"
}

# generate the report

inc_igvreport = FALSE
exc_igvreport = TRUE
# lite version
mgmt = "true"
render(report_UKHD, 
       output_format = "html_document", 
       output_dir = opt$output_dir,
       output_file = paste0(prefix,"_Rapid-CNS2_report_lite.html"))

inc_igvreport = TRUE
exc_igvreport = FALSE
# full version
render(report_UKHD,
       output_format = "html_document",
       output_dir = opt$output_dir,
       output_file = paste0(prefix,"_Rapid-CNS2_report_full.html"))
