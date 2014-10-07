# Create data for cytobands fro UCSC data

temp <- tempfile(fileext = ".txt.gz")
for (genome in c("hg18","hg19","hg38","mm10")) {
  download.file(paste0("http://hgdownload.cse.ucsc.edu/goldenPath/",genome,"/database/cytoBandIdeo.txt.gz"),temp)
  data <- read.table(temp,sep="\t")
  assign(paste0("chrom.bands.",genome),data)
  save(list=paste0("chrom.bands.",genome), file = paste0('data/chrom.bands.',genome,'.rda'), compress='xz')
}
