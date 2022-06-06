#!/usr/bin/env Rscript

library(argparse)
library(ggplot2)
library(stringr)

# Get options
parser <- ArgumentParser()

parser$add_argument("-i", "--input_dir",
                    help="Specify input directory")
parser$add_argument("-o", "--output_dir",
                    help="Specify output directory for coverage plots")
parser$add_argument("-c", "--cov_cutoff",type="double",
                    help="Specify coverage cutoff for plotting samples")

args <- parser$parse_args()

input.dir <- args$input_dir
output.dir <- args$output_dir
cutoff <- args$cov_cutoff

files <- list.files(input.dir,pattern="*coverage.csv")

percent_coverage <- as.data.frame(matrix(nrow=length(files),
                                         ncol=2))
colnames(percent_coverage) <- c("sample","coverage")

percent_coverage[,1] <- gsub("[.].*","",files)
for(i in 1:length(files)){
  file <-  read.delim(paste0(input.dir,"/",files[i]),
                      sep=",")
  zero <- nrow(file[file$totalcount < cutoff,])
  total <- nrow(file)
  percent <- 100 - (zero/total * 100)
  percent_coverage[i,2] <- percent
  
  if(percent >= 50){
    plot <- ggplot(file, aes(x = ntpos, y = log10(totalcount), group = name)) +
      geom_line() +
      geom_hline(yintercept = log10(cutoff), linetype = 2) +
      theme_bw() +
      ylab(str_glue("Log10 Raw Coverage\n(dashed line: log10 coverage of {cutoff}X)")) +
      xlab("Nucleotide Position") +
      facet_grid(name ~ ., scales = "free")
    
    ggsave(plot, filename = str_glue("{percent_coverage[i,1]}.covplot.pdf"), path = output.dir,
           width = 8, height = 3, limitsize = FALSE, useDingbats = FALSE)
  }
}

write.table(percent_coverage,
            str_glue("{output.dir}/percentcov_{cutoff}x.tsv"),
            row.names = F,
            col.names = T,
            quote = F,
            sep= "\t")

