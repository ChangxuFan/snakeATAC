#!/opt/apps/R/3.6.1/bin/Rscript
suppressMessages(library(dplyr))
suppressMessages(library(liteRnaSeqFanc))
suppressMessages(library("argparse"))
options(stringsAsFactors = F)

parser <- argparse::ArgumentParser(add_help = T)

parser$add_argument("-d", "--snake_dir", required = F, default="./")
parser$add_argument("-b", "--bed",required = T)
parser$add_argument("-s", "--subset",required = F, default = 42.001, type = "integer")
parser$add_argument("-t", "--threads",required = F, default = 12, type = "integer")
parser$add_argument("-S", "--scale_fac",required = F, default = 3, type = "integer")
parser$add_argument("-e", "--ext",required = F, default = 1000, type = "integer")
parser$add_argument("-k", "--skip_subset",required = F, action = "store_true")

args <- parser$parse_args()

