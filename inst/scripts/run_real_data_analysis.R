#!/usr/bin/env Rscript

suppressPackageStartupMessages(library(itsne))

args <- commandArgs(trailingOnly = TRUE)
output_dir <- if (length(args)) args[[1]] else "results/real-data"

run_real_data_analysis(output_dir = output_dir)
