#!/usr/bin/env Rscript

suppressPackageStartupMessages(library(itsne))

args <- commandArgs(trailingOnly = TRUE)
output_dir <- if (length(args)) args[[1]] else "results/simulations"

run_simulation_studies(output_dir = output_dir)
