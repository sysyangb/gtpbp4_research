#!/usr/bin/env Rscript
args<-commandArgs(trailingOnly=TRUE)
root<-normalizePath(if(length(args))args[1] else getwd())
Sys.setenv(RBP_FIGURE_FC="2",RBP_FIGURE_STARS="1",RBP_FIGURE_PREFIX="Figure1A",
           RBP_FIGURE_DIR=file.path(root,"outputs/Figure1A_FC2"))
dest<-file.path(root,"outputs/Figure1A_FC2")
builder<-file.path(root,"outputs/TNBC_DESeq2_RBP_20261026_FC1p5/scripts/plot_rbp_percent_stacked_Holm_14pt.R")
local_builder<-file.path(dest,"Figure1A_plot_builder.R")
if(file.exists(builder)) file.copy(builder,local_builder,overwrite=TRUE)
source(local_builder)
