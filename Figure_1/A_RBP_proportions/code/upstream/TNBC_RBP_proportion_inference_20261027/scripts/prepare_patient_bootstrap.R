#!/usr/bin/env Rscript
# Joint patient bootstrap preparation for the FC1.5 RBP composition statistic.
# This is a one-step, fixed-nuisance approximation to patient-level resampling.
Sys.setenv(OPENBLAS_NUM_THREADS="2",OMP_NUM_THREADS="2",VECLIB_MAXIMUM_THREADS="2")
root <- normalizePath(getwd())
old <- file.path(root,"outputs/TNBC_DESeq2_RBP_20261026_FC1p5")
out <- file.path(root,"outputs/TNBC_RBP_proportion_inference_20261027")
.libPaths(c(file.path(root,"TCGA/.Rlib"),.libPaths()))
suppressPackageStartupMessages({library(DESeq2);library(data.table)})
setDTthreads(1L)
write_tsv <- function(x,name) write.table(as.data.frame(x),file.path(out,"QC",name),
  sep="\t",row.names=FALSE,quote=TRUE,na="NA",qmethod="double")
groups <- c("TN","HER2","HR")
raw <- lapply(groups,function(g) as.data.table(read.delim(
  gzfile(file.path(old,"results",paste0(g,"_vs_NORMAL_all_tested_genes.tsv.gz"))),check.names=FALSE)))
all_ids <- raw[[1]]$gene_id
stopifnot(all(vapply(raw,function(d)setequal(d$gene_id,all_ids),logical(1))))
raw <- lapply(raw,function(d)d[match(all_ids,gene_id)])
eligible <- Reduce(`|`,lapply(raw,function(d)!is.na(d$padj)))
ids <- all_ids[eligible]
beta <- sapply(raw,function(d)d$log2FoldChange[eligible]);colnames(beta)<-groups
se <- sapply(raw,function(d)d$lfcSE[eligible]);colnames(se)<-groups
valid <- sapply(raw,function(d)!is.na(d$padj[eligible]));colnames(valid)<-groups
pval <- sapply(raw,function(d)d$pvalue[eligible]);colnames(pval)<-groups
original_padj <- sapply(raw,function(d)d$padj[eligible]);colnames(original_padj)<-groups
rbp <- fread(file.path(old,"metadata/RBP_catalog_and_coverage.tsv"))
is_rbp <- ids %in% rbp$gene_id
obs <- fread(file.path(old,"results/RBP_up_down_summary.tsv"))[match(paste0(groups,"_vs_NORMAL"),comparison)]
observed <- as.matrix(obs[,.(RBP_up,RBP_down,RBP_total)])
rownames(observed)<-groups

message("Loading original four-group fitted model: ",format(Sys.time()))
dds <- readRDS(file.path(old,"QC/fresh_fitted_DESeq2.rds"))
idx <- match(ids,rownames(dds));stopifnot(!anyNA(idx))
meta <- fread(file.path(old,"metadata/selected_samples.tsv"))[match(colnames(dds),sample_id)]
stopifnot(identical(meta$sample_id,colnames(dds)),all(mcols(dds)$betaConv))
patient <- unique(meta$patient_id)
patient_index <- match(meta$patient_id,patient)
patient_pattern <- vapply(patient,function(p)paste(sort(meta[patient_id==p,condition]),collapse="+"),character(1))
strata <- split(seq_along(patient),patient_pattern)
stopifnot(all(lengths(strata)>1L))
write_tsv(data.table(patient_id=patient,stratum=patient_pattern),"patient_strata.tsv")
write_tsv(meta[,.(sample_id,patient_id,condition)],"sample_cluster_mapping.tsv")
disp <- dispersions(dds)[idx]
grp_order <- c("NORMAL",groups)
influence <- list(); information <- matrix(NA_real_,length(ids),4,dimnames=list(ids,grp_order))
fd_checks <- list()
for(g in grp_order) {
  s <- which(meta$condition==g)
  mu <- assay(dds,"mu")[idx,s,drop=FALSE]
  y <- counts(dds)[idx,s,drop=FALSE]
  # NB log-mean score and expected Fisher information. The 0.5 lower bound
  # matches the minmu convention of the fitted bulk DESeq2 Wald model.
  H <- rowSums(pmax(mu,0.5)/(1+disp*pmax(mu,0.5)))
  score <- (y-mu)/(1+disp*mu)
  # Observed score Jacobian gives the derivative of the weighted likelihood
  # estimator. Expected information is retained separately for the DESeq2 SE check.
  Hobs_all <- rowSums(mu*(1+disp*y)/(1+disp*mu)^2)
  stopifnot(all(is.finite(Hobs_all)),all(Hobs_all>0))
  influence[[g]] <- score/(Hobs_all*log(2))
  information[,g] <- H
  # Finite-difference check on well-expressed genes, independent of model SEs.
  candidates <- which(rowMeans(mu)>100 & rowMeans(y)>100)
  set.seed(20261027+match(g,grp_order))
  selected <- sample(candidates,min(12L,length(candidates)))
  for(j in selected) {
    yy<-as.numeric(y[j,]);mm<-as.numeric(mu[j,]);a<-disp[j]
    # Solve the actual weighted NB score before and after a tiny weight change.
    scorefun <- function(t,w)sum(w*(yy-mm*exp(t))/(1+a*mm*exp(t)))
    base <- uniroot(function(t)scorefun(t,rep(1,length(s))),c(-1,1),tol=1e-12)$root
    mm0 <- mm*exp(base)
    k <- which.max(abs((yy-mm0)/(1+a*mm0)))
    eps <- 1e-4;w<-rep(1,length(s));w[k]<-1+eps
    changed <- uniroot(function(t)scorefun(t,w),c(-1,1),tol=1e-12)$root
    numerical <- (changed-base)/eps/log(2)
    # Check the same observed-Jacobian derivative used by the bootstrap.
    Hobs <- sum(mm0*(1+a*yy)/(1+a*mm0)^2)
    derivative <- (yy[k]-mm0[k])/(1+a*mm0[k])/Hobs/log(2)
    rel <- abs(numerical-derivative)/max(abs(derivative),1e-12)
    stopifnot(rel<0.005)
    fd_checks[[length(fd_checks)+1L]] <- data.table(group=g,gene_id=ids[j],
      finite_difference=numerical,score_derivative=derivative,relative_error=rel)
  }
  rm(y,mu,score,Hobs_all);invisible(gc())
  message("Prepared ",g,": ",length(s)," samples; ",length(ids)," genes")
}
write_tsv(rbindlist(fd_checks),"score_finite_difference_validation.tsv")
# Check the information normalization against DESeq2's exported coefficient SE.
info_se <- sapply(groups,function(g)sqrt(1/information[,"NORMAL"]+1/information[,g])/log(2))
se_ratio <- info_se/se
write_tsv(data.table(group=groups,median_ratio=apply(se_ratio,2,median,na.rm=TRUE),
  q01_ratio=apply(se_ratio,2,quantile,probs=.01,na.rm=TRUE),
  q99_ratio=apply(se_ratio,2,quantile,probs=.99,na.rm=TRUE)),"fisher_information_SE_check.tsv")
stopifnot(all(abs(apply(se_ratio,2,median,na.rm=TRUE)-1)<0.01))
# Baseline FDR and FC decisions must reproduce every bar in the original figure.
for(k in seq_along(groups)) {
  q <- rep(NA_real_,nrow(beta));q[valid[,k]]<-p.adjust(pval[valid[,k],k],"BH")
  stopifnot(max(abs(q-original_padj[,k]),na.rm=TRUE)<1e-12)
  up <- sum(is_rbp & !is.na(q) & q<1e-10 & beta[,k]>=log2(1.5))
  down <- sum(is_rbp & !is.na(q) & q<1e-10 & beta[,k]<=-log2(1.5))
  stopifnot(up==observed[k,"RBP_up"],down==observed[k,"RBP_down"])
}
rm(dds);invisible(gc())
model <- list(ids=ids,groups=groups,beta=beta,se=se,valid=valid,is_rbp=is_rbp,
  influence=influence,information=information,meta=as.data.frame(meta),patient=patient,
  patient_index=patient_index,strata=strata,observed=observed,alpha=1e-10,fc=1.5,
  n_original_genes=nrow(raw[[1]]),seed=20261027L)
saveRDS(model,file.path(out,"QC/patient_bootstrap_model.rds"),compress=FALSE)
write_tsv(data.table(metric=c("samples","patients","paired_patients","genes_in_BH_family_union","RBP_in_family"),
  value=c(nrow(meta),length(patient),sum(table(meta$patient_id)>1),length(ids),sum(is_rbp))),"model_summary.tsv")
writeLines(capture.output(sessionInfo()),file.path(out,"QC/R_sessionInfo.txt"))
message("Prepared and validated bootstrap model: ",format(Sys.time()))
