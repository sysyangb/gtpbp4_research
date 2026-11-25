#!/usr/bin/env Rscript
# Check the linearized patient bootstrap against exactly solved, fixed-dispersion
# weighted NB score equations for RBP genes near the FC1.5 boundary.
Sys.setenv(OPENBLAS_NUM_THREADS="1",OMP_NUM_THREADS="1",VECLIB_MAXIMUM_THREADS="1")
root<-normalizePath(getwd());out<-file.path(root,"outputs/TNBC_RBP_proportion_inference_20261027")
.libPaths(c(file.path(root,"TCGA/.Rlib"),.libPaths()))
suppressPackageStartupMessages({library(DESeq2);library(data.table)})
setDTthreads(1L)
m<-readRDS(file.path(out,"QC/patient_bootstrap_model.rds"))
fc_cut<-as.numeric(Sys.getenv("RBP_FC","1.5"));lfc_cut<-log2(fc_cut)
suffix<-if(fc_cut==1.5) "" else paste0("_FC",gsub("\\.","p",as.character(fc_cut)))
pool<-which(m$is_rbp & rowSums(m$valid)==3)
sel<-sort(unique(unlist(lapply(1:3,function(k)head(pool[order(abs(abs(m$beta[pool,k])-lfc_cut))],12)))))
d<-readRDS(file.path(root,"outputs/TNBC_DESeq2_RBP_20261026_FC1p5/QC/fresh_fitted_DESeq2.rds"))
idx<-match(m$ids[sel],rownames(d));y<-counts(d)[idx,,drop=FALSE]
mu<-assay(d,"mu")[idx,,drop=FALSE];alpha<-dispersions(d)[idx]
rm(d);invisible(gc())
glist<-lapply(c("NORMAL",m$groups),function(g)which(m$meta$condition==g))
names(glist)<-c("NORMAL",m$groups)
solve_shift<-function(yy,mm,a,w) {
  fun<-function(t)sum(w*(yy-mm*exp(t))/(1+a*mm*exp(t)))
  interval<-c(-2,2)
  if(fun(-2)*fun(2)>0)interval<-c(-10,10)
  uniroot(fun,interval,tol=1e-10)$root/log(2)
}
base<-sapply(glist,function(s)vapply(seq_along(sel),function(j)
  solve_shift(y[j,s],mu[j,s],alpha[j],rep(1,length(s))),numeric(1)))
set.seed(20261082);M<-40L;records<-list()
for(b in seq_len(M)) {
  wpat<-integer(length(m$patient))
  for(ids in m$strata)wpat[ids]<-tabulate(sample.int(length(ids),length(ids),replace=TRUE),nbins=length(ids))
  one<-exact<-matrix(NA_real_,length(sel),4)
  for(k in 1:4) {
    s<-glist[[k]];w<-wpat[m$patient_index[s]]
    one[,k]<-as.numeric(m$influence[[k]][sel,,drop=FALSE] %*% (w-1))
    exact[,k]<-vapply(seq_along(sel),function(j)solve_shift(y[j,s],mu[j,s],alpha[j],w),numeric(1))-base[,k]
  }
  for(k in 1:3) {
    a<-m$beta[sel,k]+one[,k+1]-one[,1]
    e<-m$beta[sel,k]+exact[,k+1]-exact[,1]
    records[[length(records)+1L]]<-data.table(iteration=b,group=m$groups[k],gene_id=m$ids[sel],
      one_step_log2FC=a,exact_fixed_nuisance_log2FC=e,error_log2FC=a-e,
      FC_boundary_decision_agrees=(abs(a)>=lfc_cut)==(abs(e)>=lfc_cut))
  }
}
r<-rbindlist(records)
summary<-r[,.(n=.N,median_abs_error=median(abs(error_log2FC)),
  q95_abs_error=quantile(abs(error_log2FC),.95),max_abs_error=max(abs(error_log2FC)),
  FC_boundary_agreement=mean(FC_boundary_decision_agrees)),by=group]
write.table(summary,file.path(out,"QC",paste0("exact_score_calibration_summary",suffix,".tsv")),sep="\t",row.names=FALSE,quote=TRUE)
con<-gzfile(file.path(out,"QC",paste0("exact_score_calibration_details",suffix,".tsv.gz")),"wt")
write.table(r,con,sep="\t",row.names=FALSE,quote=TRUE);close(con)
print(summary)
cat("Calibration:",length(sel),"boundary-adjacent RBP genes x",M,"joint patient draws\n")
