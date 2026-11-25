#!/usr/bin/env Rscript
# Patient-stratified, joint one-step NB-score bootstrap for DE-RBP composition.
# All patient profiles move together across genes; matched normal/tumor rows
# and all three contrasts receive the same patient weights in every draw.
# Conditional approximation: the fitted dispersion, size factors, Wald SEs,
# prefilter and Cook/independent-filter eligibility are fixed at observed values.
Sys.setenv(OPENBLAS_NUM_THREADS="2",OMP_NUM_THREADS="2",VECLIB_MAXIMUM_THREADS="2")
root<-normalizePath(getwd());out<-file.path(root,"outputs/TNBC_RBP_proportion_inference_20261027")
.libPaths(c(file.path(root,"TCGA/.Rlib"),.libPaths()))
suppressPackageStartupMessages({library(Rcpp);library(data.table)})
setDTthreads(1L)
Rcpp::sourceCpp(file.path(out,"scripts/bootstrap_decisions.cpp"),rebuild=FALSE)
model<-readRDS(file.path(out,"QC/patient_bootstrap_model.rds"))
requested_fc<-as.numeric(Sys.getenv("RBP_FC","1.5"))
stopifnot(is.finite(requested_fc),requested_fc>1)
tag<-if(requested_fc==1.5) "FC1p5" else paste0("FC",gsub("\\.","p",as.character(requested_fc)))
suffix<-if(requested_fc==1.5) "" else paste0("_",tag)
if(requested_fc!=model$fc) {
  source_dir<-file.path(root,"outputs/TNBC_DESeq2_RBP_20261026_FC1p5/results")
  selected<-list();summary_rows<-list()
  for(g in model$groups) {
    tab<-fread(file.path(source_dir,paste0(g,"_vs_NORMAL_all_tested_RBP.tsv")))
    tab[,direction:=fcase(!is.na(padj)&padj<model$alpha&log2FoldChange>=log2(requested_fc),"Up",
                          !is.na(padj)&padj<model$alpha&log2FoldChange<=-log2(requested_fc),"Down",
                          default="Not_significant")]
    sig<-tab[direction!="Not_significant"];selected[[g]]<-sig
    write.table(sig,file.path(out,"results",paste0(tag,"_",g,"_vs_NORMAL_DE_RBP.tsv")),
                sep="\t",row.names=FALSE,quote=TRUE,na="NA")
    summary_rows[[g]]<-data.table(comparison=paste0(g,"_vs_NORMAL"),
      RBP_up=sum(sig$direction=="Up"),RBP_down=sum(sig$direction=="Down"),RBP_total=nrow(sig))
  }
  observed_table<-rbindlist(summary_rows)
  write.table(observed_table,file.path(out,"results",paste0(tag,"_RBP_up_down_summary.tsv")),
              sep="\t",row.names=FALSE,quote=TRUE)
  write.table(rbindlist(selected),file.path(out,"results",paste0(tag,"_three_contrasts_significant_RBP.tsv")),
              sep="\t",row.names=FALSE,quote=TRUE,na="NA")
  model$observed<-as.matrix(observed_table[,.(RBP_up,RBP_down,RBP_total)])
  rownames(model$observed)<-model$groups;model$fc<-requested_fc
}
B<-as.integer(Sys.getenv("RBP_BOOTSTRAP_REPLICATES","19999"))
stopifnot(B>=999L)
batch_size<-64L;set.seed(model$seed)
groups<-model$groups;group_indices<-lapply(c("NORMAL",groups),function(g)which(model$meta$condition==g))
names(group_indices)<-c("NORMAL",groups)
zero<-replicate(4,matrix(0,nrow(model$beta),1),simplify=FALSE)
baseline<-bootstrap_counts(zero,model$beta,model$se,model$valid,model$is_rbp,model$alpha,log2(model$fc))
for(k in seq_along(groups)) stopifnot(baseline[1,3*k-2]==model$observed[k,"RBP_up"],
                                     baseline[1,3*k-1]==model$observed[k,"RBP_down"])
cat("Baseline decisions match original graph exactly\n")
npatient<-length(model$patient);draws<-matrix(NA_real_,B,9)
colnames(draws)<-as.vector(rbind(paste0(groups,"_up"),paste0(groups,"_down"),paste0(groups,"_all_DE")))
st<-proc.time()[[3]];last_progress<-0L
for(start in seq.int(1L,B,by=batch_size)) {
  end<-min(B,start+batch_size-1L);nb<-end-start+1L
  weights<-matrix(0,npatient,nb)
  for(ids in model$strata) {
    n<-length(ids)
    for(b in seq_len(nb)) weights[ids,b]<-tabulate(sample.int(n,n,replace=TRUE),nbins=n)
  }
  # Exact patient resampling counts in the seven observed sample-pattern strata.
  stopifnot(all(colSums(weights)==npatient))
  delta<-lapply(seq_along(group_indices),function(k) {
    s<-group_indices[[k]]
    w<-weights[model$patient_index[s],,drop=FALSE]
    stopifnot(all(colSums(w)==length(s)))
    model$influence[[names(group_indices)[k]]] %*% (w-1)
  })
  draws[start:end,]<-bootstrap_counts(delta,model$beta,model$se,model$valid,model$is_rbp,
                                    model$alpha,log2(model$fc))
  if(end-last_progress>=1024L || end==B) {
    cat("Completed ",end,"/",B,"; elapsed ",round(proc.time()[[3]]-st,1)," sec\n",sep="")
    flush.console();last_progress<-end
    saveRDS(list(completed=end,draws=draws[seq_len(end),,drop=FALSE],seed=.Random.seed),
            file.path(out,"QC",paste0("bootstrap_progress",suffix,".rds")),compress=FALSE)
  }
}
stopifnot(!anyNA(draws))
props<-sapply(seq_along(groups),function(k)draws[,3*k-2]/(draws[,3*k-2]+draws[,3*k-1]))
colnames(props)<-groups
stopifnot(all(is.finite(props)))
observed_prop<-model$observed[,"RBP_up"]/rowSums(model$observed[,c("RBP_up","RBP_down"),drop=FALSE])
pairs<-combn(seq_along(groups),2L)
results<-lapply(seq_len(ncol(pairs)),function(j) {
  a<-pairs[1,j];b<-pairs[2,j];d<-observed_prop[a]-observed_prop[b]
  boot_d<-props[,a]-props[,b];centered<-boot_d-d
  exceed<-sum(abs(centered)>=abs(d));p<-(exceed+1)/(B+1)
  mc<-binom.test(exceed,B)$conf.int
  basic_ci<-2*d-rev(quantile(boot_d,c(.025,.975),names=FALSE))
  data.table(comparison=paste(groups[c(a,b)],collapse=" vs "),group1=groups[a],group2=groups[b],
    up_fraction_group1=observed_prop[a],up_fraction_group2=observed_prop[b],
    difference_percentage_points=100*d,bootstrap_SE_percentage_points=100*sd(boot_d),
    bootstrap_bias_percentage_points=100*(mean(boot_d)-d),
    CI95_basic_low_pp=100*basic_ci[1],CI95_basic_high_pp=100*basic_ci[2],
    p_raw=p,tail_exceedances=exceed,B=B,MC_p_low=mc[1],MC_p_high=mc[2],
    normal_approx_p=2*pnorm(-abs(d)/sd(boot_d)))
})
ans<-rbindlist(results);ans[,`:=`(p_Holm=p.adjust(p_raw,"holm"),p_BH=p.adjust(p_raw,"BH"))]
ans[,significant_Holm_005:=p_Holm<.05]
write_tsv<-function(d,p){con<-if(grepl("\\.gz$",p))gzfile(p,"wt")else file(p,"wt");
  on.exit(close(con));write.table(as.data.frame(d),con,sep="\t",row.names=FALSE,quote=TRUE,na="NA",qmethod="double")}
write_tsv(ans,file.path(out,"results",paste0(tag,"_pairwise_patient_bootstrap.tsv")))
draw_table<-as.data.table(draws);draw_table[,iteration:=.I]
draw_table[,c("TN_up_fraction","HER2_up_fraction","HR_up_fraction"):=as.data.table(props)]
write_tsv(draw_table,file.path(out,"results",paste0(tag,"_joint_patient_bootstrap_draws.tsv.gz")))
write_tsv(data.table(group=groups,observed_up_fraction=observed_prop,
  bootstrap_mean=colMeans(props),bootstrap_sd=apply(props,2,sd),
  CI95_percentile_low=apply(props,2,quantile,probs=.025),
  CI95_percentile_high=apply(props,2,quantile,probs=.975)),
  file.path(out,"results",paste0(tag,"_proportion_uncertainty.tsv")))
write_tsv(as.data.table(cor(props),keep.rownames="group"),
          file.path(out,"QC",paste0("joint_proportion_correlation",suffix,".tsv")))
saveRDS(list(results=ans,draws=draws,proportions=props,observed=observed_prop,
             B=B,seed=model$seed,fc=model$fc,observed_counts=model$observed,
             elapsed_seconds=proc.time()[[3]]-st),
        file.path(out,"QC",paste0("completed_bootstrap",suffix,".rds")),compress=FALSE)
writeLines(capture.output(sessionInfo()),file.path(out,"QC",paste0("bootstrap_R_sessionInfo",suffix,".txt")))
print(ans,digits=8)
cat("FINISHED: ",format(Sys.time()),"; ",B," joint patient draws\n",sep="")
