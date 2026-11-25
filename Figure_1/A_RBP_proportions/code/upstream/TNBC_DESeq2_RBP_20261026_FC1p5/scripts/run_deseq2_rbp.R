#!/usr/bin/env Rscript
# Refit raw TCGA-BRCA counts, then intersect each tumor-vs-normal result with RBP.
Sys.setenv(OMP_NUM_THREADS="1", OPENBLAS_NUM_THREADS="1",
           MKL_NUM_THREADS="1", VECLIB_MAXIMUM_THREADS="1", TZ="Asia/Shanghai")
args <- commandArgs(trailingOnly=TRUE)
root <- normalizePath(if(length(args)) args[1] else getwd())
out <- file.path(root, "outputs", "TNBC_DESeq2_RBP_20261026_FC1p5")
input <- file.path(root, "TCGA")
.libPaths(c(file.path(input, ".Rlib"), .libPaths()))
suppressPackageStartupMessages({library(data.table); library(DESeq2); library(ggplot2)})
stopifnot(requireNamespace("ragg", quietly=TRUE), requireNamespace("svglite", quietly=TRUE),
          requireNamespace("scales", quietly=TRUE))
setDTthreads(1L); options(mc.cores=1L, width=140); set.seed(20261026)
for(d in c("results", "metadata", "QC", "figures")) dir.create(file.path(out,d), recursive=TRUE, showWarnings=FALSE)
# Base R preserves subnormal probabilities (< 2.23e-308); the installed fwrite
# version mis-serializes these values. Do not clamp or replace very small P values.
write_tsv <- function(x, path) {
  target <- file.path(out,path)
  con <- if(grepl("\\.gz$",target)) gzfile(target,"wt") else file(target,"wt")
  on.exit(close(con))
  write.table(as.data.frame(x),con,sep="\t",row.names=FALSE,quote=TRUE,na="NA",qmethod="double")
}
strip_version <- function(x) sub("\\.[0-9]+$", "", x)
alpha <- 1e-10; fc <- 1.5; lfc <- log2(fc)
files <- c(counts="TCGA-BRCA_counts.txt", phenotype="TCGA-BRCA_HR_HER2_NORMAL.condition_v2.tsv",
           RBP="Homo_sapiens.RBP.txt", clinical="TCGA-BRCA_clinical.txt",
           clinical_biotab="external_receptor_sources/nationwidechildrens.org_clinical_patient_brca.txt",
           annotation="DESeq2_TN_RBP_analysis/KEGG_enrichment_top500_TN_all_genes/reference_snapshot/Ensembl_BioMart_hsapiens_gene_annotation.tsv")
paths <- file.path(input,files); names(paths) <- names(files)
stopifnot(all(file.exists(paths)))
write_tsv(data.table(input=names(files),relative_path=unname(files),
                    bytes=file.info(paths)$size, md5=unname(tools::md5sum(paths))), "QC/input_manifest.tsv")
export_only <- "--export-only" %in% args
message("Started ", format(Sys.time()), "; DESeq2 ", packageVersion("DESeq2"),
        if(export_only) "; export from this run's fresh fit" else "; new fit, no old checkpoint")

# Preserve the established receptor-based cohort; do not infer subtypes from RNA.
meta <- fread(paths["phenotype"], colClasses=c(sample_type_code="character"), na.strings=NULL)
stopifnot(!anyDuplicated(meta$sample_id))
header <- strsplit(readLines(paths["clinical_biotab"],n=1),"\t",fixed=TRUE)[[1]]
clinical <- fread(paths["clinical_biotab"], skip=3, header=FALSE, col.names=header)
sex <- unique(clinical[,.(patient_id=bcr_patient_barcode,sex=toupper(gender))])
simple <- suppressWarnings(fread(paths["clinical"]))
sex <- rbind(sex, unique(simple[!submitter_id %in% sex$patient_id,.(patient_id=submitter_id,sex=toupper(gender))]))
stopifnot(!anyDuplicated(sex$patient_id))
meta <- merge(meta,sex,by="patient_id",all.x=TRUE,sort=FALSE)
eligible <- meta[sex=="FEMALE" & ((condition=="NORMAL" & sample_type=="SOLID_TISSUE_NORMAL") |
                                  (condition %in% c("TN","HR","HER2") & sample_type=="PRIMARY_TUMOR"))]
dt <- fread(paths["counts"],select=c("V1",eligible$sample_id),check.names=FALSE,showProgress=FALSE)
setnames(dt,"V1","gene_id_version")
stopifnot(all(eligible$sample_id %in% names(dt)))
eligible[,library_size:=vapply(sample_id,function(id) sum(dt[[id]]),numeric(1))]
setorder(eligible,condition,patient_id,-library_size,sample_id)
selected <- eligible[, .SD[1L], by=.(condition,patient_id)]
selected[,order_internal:=match(condition,c("NORMAL","TN","HR","HER2"))]
setorder(selected,order_internal,patient_id); selected[,order_internal:=NULL]
excluded <- meta[!sample_id %in% selected$sample_id]
excluded[,reason:=fcase(sex=="MALE","male",is.na(sex),"missing_sex",
  sample_type=="METASTATIC","metastatic",condition=="UNKNOWN","unknown_subtype",
  sample_id %in% eligible$sample_id,"duplicate_aliquot",default="ineligible_sample_type")]
write_tsv(selected,"metadata/selected_samples.tsv"); write_tsv(excluded,"metadata/excluded_samples.tsv")
group_n <- selected[,.(N=.N),by=condition][match(c("TN","HR","HER2","NORMAL"),condition)]
stopifnot(all(group_n$N >= 3), !anyDuplicated(selected$sample_id),
          !anyDuplicated(selected[,.(condition,patient_id)]))
write_tsv(group_n,"metadata/sample_group_counts.tsv"); print(group_n)
old <- fread(file.path(input,"DESeq2_four_group_marker_analysis/metadata/selected_samples.tsv"))
stopifnot(setequal(selected$sample_id,old$sample_id))

ids <- selected$sample_id
mat <- as.matrix(dt[,..ids]); gene_id <- strip_version(dt$gene_id_version)
stopifnot(is.numeric(mat), !anyNA(mat), all(is.finite(mat)), all(mat>=0),
          all(mat==floor(mat)), max(mat)<=.Machine$integer.max)
storage.mode(mat) <- "integer"
raw_n <- nrow(mat); duplicate_n <- sum(duplicated(gene_id))
if(duplicate_n) mat <- rowsum(mat,group=gene_id,reorder=FALSE) else rownames(mat) <- gene_id
stopifnot(max(mat)<=.Machine$integer.max); storage.mode(mat)<-"integer"
rm(dt); invisible(gc())
keep <- rowSums(mat>=10L)>=10L
rbp <- fread(paths["RBP"],header=FALSE,sep="\t",colClasses="character",
  col.names=c("gene_id","gene_symbol","entrez_id","gene_description","aliases"))
rbp[,gene_id:=strip_version(gene_id)]
stopifnot(nrow(rbp)==2961L,!anyDuplicated(rbp$gene_id))
rbp[,status:=fcase(!gene_id %in% rownames(mat),"absent_from_counts",
                  gene_id %in% rownames(mat)[keep],"tested",default="low_count_filtered")]
write_tsv(rbp,"metadata/RBP_catalog_and_coverage.tsv")
write_tsv(data.table(metric=c("input_gene_rows","duplicate_unversioned_ids","unique_gene_ids",
  "genes_after_prefilter","genes_removed_low_counts","RBP_catalog","RBP_in_counts","RBP_after_prefilter"),
  value=c(raw_n,duplicate_n,nrow(mat),sum(keep),sum(!keep),nrow(rbp),
  sum(rbp$status!="absent_from_counts"),sum(rbp$status=="tested"))),"QC/gene_filter_summary.tsv")
coldata <- as.data.frame(selected); rownames(coldata)<-coldata$sample_id
coldata$condition <- factor(coldata$condition,levels=c("NORMAL","TN","HR","HER2"))
stopifnot(identical(colnames(mat),rownames(coldata)))
dds <- DESeqDataSetFromMatrix(mat[keep,,drop=FALSE],coldata,design=~condition)
rm(mat); invisible(gc())
message("Fitting ", nrow(dds), " genes across ", ncol(dds), " samples")
# Match the existing project analysis: no automatic replacement of outlier counts.
# Cook's-distance filtering remains active in results().
if(export_only) {
  fresh <- readRDS(file.path(out,"QC/fresh_fitted_DESeq2.rds"))
  stopifnot(identical(rownames(fresh),rownames(dds)),identical(colnames(fresh),colnames(dds)),
            identical(as.character(fresh$condition),as.character(dds$condition)),
            identical(counts(fresh),counts(dds)),identical(design(fresh),design(dds)))
  dds <- fresh; rm(fresh); invisible(gc())
  message("Validated original counts and metadata against this run's fresh fit")
} else {
  dds <- DESeq(dds,fitType="parametric",betaPrior=FALSE,parallel=FALSE,
               minReplicatesForReplace=Inf,quiet=FALSE)
  message("Fresh DESeq2 fit finished: ", format(Sys.time()))
  saveRDS(dds,file.path(out,"QC/fresh_fitted_DESeq2.rds"),compress=FALSE)
}
write_tsv(data.table(sample_id=colnames(dds),condition=as.character(dds$condition),
                    size_factor=sizeFactors(dds)),"QC/size_factors.tsv")
write_tsv(data.table(gene_id=rownames(dds),dispersion=dispersions(dds),
                    beta_converged=mcols(dds)$betaConv),"QC/model_diagnostics.tsv")

annot <- fread(paths["annotation"],colClasses="character")
annot[is.na(hgnc_symbol),hgnc_symbol:=""]
annot[,nonempty:=as.integer(nzchar(hgnc_symbol))]; setorder(annot,ensembl_gene_id,-nonempty)
annot <- unique(annot,by="ensembl_gene_id"); annot[,nonempty:=NULL]
summary_list<-list(); rbp_all<-list(); rbp_sig<-list(); filtering<-list()
summary_noFC <- list(); rbp_sig_noFC <- list()
for(g in c("TN","HR","HER2")) {
  contrast_name <- paste0(g,"_vs_NORMAL")
  message("Extracting ",contrast_name,"; padj < ",alpha,"; |log2FC| >= ",lfc)
  res <- results(dds,contrast=c("condition",g,"NORMAL"),alpha=alpha,
                 pAdjustMethod="BH",independentFiltering=TRUE,cooksCutoff=TRUE,lfcThreshold=0)
  tab <- as.data.table(as.data.frame(res),keep.rownames="gene_id")
  tab <- merge(tab,annot,by.x="gene_id",by.y="ensembl_gene_id",all.x=TRUE,sort=FALSE)
  idx <- match(tab$gene_id,rbp$gene_id)
  tab[,`:=`(contrast=contrast_name,is_RBP=!is.na(idx),
      RBP_symbol=rbp$gene_symbol[idx],RBP_description=rbp$gene_description[idx])]
  tab[,gene_symbol:=fifelse(!is.na(RBP_symbol),RBP_symbol,hgnc_symbol)]
  tab[,direction:=fcase(!is.na(padj)&padj<alpha&!is.na(log2FoldChange)&log2FoldChange>=lfc,"Up",
                       !is.na(padj)&padj<alpha&!is.na(log2FoldChange)&log2FoldChange<=-lfc,"Down",
                       default="Not_significant")]
  tab[,direction_noFC:=fcase(!is.na(padj)&padj<alpha&!is.na(log2FoldChange)&log2FoldChange>0,"Up",
                            !is.na(padj)&padj<alpha&!is.na(log2FoldChange)&log2FoldChange<0,"Down",
                            default="Not_significant")]
  setcolorder(tab,c("contrast","gene_id","gene_symbol","baseMean","log2FoldChange","lfcSE",
    "stat","pvalue","padj","direction","is_RBP","hgnc_symbol","entrezgene_id","gene_biotype",
    "RBP_symbol","RBP_description"))
  setorder(tab,padj,-log2FoldChange,na.last=TRUE)
  write_tsv(tab,paste0("results/",contrast_name,"_all_tested_genes.tsv.gz"))
  write_tsv(tab[direction!="Not_significant"],paste0("results/",contrast_name,"_DEGs.tsv"))
  rt <- tab[is_RBP==TRUE]; sig <- rt[direction!="Not_significant"]
  sig0 <- copy(rt[direction_noFC!="Not_significant"]); sig0[,direction:=direction_noFC]
  write_tsv(rt,paste0("results/",contrast_name,"_all_tested_RBP.tsv"))
  write_tsv(sig,paste0("results/",contrast_name,"_DE_RBP.tsv"))
  write_tsv(sig[direction=="Up"],paste0("results/",contrast_name,"_RBP_Up.tsv"))
  write_tsv(sig[direction=="Down"],paste0("results/",contrast_name,"_RBP_Down.tsv"))
  write_tsv(sig0,paste0("results/",contrast_name,"_DE_RBP_noFC.tsv"))
  write_tsv(sig0[direction=="Up"],paste0("results/",contrast_name,"_RBP_Up_noFC.tsv"))
  write_tsv(sig0[direction=="Down"],paste0("results/",contrast_name,"_RBP_Down_noFC.tsv"))
  stopifnot(!anyDuplicated(sig$gene_id),all(sig$padj<alpha),all(abs(sig$log2FoldChange)>=lfc))
  rbp_all[[g]]<-rt; rbp_sig[[g]]<-sig
  rbp_sig_noFC[[g]] <- sig0
  summary_list[[g]] <- data.table(comparison=contrast_name,tumor_n=group_n[condition==g,N],
    normal_n=group_n[condition=="NORMAL",N],RBP_up=sum(sig$direction=="Up"),RBP_down=sum(sig$direction=="Down"),
    RBP_total=nrow(sig),all_gene_up=sum(tab$direction=="Up"),all_gene_down=sum(tab$direction=="Down"),
    tested_RBP=nrow(rt),RBP_padj_NA=sum(is.na(rt$padj)),all_gene_padj_NA=sum(is.na(tab$padj)))
  summary_noFC[[g]] <- data.table(comparison=contrast_name,tumor_n=group_n[condition==g,N],
    normal_n=group_n[condition=="NORMAL",N],RBP_up=sum(sig0$direction=="Up"),RBP_down=sum(sig0$direction=="Down"),
    RBP_total=nrow(sig0),tested_RBP=nrow(rt),RBP_padj_NA=sum(is.na(rt$padj)))
  filtering[[g]]<-data.table(contrast=contrast_name,independent_filter_baseMean=unname(metadata(res)$filterThreshold),
    fdr_alpha=alpha,BH_nonmissing_padj=sum(!is.na(res$padj)),NA_pvalue=sum(is.na(res$pvalue)))
}
summ <- rbindlist(summary_list); sig_all<-rbindlist(rbp_sig); all_rbp<-rbindlist(rbp_all)
summ0 <- rbindlist(summary_noFC)
setorder(summ,-RBP_up,comparison); setorder(summ0,-RBP_up,comparison)
write_tsv(summ,"results/RBP_up_down_summary.tsv")
write_tsv(summ0,"results/RBP_up_down_summary_noFC.tsv")
write_tsv(sig_all,"results/three_contrasts_significant_RBP.tsv")
write_tsv(rbindlist(rbp_sig_noFC),"results/three_contrasts_significant_RBP_noFC.tsv")
write_tsv(rbindlist(filtering),"QC/independent_filtering.tsv")
params <- data.table(parameter=c("date","method","DESeq2_version","design","reference","cohort",
  "subtypes","deduplication","low_count_filter","FDR","FDR_cutoff","FC_cutoff","up_rule","down_rule",
  "FDR_scope","independent_filtering","outlier_count_replacement","Cooks_filter","effect_estimate","model_refit",
  "interpretation","colors","dependence"),value=c("2026-10-26","negative-binomial two-sided Wald test",
  as.character(packageVersion("DESeq2")),"~ condition; one four-group fit","NORMAL",
  "Female primary tumors and solid tissue normal; established project cohort",
  "HER2: HER2+; HR: HER2- and ER+ or PR+; TN: ER-/PR-/HER2-",
  "Highest library size per patient and condition; ties by sample ID","count >= 10 in >= 10 samples",
  "BH-adjusted pvalue (padj)",as.character(alpha),as.character(fc),"padj < 1e-10 AND log2FC >= log2(1.5)",
  "padj < 1e-10 AND log2FC <= -log2(1.5)","Genome-wide per contrast, before RBP intersection",
  "TRUE; results alpha = 1e-10","Disabled (minReplicatesForReplace=Inf)","TRUE",
  "Unshrunken DESeq2 maximum-likelihood log2FC","Fresh DESeq() fit from raw counts; no old checkpoint loaded",
  "Three separate tumor-vs-normal lists; overlapping genes across comparisons are retained",
  "Up #AA4643 (170,70,67); Down #4572A7 (69,114,167); user-specified Excel palette",
  "Unpaired cohort-level contrasts; some patients contribute both normal and tumor; no patient/batch covariates"))
params <- rbind(params,data.table(parameter=c("no_FC_up_rule","no_FC_down_rule","plot_order"),
  value=c("padj < 1e-10 AND log2FC > 0; no effect-size cutoff",
          "padj < 1e-10 AND log2FC < 0; no effect-size cutoff",
          "Descending upregulated RBP count independently in both figures")))
write_tsv(params,"QC/analysis_parameters.tsv")
writeLines(capture.output(sessionInfo()),file.path(out,"QC/R_sessionInfo.txt"))

# Standalone R plotter uses the user-specified red/blue palette and verifies all six bars.
source(file.path(out,"scripts","plot_rbp_counts.R"), local=TRUE)
draw_rbp_counts(out,no_fc=FALSE)
draw_rbp_counts(out,no_fc=TRUE)
print(summ)
print(summ0)
message("Analysis and figure exports complete: ", format(Sys.time()))
