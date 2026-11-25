#!/usr/bin/env Rscript
args<-commandArgs(trailingOnly=TRUE);positional<-args[!grepl("^--",args)]
root<-normalizePath(if(length(positional)) positional[1] else getwd())
fc_arg<-grep("^--fc=",args,value=TRUE)
fc_cut<-if(length(fc_arg)) as.numeric(sub("^--fc=","",fc_arg[1])) else as.numeric(Sys.getenv("RBP_FIGURE_FC","1.5"))
use_stars<-"--stars" %in% args || Sys.getenv("RBP_FIGURE_STARS","0")=="1"
stopifnot(fc_cut %in% c(1.5,2))
tag<-if(fc_cut==1.5) "FC1p5" else "FC2"
out<-file.path(root,"outputs/TNBC_DESeq2_RBP_20261026_FC1p5")
inference<-file.path(root,"outputs/TNBC_RBP_proportion_inference_20261027")
.libPaths(c(file.path(root,"TCGA/.Rlib"),.libPaths()))
suppressPackageStartupMessages({library(data.table);library(ggplot2)})
setDTthreads(1L)
font_pt<-14
summ<-fread(if(fc_cut==1.5) file.path(out,"results/RBP_up_down_summary.tsv") else
             file.path(inference,"results/FC2_RBP_up_down_summary.tsv"))
groups<-c("TN","HER2","HR")
summ<-summ[match(paste0(groups,"_vs_NORMAL"),comparison)]
stopifnot(!anyNA(summ$comparison),all(summ$RBP_total==summ$RBP_up+summ$RBP_down))
dt<-melt(summ[,.(comparison,RBP_up,RBP_down,RBP_total)],
  id.vars=c("comparison","RBP_total"),variable.name="direction",value.name="count")
dt[,fraction:=count/RBP_total]
dt[,label:=sprintf("%.1f%%",100*fraction)]
dt[,direction:=factor(direction,levels=c("RBP_up","RBP_down"),labels=c("Up","Down"))]
dt[,subtype:=factor(comparison,levels=paste0(groups,"_vs_NORMAL"),labels=groups)]
setorder(dt,subtype,direction)
palette<-c(Up="#AA4643",Down="#4572A7")
tests<-fread(file.path(inference,"results",paste0(tag,"_pairwise_patient_bootstrap.tsv")))
tests<-tests[match(c("TN vs HER2","HER2 vs HR","TN vs HR"),comparison)]
stopifnot(!anyNA(tests$comparison),all(tests$B==19999L))
tests[,`:=`(x=match(group1,groups),xend=match(group2,groups),
             y=c(1.055,1.155,1.255),label=sprintf("P = %.5f",p_Holm))]
if(use_stars) tests[,label:=fcase(p_Holm<.001,"***",p_Holm<.01,"**",p_Holm<.05,"*",default="ns")]
tips<-rbind(tests[,.(x,xend=x,y,yend=y-.018)],tests[,.(x=xend,xend,y,yend=y-.018)])
p<-ggplot(dt,aes(subtype,fraction,fill=direction))+
  geom_col(width=.58,position=position_stack(reverse=TRUE))+
  geom_text(aes(label=label),position=position_stack(vjust=.5,reverse=TRUE),
            family="Arial",size=font_pt,size.unit="pt",color="white")+
  geom_segment(data=tests,aes(x=x,xend=xend,y=y,yend=y),inherit.aes=FALSE,
               color="black",linewidth=.35)+
  geom_segment(data=tips,aes(x=x,xend=xend,y=y,yend=yend),inherit.aes=FALSE,
               color="black",linewidth=.35)+
  geom_text(data=tests,aes(x=(x+xend)/2,y=y+.008,label=label),inherit.aes=FALSE,
            family="Arial",size=font_pt,size.unit="pt",color="black",vjust=0)+
  scale_fill_manual(values=palette,breaks=c("Up","Down"))+
  scale_y_continuous(limits=c(0,1.35),breaks=c(0,.25,.5,.75,1),
                     labels=scales::label_percent(accuracy=1),
                     expand=expansion(mult=c(0,0)),guide=guide_axis(cap="both"))+
  scale_x_discrete(expand=expansion(add=.55))+
  labs(x=NULL,y=NULL,fill=NULL)+
  theme_classic(base_size=font_pt,base_family="Arial")+
  theme(text=element_text(size=font_pt,family="Arial"),
        axis.text=element_text(size=font_pt,color="black"),
        axis.text.x=element_text(size=font_pt,margin=margin(t=4)),
        axis.line=element_line(linewidth=.35),axis.ticks=element_line(linewidth=.35),
        axis.ticks.x=element_blank(),
        legend.position="top",legend.justification="center",legend.text=element_text(size=font_pt),
        legend.key.size=grid::unit(12,"pt"),legend.spacing.x=grid::unit(2,"mm"),
        legend.margin=margin(0,0,3,0),
        plot.margin=margin(5,5,4,4))
figure_dir<-Sys.getenv("RBP_FIGURE_DIR",file.path(out,"figures"))
figure_prefix<-Sys.getenv("RBP_FIGURE_PREFIX",paste0("RBP_percent_stacked_",tag,"_Holm_14pt"))
dir.create(figure_dir,recursive=TRUE,showWarnings=FALSE)
stem<-file.path(figure_dir,figure_prefix)
ragg::agg_png(paste0(stem,".png"),width=115,height=115,units="mm",res=300);print(p);dev.off()
ragg::agg_tiff(paste0(stem,".tiff"),width=115,height=115,units="mm",res=600,compression="lzw");print(p);dev.off()
svglite::svglite(paste0(stem,".svg"),width=115/25.4,height=115/25.4);print(p);dev.off()
grDevices::cairo_pdf(paste0(stem,".pdf"),width=115/25.4,height=115/25.4,family="Arial");print(p);dev.off()
# Explicitly check preservation of source proportions and all annotation mappings.
if(fc_cut==1.5) {
  original<-fread(file.path(out,"figures/RBP_percent_stacked_FC1p5_source_data.tsv"))
  stopifnot(identical(as.character(dt$subtype),as.character(original$subtype)),
            identical(as.character(dt$direction),as.character(original$direction)),
            identical(dt$count,original$count),max(abs(dt$fraction-original$fraction))<1e-12)
}
for(k in seq_along(groups)) {
  tab<-fread(file.path(out,"results",paste0(groups[k],"_vs_NORMAL_all_tested_RBP.tsv")))
  keep<-!is.na(tab$padj)&tab$padj<1e-10
  stopifnot(sum(keep&tab$log2FoldChange>=log2(fc_cut),na.rm=TRUE)==summ$RBP_up[k],
            sum(keep&tab$log2FoldChange<=-log2(fc_cut),na.rm=TRUE)==summ$RBP_down[k])
}
bars<-ggplot_build(p)$data[[1]]
stopifnot(nrow(bars)==6L,all(abs(bars$ymax-bars$ymin-dt$fraction)<1e-12),
          identical(toupper(bars$fill),unname(toupper(palette[as.character(dt$direction)]))),
          all(tests$x<tests$xend))
fwrite(tests[,.(comparison,group1,group2,p_Holm,label,x,xend,y)],
       paste0(stem,"_annotations.tsv"),sep="\t")
fwrite(dt[,.(comparison,subtype,direction,count,total=RBP_total,fraction,label)],
       paste0(stem,"_source_data.tsv"),sep="\t")
fwrite(data.table(parameter=c("FC_cutoff","gene_FDR_threshold","bootstrap_resamples","multiple_testing","group_order","text_size_pt","canvas_mm"),
  value=c(as.character(fc_cut),"1e-10","19999","Holm","TN,HER2,HR","14","115")),
  paste0(stem,"_plot_parameters.tsv"),sep="\t")
print(tests[,.(comparison,p_Holm,label)])
cat("Verified: original proportions and colors; three Holm annotations; 14 pt; square canvas\n")
