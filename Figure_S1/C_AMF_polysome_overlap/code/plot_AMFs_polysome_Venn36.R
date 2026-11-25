#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(ragg)
  library(svglite)
  library(readxl)
})
data.table::setDTthreads(1L)
options(stringsAsFactors=FALSE, width=140)
compact_mode <- '--compact' %in% commandArgs(trailingOnly=TRUE)
script_arg <- grep('^--file=', commandArgs(FALSE), value=TRUE)
script_path <- normalizePath(sub('^--file=','',script_arg[1]), mustWork=TRUE)
out_dir <- dirname(script_path)
root <- dirname(dirname(out_dir))
for(d in c('tables','QC')) dir.create(file.path(out_dir,d),recursive=TRUE,showWarnings=FALSE)
inputs <- c(
  AMF_catalogue=file.path(root,'08polysome','polysome_up_502_TN131_AMFs215_Venn_20260902','source_data','AMFs_215_source_rows.tsv'),
  FOT_limma=file.path(root,'08polysome','FOT_pipeline','results','tables','FOT_limma_all_tested.csv'),
  manuscript_AMF36=file.path(root,'outputs','01a1773c-9bfe-70a2-bd0c-558728c6497d','AMFs36_FOT_FDR_sorted_checked.xlsx')
)
stopifnot(all(file.exists(inputs)))
clean_symbol <- function(x) toupper(trimws(as.character(x)))
catalogue <- fread(inputs[['AMF_catalogue']])
stats <- fread(inputs[['FOT_limma']])
reference36 <- read_excel(inputs[['manuscript_AMF36']],sheet='AMFs_36_FDR_sorted')
catalogue[,gene_symbol:=clean_symbol(gene_symbol_clean)]
stats[,gene_symbol:=clean_symbol(Gene_name)]
stopifnot(nrow(stats)==2674,uniqueN(catalogue$gene_symbol)==215,
          !anyNA(catalogue$gene_symbol),!any(catalogue$gene_symbol==''))
amf_rows <- catalogue[gene_symbol!='RPL7']
amfs <- sort(unique(amf_rows$gene_symbol))
up_rows <- stats[!is.na(BH_FDR) & BH_FDR<0.05 &
                  !is.na(log2FC_Polysome_vs_Free) & log2FC_Polysome_vs_Free>0]
up <- sort(unique(up_rows$gene_symbol))
shared <- sort(intersect(amfs,up))
amf_only <- sort(setdiff(amfs,up))
up_only <- sort(setdiff(up,amfs))
stopifnot(length(amfs)==214,length(up)==291,length(shared)==36,
          length(amf_only)==178,length(up_only)==255,
          nrow(up_rows)==length(up),setequal(shared,reference36$Gene),
          'RPL7' %in% up,!('RPL7' %in% amfs))

# Recompute the original BH adjustment from all tested ordinary-limma P values.
stopifnot(isTRUE(all.equal(p.adjust(stats$P_Value,method='BH'),stats$BH_FDR,tolerance=1e-12)))
membership <- data.table(gene_symbol=sort(union(amfs,up)))
membership[,`:=`(AMF=gene_symbol %in% amfs,Polysome_upregulated=gene_symbol %in% up)]
membership[,region:=fifelse(AMF & Polysome_upregulated,'Intersection',fifelse(AMF,'AMF only','Polysome-up only'))]
stopifnot(nrow(membership)==469,sum(membership$AMF)==214,
          sum(membership$Polysome_upregulated)==291,
          sum(membership$AMF & membership$Polysome_upregulated)==36)
region_counts <- data.table(region=c('AMF only','Intersection','Polysome-up only'),gene_count=c(178,36,255))
set_summary <- data.table(
  metric=c('Original AMF catalogue','Excluded AMF/RP overlap RPL7','AMFs after RP exclusion',
           'All proteins tested by FOT limma','Polysome-upregulated proteins','Shared AMFs',
           'AMF only','Polysome-up only','Union'),
  gene_count=c(215,1,214,2674,291,36,178,255,469))
shared_detail <- merge(
  up_rows[gene_symbol %in% shared,.(gene_symbol,protein_accession=Protein_accession,
    unique_peptides=Unique_peptides,Free_n=n_observed_Free,Polysome_n=n_observed_Polysome,
    log2FC_Polysome_vs_Free,FC_Polysome_vs_Free=2^log2FC_Polysome_vs_Free,
    ordinary_P=P_Value,ordinary_BH_FDR=BH_FDR)],
  amf_rows[,.(gene_symbol,AMF_function=functional_annotation,core_AMF)],by='gene_symbol')
setorder(shared_detail,ordinary_BH_FDR)
shared_detail[,FDR_rank:=seq_len(.N)]
setcolorder(shared_detail,c('FDR_rank',setdiff(names(shared_detail),'FDR_rank')))
write_table <- function(dt,name) fwrite(dt,file.path(out_dir,'tables',paste0(name,'.tsv')),sep='\t',na='')
write_table(set_summary,'set_summary')
write_table(region_counts,'Venn_region_counts')
write_table(membership,'all_gene_membership_469')
write_table(amf_rows,'AMFs_214')
write_table(up_rows,'polysome_up_proteins_291')
write_table(shared_detail,'intersection_AMFs_36_FDR_sorted')
write_table(data.table(gene_symbol=amf_only),'AMF_only_178')
write_table(data.table(gene_symbol=up_only),'polysome_up_only_255')

make_circle <- function(cx,group) {
  theta <- seq(0,2*pi,length.out=1200)
  data.frame(x=cx+1.8*cos(theta),y=1.8*sin(theta),group=group)
}
circles <- rbind(make_circle(-1.05,'AMF'),make_circle(1.05,'Polysome'))
p <- ggplot()+
  geom_polygon(data=circles,aes(x,y,group=group,fill=group),alpha=.25,colour=NA)+
  geom_path(data=circles,aes(x,y,group=group,colour=group),linewidth=.65)+
  scale_fill_manual(values=c(AMF='#4572A7',Polysome='#AA4643'))+
  scale_colour_manual(values=c(AMF='#4572A7',Polysome='#AA4643'))+
  annotate('text',x=-1.65,y=if(compact_mode)2.47 else 2.35,label='AMFs',family='Arial',fontface='bold',size=if(compact_mode)18/ggplot2::.pt else 5,colour='#345D8B')+
  annotate('text',x=1.65,y=if(compact_mode)2.47 else 2.35,label='Polysome-upregulated',family='Arial',fontface='bold',size=if(compact_mode)18/ggplot2::.pt else 4.8,colour='#91413F')+
  annotate('text',x=-1.65,y=if(compact_mode)2.06 else 1.96,label='n = 214',family='Arial',size=if(compact_mode)16/ggplot2::.pt else 4.4,colour='#345D8B')+
  annotate('text',x=1.65,y=if(compact_mode)2.06 else 1.96,label='n = 291',family='Arial',size=if(compact_mode)16/ggplot2::.pt else 4.4,colour='#91413F')+
  annotate('text',x=-1.65,y=0,label=length(amf_only),family='Arial',fontface='bold',size=10,colour='#263E59')+
  annotate('text',x=0,y=0,label=length(shared),family='Arial',fontface='bold',size=11,colour='#29212D')+
  annotate('text',x=1.65,y=0,label=length(up_only),family='Arial',fontface='bold',size=10,colour='#66302F')+
  coord_fixed(xlim=if(compact_mode)c(-3.1,3.38) else c(-3.35,3.35),
              ylim=if(compact_mode)c(-1.92,2.78) else c(-3.35,2.85),expand=FALSE,clip='off')+
  theme_void(base_family='Arial')+
  theme(legend.position='none',plot.margin=margin(5,5,3,5),plot.background=element_rect(fill='white',colour=NA))
if(!compact_mode) p <- p+
  annotate('text',x=0,y=-2.07,label='36 polysome-enriched AMFs',family='Arial',fontface='bold',size=4.5)+
  annotate('text',x=0,y=-2.47,label='AMF catalogue excludes RPL7 (RP).',family='Arial',size=3.3,colour='#555555')+
  annotate('text',x=0,y=-2.79,label='Polysome vs Free; ordinary limma BH FDR < 0.05, log2FC > 0',family='Arial',size=3.1,colour='#666666')+
  annotate('text',x=0,y=-3.10,label='Schematic circles; areas are not proportional to set sizes.',family='Arial',size=2.8,colour='#777777')

# Use a compact near-square canvas so circle geometry stays circular.
width_mm <- if(compact_mode)170 else 150
height_mm <- if(compact_mode)126 else 140
stem <- file.path(out_dir,if(compact_mode)'AMFs_vs_polysome_up_Venn36_compact' else 'AMFs_vs_polysome_up_Venn36')
if(!compact_mode) {
  svglite(paste0(stem,'.svg'),width=width_mm/25.4,height=height_mm/25.4);print(p);dev.off()
}
cairo_pdf(paste0(stem,'.pdf'),width=width_mm/25.4,height=height_mm/25.4,family='Arial');print(p);dev.off()
agg_png(paste0(stem,'.png'),width=width_mm,height=height_mm,units='mm',res=600,background='white');print(p);dev.off()
if(!compact_mode) {
  agg_tiff(paste0(stem,'.tiff'),width=width_mm,height=height_mm,units='mm',res=600,compression='lzw',background='white');print(p);dev.off()
}
agg_png(file.path(out_dir,'QC',if(compact_mode)'Venn36_compact_preview.png' else 'Venn36_preview.png'),width=900,height=round(900*height_mm/width_mm),units='px',res=900/(width_mm/25.4),background='white');print(p);dev.off()
fwrite(data.table(parameter=c("AMF_set_size","polysome_up_set_size","intersection_size","AMF_only","polysome_only","union_size"),
 value=c(length(amfs),length(up),length(shared),length(amf_only),length(up_only),length(union(amfs,up)))),
 file.path(out_dir,"analysis_parameters.tsv"),sep="\t")
fwrite(data.table(source=names(inputs),md5=unname(tools::md5sum(inputs))),file.path(out_dir,'QC','input_md5.tsv'),sep='\t')
fwrite(data.table(check=c('AMF_214_unique_symbols','Polysome_up_291_unique_symbols',
  '36_intersection_matches_manuscript_workbook','178_plus_36_equals_214',
  '255_plus_36_equals_291','178_plus_36_plus_255_equals_469','original_BH_FDR_independently_verified'),result='PASS'),
  file.path(out_dir,'QC','data_validation.tsv'),sep='\t')
writeLines(capture.output(sessionInfo()),file.path(out_dir,'QC','R_sessionInfo.txt'))
cat('AMFs =',length(amfs),'; Polysome-up =',length(up),'; intersection =',length(shared),
    '; exclusive regions =',length(amf_only),'and',length(up_only),'\n')
