#!/usr/bin/env Rscript
# GO/KEGG enrichment (g:Profiler) for the DESeq2 contrasts. Arabidopsis thaliana.
suppressMessages({library(gprofiler2); library(ggplot2)})
OUT <- path.expand("~/deseq2/go"); dir.create(OUT, showWarnings=FALSE, recursive=TRUE)
at <- function(tx) unique(na.omit(regmatches(tx, regexpr("AT[1-5CMG]G[0-9]{5}", tx, ignore.case=TRUE))))
do_go <- function(csv, name){
  csv <- path.expand(csv)
  d <- read.csv(csv, row.names=1)
  d <- d[!is.na(d$padj) & d$padj<0.05 & abs(d$log2FoldChange)>1, ]
  genes <- at(rownames(d))
  cat(name, "-> DE transcripts:", nrow(d), "| AT genes:", length(genes), "\n")
  if (length(genes) < 5) { cat("  too few genes, skip\n"); return(invisible()) }
  g <- tryCatch(gost(genes, organism="athaliana", significant=TRUE,
                     sources=c("GO:BP","GO:MF","GO:CC","KEGG")), error=function(e){cat("  gost err:",conditionMessage(e),"\n");NULL})
  if (is.null(g) || is.null(g$result) || nrow(g$result)==0) { cat("  no enrichment\n"); return(invisible()) }
  write.csv(apply(g$result, 2, as.character), file.path(OUT, paste0(name,"_GO.csv")), row.names=FALSE)
  r <- g$result[order(g$result$p_value), ]; r <- head(r, 15)
  r$term_name <- factor(r$term_name, levels=rev(r$term_name))
  p <- ggplot(r, aes(-log10(p_value), term_name, fill=source)) + geom_col() +
    labs(title=paste0("GO/KEGG enrichment - ", name), subtitle=paste0(length(genes)," genes"), x="-log10(p-value)", y=NULL) +
    theme_bw(base_size=14) + theme(plot.title=element_text(face="bold"), legend.position="bottom")
  ggsave(file.path(OUT, paste0(name,"_GO.png")), p, width=11, height=6.5, dpi=320)
  ggsave(file.path(OUT, paste0(name,"_GO.pdf")), p, width=11, height=6.5)
  cat("  ", nrow(g$result), "terms ->", name, "\n")
}
do_go("~/deseq2/out/D5_volcano_interaction.csv", "interaction_ANAC017dep")
do_go("~/deseq2/out/D3_volcano_AA_in_WT.csv", "AA_in_WT")
do_go("~/deseq2/out/D4_volcano_genotype.csv", "genotype_anac017_vs_WT")
cat("GO_DONE\n")
