suppressPackageStartupMessages({
  library(clusterProfiler)
  library(org.Hs.eg.db)
  library(ggplot2)
})

DEG_all <- readRDS("data/processed/DEG_all.rds")
DEG_sig <- readRDS("data/processed/DEG_sig.rds")

up <- rownames(subset(DEG_sig, change == "up"))
down <- rownames(subset(DEG_sig, change == "down"))


run_ORA <- function(DEG_list, bg_list) {
  DEG_entrez <- bitr(
    DEG_list,
    fromType = "SYMBOL",
    toType = "ENTREZID",
    OrgDb = org.Hs.eg.db
  )$ENTREZID

  bg_entrez <- bitr(
    bg_list,
    fromType = "SYMBOL",
    toType = "ENTREZID",
    OrgDb = org.Hs.eg.db
  )$ENTREZID

  res_GO <- enrichGO(
    gene = DEG_entrez,
    universe = bg_entrez,
    OrgDb = org.Hs.eg.db,
    ont = "ALL",
    pvalueCutoff = 0.05,
    qvalueCutoff = 1,
    readable = TRUE
  )

  res_KEGG <- enrichKEGG(
    gene = DEG_entrez,
    universe = bg_entrez,
    pvalueCutoff = 0.05,
    qvalueCutoff = 1
  )

  dotplot(
    res_up[[1]],
    showCategory = 10,
    orderBy = "GeneRatio",
    split = "ONTOLOGY",
  ) +
    facet_grid(ONTOLOGY ~ ., scales = 'free_y')

  dotplot(
    res_up[[2]],
    showCategory = 10,
    orderBy = "GeneRatio"
  )

  return(list(res_GO, res_KEGG))
}

res_up <- run_ORA(up, rownames(DEG_all))

res_down <- run_ORA(down, rownames(DEG_all))
