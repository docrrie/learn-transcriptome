suppressPackageStartupMessages({
  library(GEOquery)
  library(hugene10sttranscriptcluster.db)
  library(limma)
  library(FactoMineR)
  library(factoextra)
  library(sva)
  library(patchwork)
})

# get GEO data
read_GEO <- function(series) {
  eSet <- getGEO(series, destdir = 'data/raw', getGPL = FALSE)[[1]]

  expr <- exprs(eSet)
  pd <- pData(eSet)

  expr_anno <- mapIds(
    hugene10sttranscriptcluster.db,
    keys = rownames(expr),
    column = "SYMBOL",
    keytype = "PROBEID",
    multiVals = "first"
  )

  symbol <- unname(expr_anno[rownames(expr)])
  keep <- !is.na(symbol)

  expr_g <- avereps(
    expr[keep, , drop = FALSE],
    ID = symbol[keep]
  )

  return(list(expr_g, pd))
}

GSE_32537 <- read_GEO("GSE32537")
expr_32537 <- GSE_32537[[1]]
pd_32537 <- GSE_32537[[2]]

GSE_110147 <- read_GEO("GSE110147")
expr_110147 <- GSE_110147[[1]]
pd_110147 <- GSE_110147[[2]]


# organize data
common <- intersect(
  rownames(expr_32537),
  rownames(expr_110147)
)

expr_32537_common <- expr_32537[common, , drop = FALSE]
expr_110147_common <- expr_110147[common, , drop = FALSE]

expr_combined <- cbind(
  expr_32537_common,
  expr_110147_common
)

batch <- factor(c(
  rep("GSE32537", ncol(expr_32537_common)),
  rep("GSE110147", ncol(expr_110147_common))
))


lookup <- c(
  "IPF/UIP" = "case",
  "Idiopathic pulmonary fibrosis" = "case",
  "control" = "control",
  "Normal control" = "case"
)

group <- unname(c(
  lookup[pd_32537$`final diagnosis:ch1`],
  lookup[pd_110147$`disease state:ch1`]
))

targets <- data.frame(
  row.names = colnames(expr_combined),
  batch = batch,
  group = group
)

targets <- na.omit(targets)

expr_combined <- expr_combined[,
  colnames(expr_combined) %in% rownames(targets),
  drop = FALSE
]


# generate PCA plots & correct for batch effects
plot_pca <- function(expr, targets) {
  res_pca <- PCA(t(expr), scale.unit = FALSE, graph = FALSE)

  p <- fviz_pca_ind(
    res_pca,
    habillage = targets$`batch`,
    addEllipses = TRUE,
    ellipse.level = 0.95,
    label = "none",
    title = ""
  )

  return(p)
}

p01 <- plot_pca(expr_combined, targets)

expr_combat <- ComBat(
  dat = expr_combined,
  batch = targets$batch,
)

p02 <- plot_pca(expr_combat, targets)


# save objects
saveRDS(expr_combat, file = "data/processed/expr.rds")
saveRDS(targets, file = "data/processed/targets.rds")

saveRDS(p01, file = "data/processed/p01.rds")
saveRDS(p02, file = "data/processed/p02.rds")
