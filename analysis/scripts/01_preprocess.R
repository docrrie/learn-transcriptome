suppressPackageStartupMessages({
  library(GEOquery) # get GEO data
  library(hugene10sttranscriptcluster.db) # annotate probes
  library(limma) # avereps(), normalizeBetweenArrays()
  library(FactoMineR) # PCA
  library(factoextra) # PCA plot
  library(sva) # ComBat
})

# get GEO data
read_GEO <- function(series) {
  eSet <- getGEO(series, destdir = 'data/raw', getGPL = FALSE)[[1]]

  expr <- exprs(eSet)
  pd <- pData(eSet)

  probe_symbol <- mapIds(
    hugene10sttranscriptcluster.db,
    keys = rownames(expr),
    column = "SYMBOL",
    keytype = "PROBEID",
    multiVals = "first"
  )

  symbol <- unname(probe_symbol[rownames(expr)])
  keep <- !is.na(symbol)
  expr_anno <- expr[keep, , drop = FALSE]
  rownames(expr_anno) <- unname(probe_symbol[keep])

  expr_g <- avereps(
    expr_anno,
    ID = rownames(expr_anno)
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
  "Normal control" = "control"
)

group <- unname(c(
  lookup[pd_32537$`final diagnosis:ch1`],
  lookup[pd_110147$`disease state:ch1`]
))

sample_info <- data.frame(
  row.names = colnames(expr_combined),
  batch = factor(batch, levels = c("GSE32537", "GSE110147")),
  group = factor(group, levels = c("control", "case"))
)

sample_info <- na.omit(sample_info)

expr_combined <- expr_combined[,
  colnames(expr_combined) %in% rownames(sample_info),
  drop = FALSE
]


# correction
expr_normalized <- normalizeBetweenArrays(
  expr_combined,
  method = "quantile"
)

mod <- model.matrix(~group, data = sample_info)

expr_combat <- ComBat(
  dat = expr_normalized,
  batch = sample_info$batch,
  mod = mod
)


# generate PCA plots
plot_pca <- function(expr, sample_info) {
  res_pca <- PCA(t(expr), scale.unit = FALSE, graph = FALSE)

  p <- fviz_pca_ind(
    res_pca,
    habillage = sample_info$`batch`,
    addEllipses = TRUE,
    ellipse.level = 0.95,
    label = "none",
    title = ""
  )

  return(p)
}

p01 <- plot_pca(expr_combined, sample_info)
p02 <- plot_pca(expr_combat, sample_info)


# save objects
saveRDS(expr_combat, file = "data/processed/expr.rds")
saveRDS(sample_info["group"], file = "data/processed/pd.rds")
write.csv(sample_info, file = "results/tables/sample_info.csv")

saveRDS(p01, file = "data/processed/p01_PCA_before.rds")
saveRDS(p02, file = "data/processed/p02_PCA_after.rds")
