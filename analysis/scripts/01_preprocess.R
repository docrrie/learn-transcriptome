options(check.names = FALSE)

suppressPackageStartupMessages({
  library(GEOquery)
  library(hugene10sttranscriptcluster.db)
  library(limma)
  library(FactoMineR)
  library(factoextra)
  library(sva)
  library(patchwork)
})


read_GEO <- function(series) {
  eSet <- getGEO(series, destdir = 'data/raw', getGPL = FALSE)[[1]]

  exp <- exprs(eSet)
  pd <- pData(eSet)

  exp_anno <- mapIds(
    hugene10sttranscriptcluster.db,
    keys = rownames(exp),
    column = "SYMBOL",
    keytype = "PROBEID",
    multiVals = "first"
  )

  symbol <- unname(exp_anno[rownames(exp)])
  keep <- !is.na(symbol)

  exp_g <- avereps(
    exp[keep, , drop = FALSE],
    ID = symbol[keep]
  )

  return(list(exp_g, pd))
}

GSE_32537 <- read_GEO("GSE32537")
exp_32537 <- GSE_32537[[1]]
pd_32537 <- GSE_32537[[2]]

GSE_110147 <- read_GEO("GSE110147")
exp_110147 <- GSE_110147[[1]]
pd_110147 <- GSE_110147[[2]]


common <- intersect(
  rownames(exp_32537),
  rownames(exp_110147)
)

exp_32537_common <- exp_32537[common, , drop = FALSE]
exp_110147_common <- exp_110147[common, , drop = FALSE]

exp_combined <- cbind(
  exp_32537_common,
  exp_110147_common
)

batch <- factor(c(
  rep("GSE32537", ncol(exp_32537_common)),
  rep("GSE110147", ncol(exp_110147_common))
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

phe <- data.frame(
  row.names = colnames(exp_combined),
  batch = batch,
  group = group
)

phe <- na.omit(phe)

exp_combined <- exp_combined[,
  colnames(exp_combined) %in% rownames(phe),
  drop = FALSE
]


theme_set(theme_classic(base_family = "Arial"))

plot_pca <- function(exp) {
  res_pca <- PCA(t(exp), scale.unit = FALSE, graph = FALSE)
  pct <- res_pca$eig[, 2]

  fviz_pca_ind(
    res_pca,
    habillage = phe$`batch`,
    addEllipses = TRUE,
    ellipse.level = 0.95,
    label = "none",
  ) +
    labs(
      title = "",
      x = sprintf("PCA1 (%.1f%%)", pct[1]),
      y = sprintf("PCA2 (%.1f%%)", pct[2])
    )
}

p1 <- plot_pca(exp_combined)

exp_combat <- ComBat(
  dat = as.matrix(exp_combined),
  batch = phe$batch,
)

p2 <- plot_pca(exp_combat)

combined <- (p1 + p2) +
  plot_annotation(tag_levels = "A") &
  theme(
    plot.tag = element_text(size = 16, face = "bold"),
    plot.tag.position = c(0.02, 0.98)
  )

ggsave(
  filename = "results/figures/01_PCA.tif",
  plot = combined,
  width = 12,
  height = 6,
  units = "in",
  dpi = 300,
  device = "tiff",
  compression = "lzw"
)
