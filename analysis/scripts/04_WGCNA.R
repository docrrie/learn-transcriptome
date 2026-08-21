suppressPackageStartupMessages({
  library(WGCNA)
})

enableWGCNAThreads()
set.seed(2026)

expr <- readRDS("data/processed/expr.rds")
pd <- readRDS("data/processed/pd.rds")


# preprocess
expr <- expr[rowMeans(expr) > quantile(rowMeans(expr), 0.05), ]

gene_mad <- apply(expr, 1, mad)
keep <- names(sort(gene_mad, decreasing = TRUE))[1:5000]
expr <- t(expr[keep, ])

pd <- data.frame(
  case = as.integer(pd$group == "case"),
  control = as.integer(pd$group == "control"),
  row.names = rownames(pd)
)

# QC
# gsg <- goodSamplesGenes(expr, verbose = 3)

# sample tree
sample_tree <- hclust(dist(expr), method = "average")
trait_colors <- numbers2colors(pd, signed = FALSE)

plotDendroAndColors(
  sample_tree,
  trait_colors,
  groupLabels = names(pd),
  main = ""
)

# pick soft threshold
# sft <- pickSoftThreshold(expr, powerVector = 1:30, networkType = "unsigned")

# fit <- sft$fitIndices

# plot(
#   fit[, "Power"],
#   -sign(fit[, "slope"]) * fit[, "SFT.R.sq"],
#   xlab = "Soft Threshold (power)",
#   ylab = "Scale Free Topology Model Fit, signed R^2",
#   type = "n",
#   main = "Scale independence"
# )
# text(
#   fit[, "Power"],
#   -sign(fit[, "slope"]) * fit[, "SFT.R.sq"],
#   labels = 1:30,
#   col = "red"
# )
# abline(h = 0.90, col = "red")

# # C：Mean connectivity
# plot(
#   fit[, "Power"],
#   fit[, "mean.k."],
#   xlab = "Soft Threshold (power)",
#   ylab = "Mean Connectivity",
#   type = "n",
#   main = "Mean connectivity"
# )
# text(
#   fit[, "Power"],
#   fit[, "mean.k."],
#   labels = 1:30,
#   col = "red"
# )

# network construction
net <- blockwiseModules(
  expr,
  power = 5,
  TOMType = "unsigned",
  networkType = "unsigned",
  minModuleSize = 30,
  deepSplit = 2,
  mergeCutHeight = 0.25,
  numericLabels = FALSE,
  pamRespectsDendro = FALSE
)


# heatmap
MEs <- orderMEs(net$MEs)
cor_ME <- cor(MEs, pd, use = "p")
p_ME <- corPvalueStudent(cor_ME, nrow(expr))

text_matrix <- paste0(
  signif(cor_ME, 2),
  "\n(P = ",
  signif(p_ME, 2),
  ")"
)
dim(text_matrix) <- dim(cor_ME)

labeledHeatmap(
  Matrix = cor_ME,
  xLabels = colnames(pd),
  yLabels = colnames(MEs),
  ySymbols = colnames(MEs),
  yColorLabels = TRUE,
  colors = blueWhiteRed(50),
  textMatrix = text_matrix,
  setStdMargins = FALSE,
  cex.text = 0.75,
  zlim = c(-1, 1)
)
