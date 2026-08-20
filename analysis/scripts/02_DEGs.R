suppressPackageStartupMessages({
  library(limma)
  library(EnhancedVolcano)
})

expr <- readRDS("data/processed/expr.rds")
pd <- readRDS("data/processed/pd.rds")


# differential expression analysis
design <- model.matrix(~ 0 + group, data = pd)
colnames(design) <- levels(pd$group)

fit <- lmFit(expr, design)

contrast.matrix <- makeContrasts(
  case - control,
  levels = design
)

fit2 <- contrasts.fit(fit, contrast.matrix)
fit2 <- eBayes(fit2, robust = TRUE)

res_DEG <- topTable(
  fit2,
  number = Inf,
  adjust.method = "BH",
  sort.by = "P"
)

res_DEG_sig <- subset(
  res_DEG,
  adj.P.Val < 0.05 & abs(logFC) >= log2(1.3)
)

cat("DEGs:", nrow(res_DEG_sig), "\n")
cat("up:", sum(res_DEG_sig$logFC > 0), "\n")
cat("down:", sum(res_DEG_sig$logFC < 0), "\n")


# volcano plot
EnhancedVolcano(
  res_DEG,

  lab = rep("", nrow(res_DEG)),

  x = "logFC",
  y = "adj.P.Val",

  pCutoff = 0.05,
  FCcutoff = 1.3,

  title = "",
  subtitle = "",
  caption = "",

  xlab = "log2 Fold Change",
  ylab = "-log10 adjusted P-value",

  legendLabels = c("", "", ""),
  legendPosition = "right",
  legendIconSize = 4
)


# save objects
write.csv(
  res_DEG,
  "results/tables/res_DEG_all.csv",
)

write.csv(
  res_DEG_sig,
  "results/tables/DEGs.csv",
)
