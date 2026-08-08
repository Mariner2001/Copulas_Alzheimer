#!/usr/bin/env Rscript
# =============================================================================
#  Benchmark GCGM (Gaussian copula graphical model) sobre los MISMOS datos
#  que el vine.  Es la comparacion controlada del articulo.
#
#  Uso:   Rscript gcgm_benchmark.R
#  Entrada:  data/adni_merged.csv   (generado por 01_merge.py)
#  Salida:   output/gcgm_results.csv, output/gcgm_summary.txt
# =============================================================================

if (!requireNamespace("BDgraph", quietly = TRUE))
  install.packages("BDgraph", repos = "https://cloud.r-project.org")
library(BDgraph)

set.seed(20240607)
d <- read.csv("data/adni_merged.csv")

VARS <- c("AGE","SEX_bin","EDUC","ADNI_MEM","ADNI_EF","APOE4_bin","AMY","FDG",
          "Hippocampus","Caudate","Putamen","Thalamus","PCC","Precuneus")
LAB  <- c("Age","Sex","Education","Memory","ExecFun","APOE4","Amyloid","FDG",
          "Hippocampus","Caudate","Putamen","Thalamus","PCC","Precuneus")

X <- d[, VARS]
X$AMY <- -X$AMY                       # misma orientacion que el vine: menor = mas patologia
X <- as.matrix(X)
colnames(X) <- LAB

# indices (1-based) de las variables NO continuas: sexo, educacion, APOE4
not.cont <- rep(0L, length(VARS))
not.cont[match(c("SEX_bin","EDUC","APOE4_bin"), VARS)] <- 1L

cat("Ajustando GCGM:  n =", nrow(X), " p =", ncol(X), "\n")
fit <- bdgraph(data = X, method = "gcgm", not.cont = not.cont,
               iter = 30000, burnin = 10000, g.prior = 0.2, save = TRUE)

P   <- plinks(fit)                    # probabilidades de inclusion a posteriori
K   <- fit$K_hat                      # matriz de precision estimada
PC  <- -K / sqrt(outer(diag(K), diag(K)))   # correlaciones parciales
diag(PC) <- 1

# ---- tabla larga de todos los pares ----
res <- do.call(rbind, lapply(1:(ncol(X)-1), function(i)
  do.call(rbind, lapply((i+1):ncol(X), function(j)
    data.frame(a = LAB[i], b = LAB[j], EIP = P[i,j], pcor = PC[i,j])))))
res <- res[order(-res$EIP), ]
dir.create("output", showWarnings = FALSE)
write.csv(res, "output/gcgm_results.csv", row.names = FALSE)

# ---- resumen que reproduce las cifras citadas en el articulo ----
n_edges <- sum(res$EIP > 0.5)
am      <- res[(res$a=="Age" & res$b=="Memory") | (res$a=="Memory" & res$b=="Age"), ]

sink("output/gcgm_summary.txt")
cat("GCGM benchmark\n==============\n")
cat("n =", nrow(X), "  p =", ncol(X), "  pares =", nrow(res), "\n")
cat("Aristas con EIP > 0.5 :", n_edges, "   <- articulo: 33\n\n")
cat("Arista edad-memoria:\n"); print(am)
cat("   <- articulo: EIP ~ 0.01, correlacion parcial ~ 0\n\n")
cat("15 pares con mayor EIP:\n"); print(head(res, 15), row.names = FALSE)
sink()

cat("\nHECHO. Ver output/gcgm_summary.txt\n")
cat("Aristas con EIP>0.5:", n_edges, "\n")

# -----------------------------------------------------------------------------
# SOLAPAMIENTO CON EL VINE
# Ejecutar despues de 02_bootstrap.py + 03_aggregate.py; requiere
# output/bootstrap_B1000.csv, que lista las aristas activas del vine.
# -----------------------------------------------------------------------------
if (file.exists("output/vine_edges.csv")) {
  ve <- read.csv("output/vine_edges.csv")
  key <- function(a,b) apply(cbind(a,b), 1, function(z) paste(sort(z), collapse="|"))
  gk <- key(res$a, res$b)[res$EIP > 0.5]
  vk <- key(ve$a, ve$b)
  cat("\nCompartidas :", length(intersect(gk, vk)),  "  <- articulo: 19\n")
  cat("Solo vine   :", length(setdiff(vk, gk)),      "  <- articulo: 6\n")
  cat("Solo GCGM   :", length(setdiff(gk, vk)),      "  <- articulo: 14\n")
}
