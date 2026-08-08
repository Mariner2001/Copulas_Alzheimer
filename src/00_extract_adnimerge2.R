#!/usr/bin/env Rscript
# Extrae de ADNIMERGE2_*.tar.gz las tres tablas que no vienen como CSV suelto.
# Uso:  Rscript src/00_extract_adnimerge2.R
# Salida: data/adsl.csv, data/arm.csv, data/uwnp.csv
tgz <- list.files("data", pattern = "^ADNIMERGE2", full.names = TRUE)
if (length(tgz) == 0) stop("No encuentro ADNIMERGE2_*.tar.gz en data/")
cat("Usando:", tgz[1], "\n")
untar(tgz[1], files = c("ADNIMERGE2/data/ADSL.rda",
                        "ADNIMERGE2/data/ARM.rda",
                        "ADNIMERGE2/data/UWNPSYCHSUM.rda"), exdir = "data/_tmp")
suppressWarnings({
  load("data/_tmp/ADNIMERGE2/data/ADSL.rda")
  load("data/_tmp/ADNIMERGE2/data/ARM.rda")
  load("data/_tmp/ADNIMERGE2/data/UWNPSYCHSUM.rda")
})
a <- as.data.frame(ADSL); a$APOE <- as.character(a$APOE)
write.csv(a[, c("SUBJID","ORIGPROT","AGE","SEX","EDUC","DX","APOE")],
          "data/adsl.csv", row.names = FALSE, fileEncoding = "UTF-8")
write.csv(as.data.frame(ARM)[, c("RID","ORIGPROT","COLPROT","ARM","ENROLLED")],
          "data/arm.csv", row.names = FALSE)
write.csv(as.data.frame(UWNPSYCHSUM)[, c("RID","VISCODE2","COLPROT","ORIGPROT","ADNI_MEM","ADNI_EF")],
          "data/uwnp.csv", row.names = FALSE)
unlink("data/_tmp", recursive = TRUE)
cat("OK -> data/adsl.csv, data/arm.csv, data/uwnp.csv\n")
