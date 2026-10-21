#!/usr/bin/env Rscript
# Measured snATAC per site split by experimental condition (YS, MAS, MAE) for the cardiomyocyte and fibroblast pseudobulks.
# Used to validate the age / exercise class labels with independent accessibility data. Output: Data/derived/measured_by_condition.csv
# Usage: Rscript 04_measured_by_condition.R [data_dir]
suppressPackageStartupMessages({library(data.table); library(Matrix)})
args <- commandArgs(trailingOnly = TRUE); data_dir <- if (length(args) >= 1) args[1] else "Data"
x <- readRDS(file.path(data_dir, "heart_snatac_pseudobulk_raw_counts.rds"))
lib <- colSums(x); cell_of <- sub("--.*", "", colnames(x)); samp <- sub(".*--", "", colnames(x))
cond <- ifelse(grepl("^YS", samp), "YS", ifelse(grepl("^MAS", samp), "MAS", ifelse(grepl("^MAE", samp), "MAE", NA)))
cat("samples per condition (CM):", paste(names(table(cond[cell_of == "Cardiomyocyte"])), table(cond[cell_of == "Cardiomyocyte"]), collapse = " "), "\n")
pk <- data.table(do.call(rbind, strsplit(rownames(x), "-"))); setnames(pk, c("chr", "s", "e"))
pk[, `:=`(s = as.integer(s), e = as.integer(e), peak_id = rownames(x), wkb = (as.integer(e) - as.integer(s)) / 1000)]
for (ct in c("Cardiomyocyte", "Fibroblast")) for (cd in c("YS", "MAS", "MAE")) {
  cols <- which(cell_of == ct & cond == cd & !is.na(cond))
  cpm <- sweep(x[, cols, drop = FALSE], 2, lib[cols] / 1e6, "/")
  pk[[paste0(if (ct == "Cardiomyocyte") "cm_" else "fb_", tolower(cd))]] <- as.numeric(rowMeans(cpm)) / pk$wkb
}
setkey(pk, chr, s, e)
one <- function(f, label) {
  s <- fread(file.path(data_dir, "derived", f)); q <- s[, .(site_id, chr = seqnames, s = motif_mid, e = motif_mid)]
  ov <- foverlaps(q, pk, type = "within", nomatch = NA); ov <- ov[!duplicated(site_id)]
  ov[, `:=`(site_table = label)]
  ov[, .(site_id, site_table, snatac_peak = peak_id, cm_ys, cm_mas, cm_mae, fb_ys, fb_mas, fb_mae)]
}
res <- rbind(one("cm_fos_sites.csv", "cm"), one("fb_fos_sites.csv", "fb"))
fwrite(res, file.path(data_dir, "derived", "measured_by_condition.csv"))
cat("written", nrow(res), "rows;", sum(!is.na(res$snatac_peak)), "with a peak\n")
