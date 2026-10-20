#!/usr/bin/env Rscript
# Step 3: measured snATAC accessibility at each Fos site, from the pseudobulk count matrix.
#
# For every site, the snATAC peak containing the motif midpoint is found (peaks come from the rds row names);
# accessibility = mean over the samples of the matching cell type of counts / library size (CPM) / peak width (kb).
# Sites without a containing peak get NA. Both CM and fibroblast columns are given for every site, so each site
# table can be validated against its own cell type and the other one.
#
# Usage: Rscript 03_measured_atac.R [data_dir] [out_csv]

suppressPackageStartupMessages({library(data.table); library(Matrix)})

args <- commandArgs(trailingOnly = TRUE)
data_dir <- if (length(args) >= 1) args[1] else "Data"
out_csv <- if (length(args) >= 2) args[2] else file.path(data_dir, "derived", "measured_atac_sites.csv")

x <- readRDS(file.path(data_dir, "heart_snatac_pseudobulk_raw_counts.rds"))
lib <- colSums(x)
cell_of <- sub("--.*", "", colnames(x))
cpm_kb <- function(cell_type) {
  cols <- which(cell_of == cell_type)
  cpm <- sweep(x[, cols, drop = FALSE], 2, lib[cols] / 1e6, "/")
  as.numeric(rowMeans(cpm))
}

pk <- data.table(do.call(rbind, strsplit(rownames(x), "-")))
setnames(pk, c("chr", "s", "e"))
pk[, `:=`(s = as.integer(s), e = as.integer(e), peak_id = rownames(x), width_kb = (as.integer(e) - as.integer(s)) / 1000)]
pk[, `:=`(cm_cpm_per_kb = cpm_kb("Cardiomyocyte") / width_kb, fb_cpm_per_kb = cpm_kb("Fibroblast") / width_kb)]
setkey(pk, chr, s, e)

one <- function(sites_csv, label) {
  s <- fread(file.path(data_dir, "derived", sites_csv))
  q <- s[, .(site_id, chr = seqnames, s = motif_mid, e = motif_mid)]
  ov <- foverlaps(q, pk, type = "within", nomatch = NA)
  ov <- ov[!duplicated(site_id)]
  out <- ov[, .(site_id, snatac_peak = peak_id, snatac_peak_width = e - s,
                cm_cpm_per_kb, fb_cpm_per_kb)]
  out[, site_table := label]
  cat(sprintf("%s: %d sites, %d with an snATAC peak\n", label, nrow(out), sum(!is.na(out$snatac_peak))))
  out
}
res <- rbind(one("cm_fos_sites.csv", "cm"), one("fb_fos_sites.csv", "fb"))
# peak width was overwritten by the join columns; recompute explicitly
res[, snatac_peak_width := as.integer(sub(".*-([0-9]+)$", "\\1", snatac_peak)) - as.integer(sub("^[^-]+-([0-9]+)-.*$", "\\1", snatac_peak))]
res[, `:=`(log2_cm_cpm_per_kb = log2(cm_cpm_per_kb + 1), log2_fb_cpm_per_kb = log2(fb_cpm_per_kb + 1))]
fwrite(res, out_csv)
print(summary(res[, .(cm_cpm_per_kb, fb_cpm_per_kb)]))
