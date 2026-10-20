#!/usr/bin/env Rscript
# Step 1: define AP-1 (Fos) site tables from the TOBIAS table, and cross-map fibroblast MPRA elements.
#
# All Fos motif instances are kept (no peak, score or snATAC filter). The slide counts (2,663 CM sites,
# 2,320 Stable) could not be reproduced from this table, so n here is the true n of the table.
#
# Binding labels (YS / MAS / MAE = young sedentary / aged sedentary / aged + exercise), three independent views of
# the same 7 binding patterns (a site can carry all three, so ageing, exercise and their combination can be compared):
#   age_class  (YS -> MAS):  Stable (bound in both), Age Up (MAS only), Age Down (YS only), Unbound (neither)
#   ex_class   (MAS -> MAE): Stable (bound in both), Exercise Up (MAE only), Exercise Down (MAS only), Unbound (neither)
#   binding_category: the full YS/MAS/MAE pattern, one of the 7 bound combinations (111, 110, 101, 100, 011, 010, 001)
# stable_dynamic: Stable = bound in all three conditions (111), Dynamic = anything else.
#
# Fibroblast: the MPRA elements (170 bp, mm10) and the TOBIAS fibroblast sites differ in number, location and
# called binding, so the cross-map keeps both definitions side by side with an overlap_status column.
#
# Usage: Rscript 01_define_sites.R [data_dir] [out_dir]

suppressPackageStartupMessages(library(data.table))

args <- commandArgs(trailingOnly = TRUE)
data_dir <- if (length(args) >= 1) args[1] else "Data"
out_dir <- if (length(args) >= 2) args[2] else file.path(data_dir, "derived")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

TF_NAME <- "Fos"
BINDING_CATEGORY <- c(
  "111" = "YS+MAS+MAE", "110" = "YS+MAS", "101" = "YS+MAE", "100" = "YS only",
  "011" = "MAS+MAE", "010" = "MAS only", "001" = "MAE only"
)
tob <- fread(
  cmd = paste("gunzip -c", shQuote(file.path(data_dir, "Tobias_df_Classified_Full.csv.gz"))),
  select = c("seqnames", "start", "end", "peak_chr", "peak_start", "peak_end", "motif_identifier",
             "in_YS", "in_MAS", "in_MAE", "YS_score", "MAS_score", "MAE_score", "category",
             "CellType", "name")
)

define_sites <- function(cell_type) {
  s <- tob[CellType == cell_type & name == TF_NAME]
  s[, pattern := paste0(as.integer(in_YS), as.integer(in_MAS), as.integer(in_MAE))]
  s <- s[pattern != "000"]
  s[, `:=`(
    site_id = motif_identifier,
    binding_category = BINDING_CATEGORY[pattern],
    age_class = fcase(in_YS & in_MAS, "Stable", !in_YS & in_MAS, "Age Up", in_YS & !in_MAS, "Age Down", default = "Unbound"),
    ex_class = fcase(in_MAS & in_MAE, "Stable", !in_MAS & in_MAE, "Exercise Up", in_MAS & !in_MAE, "Exercise Down", default = "Unbound"),
    stable_dynamic = ifelse(pattern == "111", "Stable", "Dynamic"),
    motif_mid = (start + end) %/% 2L,
    peak_width = peak_end - peak_start,
    peak_coord = paste(peak_chr, peak_start, peak_end, sep = "-")
  )]
  setcolorder(s, c("site_id", "seqnames", "start", "end", "motif_mid", "peak_chr", "peak_start", "peak_end",
                   "peak_coord", "peak_width", "pattern", "binding_category", "age_class", "ex_class", "stable_dynamic"))
  s
}

summarise <- function(s, label) {
  cat(sprintf("\n== %s %s sites: %d (peaks: %d) ==\n", label, TF_NAME, nrow(s), uniqueN(s$peak_coord)))
  print(s[, .N, by = .(stable_dynamic)])
  print(s[, .N, by = .(pattern, binding_category)][order(-N)])
  print(s[, .N, by = age_class][order(-N)])
  print(s[, .N, by = ex_class][order(-N)])
  print(s[, .N, by = .(age_class, ex_class)][order(-N)])
}

cm <- define_sites("Cardiomyocyte")
fb <- define_sites("Fibroblast")
fwrite(cm, file.path(out_dir, "cm_fos_sites.csv"))
fwrite(fb, file.path(out_dir, "fb_fos_sites.csv"))
summarise(cm, "Cardiomyocyte")
summarise(fb, "Fibroblast")

# ---- fibroblast MPRA <-> TOBIAS cross-map ----
mpra <- fread(file.path(data_dir, "MPRA_Results_with_TSS_distance.csv"))
mpra[, `:=`(mpra_id = Coords, mpra_group = Group, mpra_category = category, mpra_replicated = replicated)]
mpra[, mpra_stable_dynamic := ifelse(Group == "Stable", "Stable", "Dynamic")]
cat("\nMPRA chromosomes:", paste(sort(unique(mpra$chr)), collapse = " "), "\n")

m <- mpra[, .(mpra_id, chr, s = start, e = end, mpra_group, mpra_category, mpra_stable_dynamic,
              mpra_replicated, log2_alpha_Endo, log2_alpha_MotifShuff, log2_alpha_FullShuff,
              delta_vs_MotifShuff, delta_vs_FullShuff, GC_content, dist_to_TSS)]
setkey(m, chr, s, e)
t <- fb[, .(site_id, chr = seqnames, s = start, e = end, tobias_pattern = pattern,
            tobias_category = category, tobias_binding_category = binding_category,
            tobias_age_class = age_class, tobias_ex_class = ex_class, tobias_stable_dynamic = stable_dynamic)]

# motif overlapping an MPRA element (any overlap; overlap_bp recorded)
hit <- foverlaps(t, m, type = "any", nomatch = NA, by.x = c("chr", "s", "e"), by.y = c("chr", "s", "e"))
setnames(hit, c("s", "e", "i.s", "i.e"), c("mpra_start", "mpra_end", "motif_start", "motif_end"))
hit[, overlap_bp := pmin(mpra_end, motif_end) - pmax(mpra_start, motif_start) + 1L]
hit[, overlap_status := fifelse(is.na(mpra_id), "Tobias only", "Matched")]

mpra_only <- mpra[!mpra_id %in% hit$mpra_id,
                  .(mpra_id, chr, mpra_start = start, mpra_end = end, mpra_group, mpra_category,
                    mpra_stable_dynamic, mpra_replicated, log2_alpha_Endo, log2_alpha_MotifShuff,
                    log2_alpha_FullShuff, delta_vs_MotifShuff, delta_vs_FullShuff, GC_content, dist_to_TSS)]
mpra_only[, overlap_status := "MPRA only"]

xmap <- rbindlist(list(hit, mpra_only), use.names = TRUE, fill = TRUE)
fwrite(xmap, file.path(out_dir, "fb_mpra_tobias_crossmap.csv"))

cat(sprintf("\n== Fibroblast MPRA (%d elements) vs TOBIAS Fos (%d sites) ==\n", nrow(mpra), nrow(fb)))
cat("MPRA elements with >=1 Tobias Fos motif:", uniqueN(hit[!is.na(mpra_id), mpra_id]), "\n")
cat("Tobias Fos sites in an MPRA element:", uniqueN(hit[!is.na(mpra_id), site_id]), "\n")
print(xmap[, .N, by = overlap_status])
print(xmap[overlap_status == "Matched", .N, by = .(mpra_group, tobias_binding_category)][order(-N)])
