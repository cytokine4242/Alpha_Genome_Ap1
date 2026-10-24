#!/usr/bin/env Rscript
# Publication-style versions of the four fibroblast slide figures: Arial, distributions rather than bare bars,
# bootstrap confidence intervals on every summary, and the test printed on the panel it belongs to.
#
#   Rscript Scripts/slides/fb_slide_figures.R
#
# Inputs are the analysis tables, so the figures rebuild whenever the analyses are rerun:
#   Analysis/slide_figures/s_fb_transplant_per_donor.csv          (make_fb_simple_figures.py)
#   Analysis/fb_swap/tables/sequence_vs_location_ci.csv           (55_fb_swap_regression_ci.py)
#   Analysis/fb_partner_shuffle/tables/own_vs_partners_per_site.csv (52_fb_partner_shuffle_analysis.py)
# Outputs Analysis/slide_figures/r_fb_*.png

suppressPackageStartupMessages({
  library(ggplot2); library(ragg); library(grid)
})

ROOT <- "/Users/z5162865/WongLab/AlphaGenome_Ap1"
OUT  <- file.path(ROOT, "Analysis", "slide_figures")
set.seed(55)

INK <- "#111111"; MUTED <- "#4d4d4d"; RULE <- "#111111"
STABLE <- "#4F4F4F"; ORANGE <- "#C2551C"; GREY <- "#AFAEAA"; BLUE <- "#1F6FB2"

# --- house style: Arial, hairline axes, ticks out, nothing on the panel that is not data ------------------
theme_pub <- function(base = 15) {
  theme_classic(base_size = base, base_family = "Arial") +
    theme(
      axis.line        = element_line(colour = RULE, linewidth = 0.4, lineend = "square"),
      axis.ticks       = element_line(colour = RULE, linewidth = 0.4),
      axis.ticks.length = unit(3.2, "pt"),
      axis.text        = element_text(colour = INK, size = base),
      axis.title       = element_text(colour = INK, size = base),
      axis.title.x     = element_text(margin = margin(t = 7)),
      axis.title.y     = element_text(margin = margin(r = 7)),
      legend.position  = "top",
      legend.justification = "left",
      legend.title     = element_blank(),
      legend.key       = element_blank(),
      legend.key.size  = unit(11, "pt"),
      legend.text      = element_text(colour = INK, size = base - 1.5),
      legend.margin    = margin(b = -2),
      plot.margin      = margin(6, 10, 4, 6),
      plot.background  = element_rect(fill = "white", colour = NA),
      panel.background = element_rect(fill = "white", colour = NA)
    )
}

# median with a percentile bootstrap interval - the summary every panel is built on
boot_med <- function(x, n = 4000) {
  x <- x[is.finite(x)]
  b <- replicate(n, median(sample(x, length(x), replace = TRUE)))
  data.frame(y = median(x), lo = unname(quantile(b, .025)), hi = unname(quantile(b, .975)), n = length(x))
}

# P formatted the way a journal would set it
fmtP <- function(p, floor = NA) {
  if (!is.na(floor) && p < floor) return(sprintf("P < %s", format(floor, scientific = FALSE)))
  if (p == 0) return("P < 0.001")
  if (p >= 0.05) return(sprintf("P = %.2f", p))
  if (p >= 1e-4) return(sprintf("P = %s", format(signif(p, 2), scientific = FALSE)))
  e <- floor(log10(p)); m <- signif(p / 10^e, 2)
  sup <- c("⁰","¹","²","³","⁴","⁵","⁶","⁷","⁸","⁹")
  dig <- paste(sup[as.integer(strsplit(as.character(abs(e)), "")[[1]]) + 1], collapse = "")
  sprintf("P = %g × 10⁻%s", m, dig)
}

# a significance bracket: two uprights and a rule, with the test above it
bracket <- function(x1, x2, y, label, tick = NULL, size = 4.6, colour = INK, vjust = -0.35) {
  tick <- if (is.null(tick)) abs(y) * 0.03 + 0.02 else tick
  list(
    annotate("segment", x = x1, xend = x1, y = y - tick, yend = y, linewidth = 0.4, colour = colour),
    annotate("segment", x = x2, xend = x2, y = y - tick, yend = y, linewidth = 0.4, colour = colour),
    annotate("segment", x = x1, xend = x2, y = y, yend = y, linewidth = 0.4, colour = colour),
    annotate("text", x = (x1 + x2) / 2, y = y, label = label, vjust = vjust,
             family = "Arial", size = size, colour = colour)
  )
}

save_fig <- function(p, name, w, h) {
  f <- file.path(OUT, paste0(name, ".png"))
  agg_png(f, width = w, height = h, units = "in", res = 400, background = "white", scaling = 1)
  print(p); invisible(dev.off())
  cat("  wrote", f, "\n")
}

# ==========================================================================================================
# 1. elements transplanted into closed chromatin
# ==========================================================================================================
td <- read.csv(file.path(OUT, "s_fb_transplant_per_donor.csv"))
td$grp <- factor(td$grp, levels = c("Stable", "Dynamic"))
ts <- do.call(rbind, lapply(levels(td$grp), function(g) cbind(grp = g, boot_med(td$lift[td$grp == g]))))
ts$grp <- factor(ts$grp, levels = levels(td$grp))
p_tr <- wilcox.test(lift ~ grp, data = td)$p.value
labs_tr <- sprintf("%s\n(n = %d)", ts$grp, ts$n)

g1 <- ggplot(td, aes(grp, lift, colour = grp, fill = grp)) +
  geom_hline(yintercept = 0, linewidth = 0.3, colour = "#9a9993") +
  geom_point(position = position_jitter(width = 0.20, height = 0, seed = 1),
             size = 1.15, alpha = 0.28, stroke = 0, show.legend = FALSE) +
  geom_errorbar(data = ts, aes(grp, y = y, ymin = lo, ymax = hi), width = 0, linewidth = 1.5,
                inherit.aes = FALSE, colour = INK) +
  geom_point(data = ts, aes(grp, y), inherit.aes = FALSE, size = 3.4, colour = INK) +
  geom_text(data = ts, aes(grp, y, label = sprintf("%.1f×", 2^y)), inherit.aes = FALSE,
            hjust = -0.42, vjust = 0.35, family = "Arial", fontface = "bold", size = 6.4, colour = INK) +
  bracket(1, 2, 8.0, fmtP(p_tr), tick = 0.28) +
  scale_colour_manual(values = c(Stable = STABLE, Dynamic = ORANGE)) +
  scale_fill_manual(values = c(Stable = STABLE, Dynamic = ORANGE)) +
  scale_y_continuous(breaks = c(0, 2, 4, 6, 8), labels = c("1", "4", "16", "64", "256")) +
  scale_x_discrete(expand = expansion(add = 0.62), labels = labs_tr) +
  coord_cartesian(ylim = c(-0.35, 8.9), expand = FALSE, clip = "on") +
  labs(x = NULL, y = "Opening of a closed site (fold)") +
  theme_pub() + theme(legend.position = "none")
save_fig(g1, "r_fb_transplant", 6.6, 4.3)

# ==========================================================================================================
# 2. what sets a swapped site's level: the element, or the place?
# ==========================================================================================================
cif <- file.path(ROOT, "Analysis/fb_swap/tables/sequence_vs_location_ci.csv")
if (file.exists(cif)) {
  ci <- read.csv(cif)
  sl <- rbind(
    data.frame(width = ci$width, term = "What the element is",
               beta = ci$beta_donor_sequence, lo = ci$seq_lo, hi = ci$seq_hi, p = ci$p_sequence_vs_location),
    data.frame(width = ci$width, term = "Where it sits",
               beta = ci$beta_recipient_location, lo = ci$loc_lo, hi = ci$loc_hi, p = ci$p_sequence_vs_location))
  sl$width <- factor(sl$width, levels = c("170", "peak"),
                     labels = c("170 bp element\n(the MPRA window)", "Whole peak\n(median 0.8 kb)"))
  sl$term <- factor(sl$term, levels = c("What the element is", "Where it sits"))
  ptab <- unique(sl[, c("width", "p")])
  ytop <- max(sl$hi) + 0.10

  g2 <- ggplot(sl, aes(width, beta, colour = term)) +
    geom_hline(yintercept = 0, linewidth = 0.3, colour = "#9a9993") +
    geom_errorbar(aes(ymin = lo, ymax = hi), width = 0, linewidth = 1.4,
                  position = position_dodge(width = 0.44)) +
    geom_point(size = 4.2, position = position_dodge(width = 0.44)) +
    geom_text(aes(label = sprintf("%.2f", beta)), position = position_dodge(width = 0.44),
              hjust = -0.45, vjust = 0.38, family = "Arial", size = 5.0, show.legend = FALSE) +
    bracket(0.78, 1.22, ytop, fmtP(ptab$p[1], floor = 5e-4), tick = 0.035, size = 4.4) +
    bracket(1.78, 2.22, ytop, fmtP(ptab$p[2], floor = 5e-4), tick = 0.035, size = 4.4) +
    scale_colour_manual(values = c("What the element is" = BLUE, "Where it sits" = GREY)) +
    scale_y_continuous(breaks = seq(0, 0.8, 0.2)) +
    scale_x_discrete(expand = expansion(add = 0.62)) +
    coord_cartesian(ylim = c(0, ytop + 0.14), expand = FALSE, clip = "on") +
    labs(x = NULL, y = "Weight on the swapped\nsite's accessibility") +
    theme_pub()
  save_fig(g2, "r_fb_scale", 7.0, 4.3)
} else {
  cat("  skipped r_fb_scale: run Scripts/55_fb_swap_regression_ci.py first\n")
}

# ==========================================================================================================
# 3. the site's own AP-1 motif against its partners, by class
# ==========================================================================================================
J <- read.csv(file.path(ROOT, "Analysis/fb_partner_shuffle/tables/own_vs_partners_per_site.csv"))
keep <- c("Stable", "AP-1 gain", "AP-1 loss")
CL <- c(Stable = "Stable", `AP-1 gain` = "AP-1 gained\nwith age", `AP-1 loss` = "AP-1 lost\nwith age")
L <- rbind(
  data.frame(grp = J$grp, what = "Its own AP-1 motif", eff = J$own_ap1),
  data.frame(grp = J$grp, what = "Its partner motifs", eff = J$partners_mean_net))
L <- L[L$grp %in% keep & is.finite(L$eff), ]
L$grp <- factor(L$grp, levels = keep, labels = CL[keep])
L$what <- factor(L$what, levels = c("Its own AP-1 motif", "Its partner motifs"))

S3 <- do.call(rbind, lapply(split(L, list(L$grp, L$what), drop = TRUE), function(d)
  cbind(grp = levels(L$grp)[d$grp[1]], what = levels(L$what)[d$what[1]], boot_med(d$eff))))
S3$grp <- factor(S3$grp, levels = levels(L$grp)); S3$what <- factor(S3$what, levels = levels(L$what))

pv <- function(w) {
  a <- L$eff[L$grp == CL["Stable"] & L$what == w]
  b <- L$eff[L$grp == CL["AP-1 loss"] & L$what == w]
  wilcox.test(a, b)$p.value
}
p_own <- pv("Its own AP-1 motif"); p_par <- pv("Its partner motifs")
dw <- 0.42

g3 <- ggplot(L, aes(grp, eff, fill = what, colour = what)) +
  geom_hline(yintercept = 0, linewidth = 0.3, colour = "#9a9993") +
  geom_violin(position = position_dodge(width = dw), width = 0.36, linewidth = 0,
              alpha = 0.24, trim = TRUE, show.legend = FALSE) +
  geom_errorbar(data = S3, aes(grp, y = y, ymin = lo, ymax = hi, group = what), inherit.aes = FALSE,
                position = position_dodge(width = dw), width = 0, linewidth = 1.4, colour = INK) +
  geom_point(data = S3, aes(grp, y, group = what, fill = what), inherit.aes = FALSE,
             position = position_dodge(width = dw), size = 3.4, shape = 21, colour = INK, stroke = 0.6) +
  geom_text(data = S3, aes(grp, y, group = what, label = sprintf("%.2f", y)), inherit.aes = FALSE,
            position = position_dodge(width = dw), vjust = 2.1, family = "Arial", size = 4.6, colour = INK) +
  bracket(1 - dw / 2, 3 - dw / 2, 0.62, fmtP(p_own), tick = 0.10, size = 4.5, colour = ORANGE) +
  bracket(1 + dw / 2, 3 + dw / 2, 0.24, paste0(fmtP(p_par), " (n.s.)"), tick = 0.10, size = 4.5, colour = MUTED) +
  scale_fill_manual(values = c("Its own AP-1 motif" = ORANGE, "Its partner motifs" = GREY)) +
  scale_colour_manual(values = c("Its own AP-1 motif" = ORANGE, "Its partner motifs" = GREY)) +
  scale_x_discrete(expand = expansion(add = 0.45)) +
  scale_y_continuous(breaks = seq(-2.5, 0, 0.5)) +
  coord_cartesian(ylim = c(-2.75, 1.02), expand = FALSE, clip = "on") +
  labs(x = NULL, y = expression(paste("Effect of shuffling it (log"[2], ")"))) +
  theme_pub()
save_fig(g3, "r_fb_specificity", 8.4, 4.5)

# ==========================================================================================================
# 4. how much of a site's motif dependence is AP-1?
# ==========================================================================================================
J$share <- abs(J$own_ap1) / (abs(J$own_ap1) + abs(J$partners_mean_net))
J$sd2 <- factor(ifelse(J$grp == "Stable", "Stable", "Dynamic"), levels = c("Stable", "Dynamic"))
K <- J[is.finite(J$share), ]
S4 <- do.call(rbind, lapply(levels(K$sd2), function(g) cbind(sd2 = g, boot_med(K$share[K$sd2 == g]))))
S4$sd2 <- factor(S4$sd2, levels = levels(K$sd2))
p_sh <- wilcox.test(share ~ sd2, data = K)$p.value

g4 <- ggplot(K, aes(sd2, share, fill = sd2, colour = sd2)) +
  geom_violin(width = 0.72, linewidth = 0, alpha = 0.26, trim = TRUE, show.legend = FALSE) +
  geom_boxplot(width = 0.14, outlier.shape = NA, linewidth = 0.45, fill = "white",
               colour = INK, show.legend = FALSE) +
  geom_errorbar(data = S4, aes(sd2, y = y, ymin = lo, ymax = hi), inherit.aes = FALSE,
                width = 0, linewidth = 1.6, colour = INK) +
  geom_point(data = S4, aes(sd2, y, fill = sd2), inherit.aes = FALSE, size = 3.4, shape = 21,
             colour = INK, stroke = 0.6, show.legend = FALSE) +
  geom_text(data = S4, aes(sd2, y, label = sprintf("%.0f%%", 100 * y)), inherit.aes = FALSE,
            hjust = -0.95, family = "Arial", fontface = "bold", size = 6.2, colour = INK) +
  bracket(1, 2, 1.055, fmtP(p_sh), tick = 0.028, size = 4.8) +
  scale_fill_manual(values = c(Stable = STABLE, Dynamic = ORANGE)) +
  scale_colour_manual(values = c(Stable = STABLE, Dynamic = ORANGE)) +
  scale_x_discrete(expand = expansion(add = 0.62), labels = sprintf("%s\n(n = %d)", S4$sd2, S4$n)) +
  scale_y_continuous(breaks = seq(0, 1, 0.25), labels = scales::percent_format(accuracy = 1)) +
  coord_cartesian(ylim = c(-0.02, 1.19), expand = FALSE, clip = "on") +
  labs(x = NULL, y = "Share of the site's motif\ndependence carried by AP-1") +
  theme_pub() + theme(legend.position = "none")
save_fig(g4, "r_fb_share", 6.6, 4.5)

cat("done\n")
