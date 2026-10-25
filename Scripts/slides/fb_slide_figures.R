#!/usr/bin/env Rscript
# The four fibroblast slide figures. Arial, the distribution behind every summary, and the test printed on the panel.
#
#   Rscript Scripts/slides/fb_slide_figures.R
#
# Each panel uses a different form on purpose - four variations of the same dot-and-whisker read as one slide shown
# four times. Horizontal boxes, then diverging bars, then notched boxes, then a stacked composition.
#
# Classes are coloured the way Emily's deck colours them, so the two talks agree on sight:
# blue = Stable, crimson = age-dynamic, pale gold = the AP-1 core.
#
# Inputs are the analysis tables, so the figures rebuild whenever the analyses are rerun:
#   Analysis/slide_figures/s_fb_transplant_per_donor.csv            (make_fb_simple_figures.py)
#   Analysis/fb_swap/tables/sequence_vs_location_ci.csv             (55_fb_swap_regression_ci.py)
#   Analysis/fb_partner_shuffle/tables/own_vs_partners_per_site.csv (52_fb_partner_shuffle_analysis.py)
#   Analysis/fb_scan/tables/top_tile_zone_by_class.csv              (56_fb_top_tile_rank.py)
# Outputs Analysis/slide_figures/r_fb_*.png

suppressPackageStartupMessages({
  library(ggplot2); library(ragg); library(grid); library(patchwork)
})

ROOT <- "/Users/z5162865/WongLab/AlphaGenome_Ap1"
OUT  <- file.path(ROOT, "Analysis", "slide_figures")
set.seed(55)

INK <- "#111111"; MUTED <- "#4d4d4d"; RULE <- "#111111"; FAINT <- "#9a9993"
BLUE <- "#347CBA"            # Stable
ROSE <- "#C73A6B"            # Dynamic, pooled and the age classes
ROSE_DK <- "#9E1F4C"; ROSE_LT <- "#E291AC"
GOLD <- "#E3B23C"; GOLD_PALE <- "#F9EAC4"   # the AP-1 core
TEAL <- "#2F6F6F"            # "what the element is", used nowhere else
GREY <- "#BDBCB8"; GREY_DK <- "#8C8B87"

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
# categorical rows: the category axis carries no scale, so it gets no line and no ticks
no_cat_axis <- function(which = "y") {
  if (which == "y") theme(axis.line.y = element_blank(), axis.ticks.y = element_blank())
  else theme(axis.line.x = element_blank(), axis.ticks.x = element_blank())
}

# P formatted the way a journal would set it
fmtP <- function(p, floor = NA) {
  if (!is.na(floor) && p < floor) return(sprintf("P < %s", format(floor, scientific = FALSE)))
  if (p == 0) return("P < 0.001")
  if (p >= 1e-4) return(sprintf("P = %s", format(signif(p, 2), scientific = FALSE)))
  e <- floor(log10(p)); m <- signif(p / 10^e, 2)
  sup <- c("⁰","¹","²","³","⁴","⁵","⁶","⁷","⁸","⁹")
  dig <- paste(sup[as.integer(strsplit(as.character(abs(e)), "")[[1]]) + 1], collapse = "")
  sprintf("P = %g × 10⁻%s", m, dig)
}

# a significance bracket: two uprights and a rule, with the test above it
bracket <- function(x1, x2, y, label, tick, size = 4.5, colour = INK, vjust = -0.35) {
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
# 1. elements transplanted into closed chromatin - horizontal boxes, one row per class, donors behind them
# ==========================================================================================================
td <- read.csv(file.path(OUT, "s_fb_transplant_per_donor.csv"))
CLS <- c("Stable", "Age Up", "Age Down", "Exercise Up", "Exercise Down")
td <- td[td$donor_class %in% CLS, ]
td$cls <- factor(td$donor_class, levels = rev(CLS))          # rev: ggplot draws the first level at the bottom
p_tr <- wilcox.test(lift ~ grp, data = td)$p.value
med <- aggregate(lift ~ cls, td, median)
nn  <- aggregate(lift ~ cls, td, length); names(nn)[2] <- "n"
labs_tr <- setNames(sprintf("%s", levels(td$cls)), levels(td$cls))
COL1 <- c(Stable = BLUE, `Age Up` = ROSE, `Age Down` = ROSE, `Exercise Up` = ROSE, `Exercise Down` = ROSE)

VCOL <- 8.95                                                  # the median column, parked outside the data range

g1 <- ggplot(td, aes(cls, lift)) +
  geom_hline(yintercept = 0, linewidth = 0.35, colour = FAINT) +
  geom_point(aes(colour = donor_class), position = position_jitter(width = 0.17, height = 0, seed = 1),
             size = 1.2, alpha = 0.30, stroke = 0, show.legend = FALSE) +
  geom_boxplot(aes(colour = donor_class, fill = donor_class), width = 0.52, alpha = 0.16,
               outlier.shape = NA, linewidth = 0.6, show.legend = FALSE) +
  geom_text(data = med, aes(cls, VCOL, label = sprintf("%.0f×", 2^lift), colour = cls), hjust = 1,
            family = "Arial", fontface = "bold", size = 5.6, show.legend = FALSE) +
  annotate("text", x = 5.72, y = VCOL, label = "median", hjust = 1, family = "Arial", size = 4.2, colour = MUTED) +
  annotate("text", x = 5.72, y = 0, label = sprintf("Stable vs Dynamic (pooled), %s", fmtP(p_tr)),
           hjust = 0, family = "Arial", size = 4.4, colour = INK) +
  scale_colour_manual(values = COL1) + scale_fill_manual(values = COL1) +
  scale_y_continuous(breaks = c(0, 2, 4, 6, 8), labels = c("1", "4", "16", "64", "256")) +
  scale_x_discrete(labels = labs_tr, expand = expansion(add = c(0.6, 1.25))) +
  coord_flip(ylim = c(-0.55, 8.35), expand = FALSE, clip = "off") +
  labs(x = NULL, y = "Opening of a closed site (fold)") +
  theme_pub() + no_cat_axis("y") +
  theme(axis.text.y = element_text(colour = c(rep(ROSE, 4), BLUE), hjust = 0, size = 15),
        plot.margin = margin(6, 56, 4, 6))
save_fig(g1, "r_fb_transplant", 10.6, 3.85)
cat(sprintf("     per-class median fold: %s\n",
            paste(sprintf("%s %.1fx", med$cls, 2^med$lift), collapse = ", ")))

# ==========================================================================================================
# 2. the element or the place? - the two weights pulling against each other from a common zero
# ==========================================================================================================
cif <- file.path(ROOT, "Analysis/fb_swap/tables/sequence_vs_location_ci.csv")
if (file.exists(cif)) {
  ci <- read.csv(cif)
  # one bar each side of zero: the element pulls left, the location right. Values are the same betas, signed for side.
  sl <- rbind(
    data.frame(width = ci$width, term = "What the element is", s = -1,
               beta = ci$beta_donor_sequence, lo = ci$seq_lo, hi = ci$seq_hi, p = ci$p_sequence_vs_location),
    data.frame(width = ci$width, term = "Where it sits", s = 1,
               beta = ci$beta_recipient_location, lo = ci$loc_lo, hi = ci$loc_hi, p = ci$p_sequence_vs_location))
  sl$x <- sl$s * sl$beta; sl$xlo <- sl$s * sl$lo; sl$xhi <- sl$s * sl$hi
  sl$width <- factor(sl$width, levels = c("peak", "170"),
                     labels = c("Whole peak\nmedian 0.8 kb", "170 bp element\nthe MPRA window"))
  sl$term <- factor(sl$term, levels = c("What the element is", "Where it sits"))
  ptab <- unique(sl[, c("width", "p")])
  LIM <- 1.02

  g2 <- ggplot(sl, aes(width, x, fill = term)) +
    geom_col(width = 0.44, show.legend = FALSE) +
    geom_errorbar(aes(ymin = xlo, ymax = xhi, colour = term), width = 0.11, linewidth = 0.55,
                  show.legend = FALSE) +
    geom_text(aes(y = x, label = sprintf("%.2f", beta), hjust = ifelse(s < 0, 1.25, -0.25)),
              family = "Arial", fontface = "bold", size = 5.4, colour = INK, show.legend = FALSE) +
    geom_hline(yintercept = 0, linewidth = 0.8, colour = INK) +
    annotate("text", x = 2.82, y = -LIM, label = "What the element IS", hjust = 0, vjust = 0.5,
             family = "Arial", fontface = "bold", size = 5.0, colour = TEAL) +
    annotate("text", x = 2.82, y = LIM, label = "WHERE it sits", hjust = 1, vjust = 0.5,
             family = "Arial", fontface = "bold", size = 5.0, colour = "#6E6D69") +
    annotate("text", x = c(2.47, 1.47), y = 0, label = sapply(ptab$p, fmtP, floor = 5e-4),
             hjust = 0.5, vjust = 0.5, family = "Arial", size = 4.3, colour = MUTED) +
    scale_fill_manual(values = c("What the element is" = TEAL, "Where it sits" = GREY)) +
    scale_colour_manual(values = c("What the element is" = "white", "Where it sits" = "#6E6D69")) +
    scale_y_continuous(breaks = seq(-0.8, 0.8, 0.4), labels = function(b) sprintf("%.1f", abs(b))) +
    scale_x_discrete(expand = expansion(add = c(0.78, 1.05))) +
    coord_flip(ylim = c(-LIM, LIM), expand = FALSE, clip = "off") +
    labs(x = NULL, y = "Weight on the swapped site's accessibility") +
    theme_pub() + no_cat_axis("y") +
    theme(axis.text.y = element_text(colour = INK, hjust = 0.5, size = 14.5, lineheight = 0.95),
          plot.margin = margin(6, 14, 4, 6))
  save_fig(g2, "r_fb_scale", 10.4, 4.0)
} else {
  cat("  skipped r_fb_scale: run Scripts/55_fb_swap_regression_ci.py first\n")
}

# ==========================================================================================================
# 3. the site's own AP-1 motif against its partners, Dynamic split into the age classes as given
#    notched boxes: the notch is a 95% interval on the median, so the comparison is inside the box
# ==========================================================================================================
J <- read.csv(file.path(ROOT, "Analysis/fb_partner_shuffle/tables/own_vs_partners_per_site.csv"))
ORD <- c("Stable", "AP-1 gain", "other Dynamic", "AP-1 loss")        # ordered by how far they sit from Stable
CL <- c(Stable = "Stable", `AP-1 gain` = "AP-1 gained\nwith age",
        `other Dynamic` = "Dynamic, not\nage-driven", `AP-1 loss` = "AP-1 lost\nwith age")
COL3 <- setNames(c(BLUE, ROSE_LT, ROSE, ROSE_DK), CL[ORD])
L <- rbind(
  data.frame(grp = J$grp, what = "Its own AP-1 motif", eff = J$own_ap1),
  data.frame(grp = J$grp, what = "Its partner motifs", eff = J$partners_mean_net))
L <- L[L$grp %in% ORD & is.finite(L$eff), ]
L$grp <- factor(L$grp, levels = ORD, labels = CL[ORD])
L$what <- factor(L$what, levels = c("Its own AP-1 motif", "Its partner motifs"))

# box drawn from precomputed statistics so the whiskers stop at the 5th and 95th centile and nothing runs off the
# panel; the notch is the usual +-1.58 IQR / sqrt(n), a 95% interval on the median
bstat <- function(x) {
  x <- x[is.finite(x)]
  q <- unname(quantile(x, c(.05, .25, .5, .75, .95)))
  h <- 1.58 * (q[4] - q[2]) / sqrt(length(x))
  data.frame(ymin = q[1], lower = q[2], middle = q[3], upper = q[4], ymax = q[5],
             notchlower = q[3] - h, notchupper = q[3] + h, n = length(x))
}
B <- do.call(rbind, lapply(split(L, list(L$grp, L$what), drop = TRUE), function(d)
  cbind(grp = levels(L$grp)[d$grp[1]], what = levels(L$what)[d$what[1]], bstat(d$eff))))
B$grp <- factor(B$grp, levels = levels(L$grp)); B$what <- factor(B$what, levels = levels(L$what))
B$key <- paste(B$grp, B$what)
FILL3 <- c(setNames(paste0(COL3, "33"), paste(levels(L$grp), levels(L$what)[1])),
           setNames(rep("white", 4), paste(levels(L$grp), levels(L$what)[2])))

pv <- function(w, cls) {
  a <- L$eff[L$grp == CL["Stable"] & L$what == w]
  b <- L$eff[L$grp == CL[cls] & L$what == w]
  wilcox.test(a, b)$p.value
}
p_own <- sapply(ORD[-1], function(k) pv("Its own AP-1 motif", k))
p_par <- sapply(ORD[-1], function(k) pv("Its partner motifs", k))
nlab <- sprintf("%s\n(n = %d)", levels(L$grp), as.integer(table(L$grp[L$what == "Its own AP-1 motif"])))
dw <- 0.60
# the test sits over the box it belongs to, which four classes make far easier to follow than four brackets
PS <- data.frame(x = (1:4) - dw / 4, y = 0.30,
                 lab = c("reference", sapply(p_own, fmtP)))

g3 <- ggplot(B, aes(grp)) +
  geom_hline(yintercept = 0, linewidth = 0.35, colour = FAINT) +
  # with stat = "identity" the notch bounds have to be supplied; ggplot warns that it does not know the two
  # names but still carries them through to the geom, which is where they are read
  geom_boxplot(aes(ymin = ymin, lower = lower, middle = middle, upper = upper, ymax = ymax,
                   notchlower = notchlower, notchupper = notchupper, colour = grp, fill = key),
               stat = "identity", width = 0.46, position = position_dodge(width = dw),
               notch = TRUE, notchwidth = 0.58, linewidth = 0.6, show.legend = FALSE) +
  geom_text(data = PS, aes(x, y, label = lab), inherit.aes = FALSE, family = "Arial", size = 4.2,
            colour = MUTED) +
  annotate("text", x = 0.74, y = 0.56, label = "vs Stable, at the site's own AP-1 motif", hjust = 0,
           family = "Arial", size = 4.2, colour = MUTED) +
  annotate("point", x = 2.40, y = 0.95, size = 4.4, shape = 22, fill = "#1111112b", colour = INK, stroke = 0.6) +
  annotate("text", x = 2.48, y = 0.95, label = "its own AP-1 motif", hjust = 0, vjust = 0.42,
           family = "Arial", size = 4.4, colour = INK) +
  annotate("point", x = 2.40, y = 0.70, size = 4.4, shape = 22, fill = "white", colour = INK, stroke = 0.6) +
  annotate("text", x = 2.48, y = 0.70, label = "its partner motifs, net of a matched control",
           hjust = 0, vjust = 0.42, family = "Arial", size = 4.4, colour = INK) +
  scale_colour_manual(values = COL3) +
  scale_fill_manual(values = FILL3) +
  scale_x_discrete(labels = nlab, expand = expansion(add = 0.50)) +
  scale_y_continuous(breaks = seq(-2.5, 0, 0.5)) +
  coord_cartesian(ylim = c(min(B$ymin) - 0.12, 1.14), expand = FALSE, clip = "on") +
  labs(x = NULL, y = expression(paste("Effect of shuffling it (log"[2], ")"))) +
  theme_pub() + theme(legend.position = "none",
                      axis.text.x = element_text(colour = unname(COL3), size = 14, lineheight = 0.95))
save_fig(g3, "r_fb_specificity", 10.4, 4.0)
cat(sprintf("     own AP-1 vs Stable: %s\n", paste(sprintf("%s %.3g", ORD[-1], p_own), collapse = ", ")))
cat(sprintf("     partners vs Stable: %s\n", paste(sprintf("%s %.3g", ORD[-1], p_par), collapse = ", ")))

# ==========================================================================================================
# 4. how a site divides its motif dependence - a composition, split by class
#    The two parts are complementary, so a median split sums to 100% and the bar is an honest whole.
# ==========================================================================================================
J$share <- abs(J$own_ap1) / (abs(J$own_ap1) + abs(J$partners_mean_net))
K <- J[is.finite(J$share) & J$grp %in% ORD, ]
K$cls <- factor(K$grp, levels = rev(ORD), labels = rev(gsub("\n", " ", CL[ORD])))
K$sd2 <- factor(ifelse(K$grp == "Stable", "Stable", "Dynamic"), levels = c("Stable", "Dynamic"))
p_sh <- wilcox.test(share ~ sd2, data = K)$p.value
boot_med_ci <- function(x, n = 4000) {
  b <- replicate(n, median(sample(x, length(x), replace = TRUE)))
  unname(quantile(b, c(.025, .975)))
}
C4 <- do.call(rbind, lapply(levels(K$cls), function(g) {
  x <- K$share[K$cls == g]; ci <- boot_med_ci(x)
  data.frame(cls = g, n = length(x), ap1 = 100 * median(x), lo = 100 * ci[1], hi = 100 * ci[2])
}))
C4$cls <- factor(C4$cls, levels = levels(K$cls))
C4$partners <- 100 - C4$ap1
ST <- rbind(data.frame(cls = C4$cls, part = "AP-1 motif itself", pct = C4$ap1),
            data.frame(cls = C4$cls, part = "Partner motifs", pct = C4$partners))
ST$part <- factor(ST$part, levels = c("AP-1 motif itself", "Partner motifs"))

g4 <- ggplot(ST, aes(cls, pct, fill = part)) +
  geom_col(width = 0.56, colour = "white", linewidth = 1.0, position = position_stack(reverse = TRUE)) +
  geom_errorbar(data = C4, aes(cls, ymin = lo, ymax = hi), inherit.aes = FALSE,
                width = 0.17, linewidth = 0.7, colour = INK) +
  geom_text(data = C4, aes(cls, ap1 - 4.5, label = sprintf("%.0f%%", ap1)), inherit.aes = FALSE,
            hjust = 1, family = "Arial", fontface = "bold", size = 5.8, colour = "#5B4300") +
  geom_text(data = C4, aes(cls, 100, label = sprintf("%.0f%%", partners)), inherit.aes = FALSE,
            hjust = 1.3, family = "Arial", fontface = "bold", size = 5.8, colour = "#44433F") +
  geom_text(data = C4, aes(cls, 2, label = sprintf("n = %d", n)), inherit.aes = FALSE,
            hjust = 0, family = "Arial", size = 4.1, colour = "#7A5E00") +
  annotate("text", x = 4.74, y = 0, label = "AP-1 motif itself", hjust = 0, family = "Arial",
           fontface = "bold", size = 4.8, colour = "#9A7400") +
  annotate("text", x = 4.74, y = 100, label = "partner motifs", hjust = 1, family = "Arial",
           fontface = "bold", size = 4.8, colour = "#6E6D69") +
  annotate("text", x = 4.74, y = 50, label = sprintf("Stable vs Dynamic, %s", fmtP(p_sh)), hjust = 0.5,
           family = "Arial", size = 4.3, colour = INK) +
  scale_fill_manual(values = c("AP-1 motif itself" = GOLD, "Partner motifs" = "#D8D7D3")) +
  scale_y_continuous(breaks = seq(0, 100, 25), labels = function(b) paste0(b, "%")) +
  scale_x_discrete(expand = expansion(add = c(0.52, 1.05))) +
  coord_flip(ylim = c(0, 100), expand = FALSE, clip = "off") +
  labs(x = NULL, y = "Share of the site's motif dependence (median per site)") +
  theme_pub() + no_cat_axis("y") +
  theme(legend.position = "none",
        axis.text.y = element_text(colour = rev(unname(COL3)), hjust = 0, size = 14.5),
        plot.margin = margin(6, 16, 4, 6))
save_fig(g4, "r_fb_share", 10.4, 4.0)
print(C4, row.names = FALSE)

# ==========================================================================================================
# 5. extra, not in the deck: where the single most important 10 bp tile sits. A clean null between the
#    classes, so it cannot carry a "Stable vs Dynamic" slide - but it is the strongest evidence that the
#    model found the motif on its own, which is what makes the composition above worth reading.
# ==========================================================================================================
zf <- file.path(ROOT, "Analysis/fb_scan/tables/top_tile_zone_by_class.csv")
tf <- file.path(ROOT, "Analysis/fb_scan/tables/top_tile_tests.csv")
if (file.exists(zf)) {
  Z <- read.csv(zf)
  TT <- if (file.exists(tf)) read.csv(tf) else NULL
  Z <- Z[Z$grouping == "grp4", ]
  Z$class <- factor(Z$class, levels = rev(ORD), labels = rev(gsub("\n", " ", CL[ORD])))
  ZN <- c("AP-1 core", "Flank, within 50 bp", "Beyond 50 bp")
  Z$zone <- factor(Z$zone, levels = ZN)
  Z <- Z[order(Z$class, Z$zone), ]
  core <- Z[Z$zone == "AP-1 core", ]
  chance <- 100 * 2 / 60
  p_core <- if (!is.null(TT)) TT$p[TT$comparison == "top tile is the AP-1 core, Stable vs Dynamic"][1] else NA

  g5 <- ggplot(Z, aes(class, pct, fill = zone)) +
    geom_col(width = 0.56, colour = "white", linewidth = 1.0, position = position_stack(reverse = TRUE)) +
    annotate("segment", x = 0.48, xend = 4.52, y = chance, yend = chance, linewidth = 0.45,
             linetype = "22", colour = INK) +
    geom_text(data = core, aes(class, pct, label = sprintf("%.0f%%", pct)), inherit.aes = FALSE,
              hjust = 1.35, family = "Arial", fontface = "bold", size = 5.8, colour = "#5B4300") +
    annotate("text", x = 4.74, y = chance, label = "chance: 2 of 60 tiles", hjust = 0, vjust = 0.5,
             family = "Arial", size = 4.2, colour = MUTED) +
    annotate("text", x = 4.74, y = 100, label = sprintf("Stable vs Dynamic, %s (n.s.)", fmtP(p_core)),
             hjust = 1, vjust = 0.5, family = "Arial", size = 4.3, colour = INK) +
    scale_fill_manual(values = c("AP-1 core" = GOLD, "Flank, within 50 bp" = GOLD_PALE,
                                 "Beyond 50 bp" = "#EDECE8")) +
    scale_y_continuous(breaks = seq(0, 100, 25), labels = function(b) paste0(b, "%")) +
    scale_x_discrete(expand = expansion(add = c(0.52, 1.05))) +
    coord_flip(ylim = c(0, 100), expand = FALSE, clip = "off") +
    labs(x = NULL, y = "Sites whose most important 10 bp tile falls in that zone") +
    theme_pub() + no_cat_axis("y") +
    theme(legend.position = "bottom", legend.justification = "left",
          axis.text.y = element_text(colour = rev(unname(COL3)), hjust = 0, size = 14.5),
          plot.margin = margin(6, 16, 4, 6))
  save_fig(g5, "r_fb_top_tile", 10.8, 4.1)
  print(core[, c("class", "n_sites", "pct", "lo", "hi")], row.names = FALSE)
} else {
  cat("  skipped r_fb_top_tile: run Scripts/56_fb_top_tile_rank.py on Katana first\n")
}

# ==========================================================================================================
# 6. importance as a density along the sequence - where, across the window, does a site's dependence sit?
#    Two panels in the same visual language as the MPRA deck's attention figure: the profile, and underneath
#    it the Stable - Dynamic difference, because the two curves are close enough that the eye cannot judge
#    the gap on its own. Cropped to +-150 bp; under 2% of either class's importance lies beyond that.
#    Left open at the top on purpose, so the talk can annotate it.
# ==========================================================================================================
pf <- file.path(ROOT, "Analysis/fb_scan/tables/importance_density_profile.csv")
cf <- file.path(ROOT, "Analysis/fb_scan/tables/importance_density_concentration.csv")
mf <- file.path(ROOT, "Analysis/fb_scan/tables/importance_density_matched.csv")
if (file.exists(pf)) {
  PR <- read.csv(pf)
  PR <- PR[PR$grouping == "stable_dynamic" & PR$version == "corrected", ]
  PR$cls <- factor(PR$class, levels = c("Stable", "Dynamic"))
  for (v in c("median", "lo", "hi")) PR[[v]] <- 100 * PR[[v]]
  XL <- 150
  PR <- PR[abs(PR$pos) <= XL, ]
  nn <- tapply(PR$n_sites, PR$cls, max)
  cc <- if (file.exists(cf)) read.csv(cf) else NULL
  mm <- if (file.exists(mf)) read.csv(mf) else NULL
  MOT <- c(-10, 10)                                  # the two tiles that overlap the 11 bp core
  COLD <- c(Stable = BLUE, Dynamic = ROSE)
  ytop <- max(PR$hi) * 1.30

  W <- reshape(PR[, c("pos", "class", "median", "fdr")], idvar = c("pos", "fdr"),
               timevar = "class", direction = "wide")
  names(W) <- sub("median\\.", "", names(W))
  W$diff <- W$Stable - W$Dynamic
  W$sig <- W$fdr < 0.05
  dlim <- max(abs(W$diff)) * 1.55

  top <- ggplot(PR, aes(pos, median, colour = cls, fill = cls)) +
    annotate("rect", xmin = MOT[1], xmax = MOT[2], ymin = 0, ymax = ytop, fill = GOLD_PALE) +
    annotate("text", x = 0, y = ytop, label = "AP-1 core", vjust = -0.3, family = "Arial",
             fontface = "bold", size = 4.5, colour = "#9A7400") +
    geom_area(aes(group = cls), position = "identity", alpha = 0.14, linewidth = 0, show.legend = FALSE) +
    geom_ribbon(aes(ymin = lo, ymax = hi, group = cls), alpha = 0.32, linewidth = 0, show.legend = FALSE) +
    geom_line(linewidth = 1.6, show.legend = FALSE) +
    annotate("text", x = -XL + 6, y = ytop * 0.97, label = sprintf("Stable  (n = %s)", nn[["Stable"]]),
             hjust = 0, vjust = 1, family = "Arial", fontface = "bold", size = 5.0, colour = BLUE) +
    annotate("text", x = -XL + 6, y = ytop * 0.83, label = sprintf("Dynamic  (n = %s)", nn[["Dynamic"]]),
             hjust = 0, vjust = 1, family = "Arial", fontface = "bold", size = 5.0, colour = ROSE) +
    scale_colour_manual(values = COLD) + scale_fill_manual(values = COLD) +
    scale_x_continuous(breaks = seq(-XL, XL, 50)) +
    scale_y_continuous(breaks = seq(0, 10, 5), labels = function(b) paste0(b, "%")) +
    coord_cartesian(xlim = c(-XL, XL), ylim = c(0, ytop * 1.1), expand = FALSE, clip = "off") +
    labs(x = NULL, y = "Share of importance\nper 10 bp") +
    theme_pub() + theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(),
                        axis.line.x = element_blank(), plot.margin = margin(16, 14, 0, 6))

  bot <- ggplot(W, aes(pos, diff)) +
    annotate("rect", xmin = MOT[1], xmax = MOT[2], ymin = -dlim, ymax = dlim, fill = GOLD_PALE) +
    geom_hline(yintercept = 0, linewidth = 0.4, colour = FAINT) +
    geom_segment(aes(xend = pos, yend = 0, colour = sig), linewidth = 1.9, show.legend = FALSE) +
    geom_point(aes(colour = sig), size = 1.9, show.legend = FALSE) +
    annotate("text", x = 0, y = max(W$diff), label = "FDR = 0.001", vjust = -0.75, hjust = 0.5,
             family = "Arial", size = 4.0, colour = INK) +
    annotate("text", x = XL - 6, y = -dlim * 0.72, label = "every other position: FDR > 0.8", hjust = 1,
             family = "Arial", size = 4.0, colour = MUTED) +
    annotate("text", x = -XL + 6, y = dlim * 0.80, label = "Stable − Dynamic (percentage points)", hjust = 0,
             family = "Arial", size = 4.0, colour = MUTED) +
    scale_colour_manual(values = c(`TRUE` = INK, `FALSE` = "#C6C5C1")) +
    scale_x_continuous(breaks = seq(-XL, XL, 50)) +
    scale_y_continuous(breaks = c(0, 1), labels = function(b) sprintf("%+.1f", b)) +
    coord_cartesian(xlim = c(-XL, XL), ylim = c(-dlim, dlim * 1.12), expand = FALSE, clip = "off") +
    labs(x = "Position relative to the AP-1 motif (bp)", y = NULL) +
    theme_pub() + theme(plot.margin = margin(2, 14, 4, 6))

  g6 <- top / bot + plot_layout(heights = c(2.4, 1.1))
  save_fig(g6, "r_fb_density", 11.0, 3.95)
  if (!is.null(cc)) print(cc[, c("class", "n_sites", "share_within_25bp",
                                 "p_stable_vs_dynamic_share", "p_share_matched")], row.names = FALSE)
  if (!is.null(mm)) cat(sprintf("     matched: Stable %.3f vs Dynamic %.3f, paired %+.4f, P = %.3g (%d pairs)\n",
                                mm$near_stable, mm$near_dynamic, mm$paired_difference, mm$p_wilcoxon, mm$n_pairs))
  cat(sprintf("     positions differing at FDR < 0.05: %s\n", paste(W$pos[W$sig], collapse = ", ")))
} else {
  cat("  skipped r_fb_density: run Scripts/57_fb_importance_density.py on Katana first\n")
}

cat("done\n")
