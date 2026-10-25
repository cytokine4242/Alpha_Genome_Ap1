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
  library(ggplot2); library(ragg); library(grid); library(patchwork); library(jsonlite)
})

ROOT <- "/Users/z5162865/WongLab/AlphaGenome_Ap1"
OUT  <- file.path(ROOT, "Analysis", "slide_figures")
set.seed(55)

INK <- "#111111"; MUTED <- "#4d4d4d"; RULE <- "#111111"; FAINT <- "#9a9993"
BLUE <- "#347CBA"            # Stable
ROSE <- "#C73A6B"            # Dynamic, pooled
ROSE_DK <- "#9E1F4C"; ROSE_LT <- "#E291AC"

# The five site classes, used in the same order and the same colours on every panel. The names are the labels as
# given (age_class / ex_class); Up means AP-1 gained, Down means AP-1 lost. A site that is dynamic with both age and
# exercise is counted under its age class, so the two exercise groups are the sites that are NOT also age-dynamic.
CLS5 <- c("Stable", "Age Up", "Age Down", "Exercise Up", "Exercise Down")
COL5 <- setNames(c(BLUE, "#C73A6B", "#8E1540", "#A85C7C", "#73314F"), CLS5)

site_classes <- function() {
  S <- read.csv(file.path(ROOT, "Katana/data/fb_subset_scan.csv"))[,
         c("site_id", "age_class", "ex_class", "stable_dynamic")]
  S <- S[!duplicated(S$site_id), ]
  S$cls <- ifelse(S$stable_dynamic == "Stable", "Stable",
           ifelse(S$age_class == "Age Up", "Age Up",
           ifelse(S$age_class == "Age Down", "Age Down",
           ifelse(S$ex_class == "Exercise Up", "Exercise Up",
           ifelse(S$ex_class == "Exercise Down", "Exercise Down", NA)))))
  S
}
# BH across the four comparisons each panel makes against Stable
fdr <- function(p) p.adjust(p, method = "BH")

# --- overlay geometry ------------------------------------------------------------------------------------
# Some annotation is better as PowerPoint text than as part of the PNG, so it can be moved, reworded or deleted
# on the slide. For that the deck has to know where ggplot put the panel inside the image. This measures the
# gtable and writes the panel rectangle (as fractions of the image) and the data ranges, together with the
# labels to draw and the data coordinates they belong at. build_fb_deck.js turns those into text boxes.
panel_geom <- function(g, W, H) {
  gt <- ggplotGrob(g)
  fill_null <- function(u, total) {                 # grid reports null units as 0; they share what is left over
    inch <- grid::convertWidth(u, "in", valueOnly = TRUE)
    isn <- as.character(grid::unitType(u)) == "null"
    if (any(isn)) {
      amt <- as.numeric(u)[isn]
      inch[isn] <- max(0, total - sum(inch[!isn])) * amt / sum(amt)
    }
    inch
  }
  wi <- fill_null(gt$widths, W); hi <- fill_null(gt$heights, H)
  p <- gt$layout[gt$layout$name == "panel", ][1, ]
  pp <- ggplot_build(g)$layout$panel_params[[1]]
  list(w_in = W, h_in = H,
       panel = list(left = sum(wi[seq_len(p$l - 1)]) / W, right = sum(wi[seq_len(p$r)]) / W,
                    top = sum(hi[seq_len(p$t - 1)]) / H, bottom = sum(hi[seq_len(p$b)]) / H),
       x_range = pp$x.range, y_range = pp$y.range)
}

write_overlay <- function(g, name, W, H, labels = list(), swatches = list()) {
  geom <- panel_geom(g, W, H)
  geom$name <- name
  geom$labels <- labels
  geom$swatches <- swatches
  jsonlite::write_json(geom, file.path(OUT, paste0(name, "_overlay.json")),
                       auto_unbox = TRUE, digits = 6, pretty = TRUE)
  cat(sprintf("  wrote %s_overlay.json (%d labels, %d swatches)\n", name, length(labels), length(swatches)))
}
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
td <- td[td$donor_class %in% CLS5, ]
td$cls <- factor(td$donor_class, levels = rev(CLS5))         # rev: ggplot draws the first level at the bottom
p_tr <- wilcox.test(lift ~ grp, data = td)$p.value
med <- aggregate(lift ~ cls, td, median)
nn  <- aggregate(lift ~ cls, td, length); names(nn)[2] <- "n"
labs_tr <- setNames(sprintf("%s", levels(td$cls)), levels(td$cls))
COL1 <- COL5

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
  theme(axis.text.y = element_text(colour = unname(rev(COL5)), hjust = 0, size = 15),
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
    labs(x = NULL, y = "Regression weight — both are positive, drawn apart only so the pair can be compared") +
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
SI <- site_classes()
J$cls <- SI$cls[match(J$site_id, SI$site_id)]
J <- J[!is.na(J$cls), ]
L <- rbind(
  data.frame(grp = J$cls, what = "Its own AP-1 motif", eff = J$own_ap1),
  data.frame(grp = J$cls, what = "Its partner motifs", eff = J$partners_mean_net))
L <- L[is.finite(L$eff), ]
L$grp <- factor(L$grp, levels = CLS5)
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
FILL3 <- c(setNames(paste0(COL5, "33"), paste(CLS5, levels(L$what)[1])),
           setNames(rep("white", length(CLS5)), paste(CLS5, levels(L$what)[2])))

pv <- function(w, cls) {
  a <- L$eff[L$grp == "Stable" & L$what == w]
  b <- L$eff[L$grp == cls & L$what == w]
  wilcox.test(a, b)$p.value
}
# four comparisons against Stable, so the panel prints the BH-adjusted value
q_own <- fdr(sapply(CLS5[-1], function(k) pv("Its own AP-1 motif", k)))
q_par <- fdr(sapply(CLS5[-1], function(k) pv("Its partner motifs", k)))
nlab <- sprintf("%s\n(n = %d)", CLS5, as.integer(table(L$grp[L$what == "Its own AP-1 motif"])))
dw <- 0.60
# the test sits over the box it belongs to, which five classes make far easier to follow than five brackets
PS <- data.frame(x = (seq_along(CLS5)) - dw / 4, y = 0.30,
                 lab = c("reference", sub("^P", "FDR", sapply(q_own, fmtP))))

g3 <- ggplot(B, aes(grp)) +
  geom_hline(yintercept = 0, linewidth = 0.35, colour = FAINT) +
  # with stat = "identity" the notch bounds have to be supplied; ggplot warns that it does not know the two
  # names but still carries them through to the geom, which is where they are read
  geom_boxplot(aes(ymin = ymin, lower = lower, middle = middle, upper = upper, ymax = ymax,
                   notchlower = notchlower, notchupper = notchupper, colour = grp, fill = key),
               stat = "identity", width = 0.46, position = position_dodge(width = dw),
               notch = TRUE, notchwidth = 0.58, linewidth = 0.6, show.legend = FALSE) +
  # the legend, the header and the FDR row are NOT drawn here: they go on the slide as PowerPoint text boxes
  # (see write_overlay below) so they can be moved, reworded or dropped without regenerating the figure
  scale_colour_manual(values = COL5) +
  scale_fill_manual(values = FILL3) +
  scale_x_discrete(labels = nlab, expand = expansion(add = c(0.95, 0.5))) +
  scale_y_continuous(breaks = seq(-2.5, 0, 0.5)) +
  coord_cartesian(ylim = c(min(B$ymin) - 0.12, 1.14), expand = FALSE, clip = "on") +
  labs(x = NULL, y = expression(paste("Effect of shuffling it (log"[2], ")"))) +
  theme_pub() + theme(legend.position = "none",
                      axis.text.x = element_text(colour = unname(COL5), size = 13, lineheight = 0.95))
save_fig(g3, "r_fb_specificity", 10.8, 4.0)
# every piece of text that was overlaid on the panel, handed to the deck as data instead
write_overlay(g3, "r_fb_specificity", 10.8, 4.0,
  labels = c(
    list(list(id = "header", x = 0.80, y = 0.80, text = "vs Stable, at the site's own AP-1 motif (BH across the four)",
              size_pt = 12, color = "4D4D4D", anchor = "start")),
    lapply(seq_len(nrow(PS)), function(i)
      list(id = paste0("fdr", i), x = PS$x[i], y = PS$y[i], text = PS$lab[i],
           size_pt = 12, color = "4D4D4D", anchor = "middle")),
    list(list(id = "leg1", x = 3.03, y = 0.95, text = "its own AP-1 motif",
              size_pt = 12.5, color = "111111", anchor = "start"),
         list(id = "leg2", x = 3.03, y = 0.70, text = "its partner motifs, net of a matched control",
              size_pt = 12.5, color = "111111", anchor = "start"))),
  swatches = list(
    list(id = "sw1", x = 2.95, y = 0.95, fill = "D6D6D6", line = "111111"),
    list(id = "sw2", x = 2.95, y = 0.70, fill = "FFFFFF", line = "111111")))
cat(sprintf("     own AP-1 vs Stable, BH: %s\n", paste(sprintf("%s %.3g", CLS5[-1], q_own), collapse = ", ")))
cat(sprintf("     partners vs Stable, BH: %s\n", paste(sprintf("%s %.3g", CLS5[-1], q_par), collapse = ", ")))

# ==========================================================================================================
# 4. how a site divides its motif dependence - a composition, split by class
#    The two parts are complementary, so a median split sums to 100% and the bar is an honest whole.
# ==========================================================================================================
J$share <- abs(J$own_ap1) / (abs(J$own_ap1) + abs(J$partners_mean_net))
K <- J[is.finite(J$share), ]
K$cls <- factor(K$cls, levels = rev(CLS5))
K$sd2 <- factor(ifelse(K$cls == "Stable", "Stable", "Dynamic"), levels = c("Stable", "Dynamic"))
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
  geom_text(data = C4, aes(cls, ap1 - 7.5, label = sprintf("%.0f%%", ap1)), inherit.aes = FALSE,
            hjust = 1, family = "Arial", fontface = "bold", size = 5.8, colour = "#5B4300") +
  geom_text(data = C4, aes(cls, 100, label = sprintf("%.0f%%", partners)), inherit.aes = FALSE,
            hjust = 1.3, family = "Arial", fontface = "bold", size = 5.8, colour = "#44433F") +
  geom_text(data = C4, aes(cls, 2, label = sprintf("n = %d", n)), inherit.aes = FALSE,
            hjust = 0, family = "Arial", size = 4.1, colour = "#7A5E00") +
  annotate("text", x = 5.74, y = 0, label = "AP-1 motif itself", hjust = 0, family = "Arial",
           fontface = "bold", size = 4.8, colour = "#9A7400") +
  annotate("text", x = 5.74, y = 100, label = "partner motifs", hjust = 1, family = "Arial",
           fontface = "bold", size = 4.8, colour = "#6E6D69") +
  annotate("text", x = 5.74, y = 50, label = sprintf("every class below Stable, %s pooled", fmtP(p_sh)),
           hjust = 0.5, family = "Arial", size = 4.3, colour = INK) +
  scale_fill_manual(values = c("AP-1 motif itself" = GOLD, "Partner motifs" = "#D8D7D3")) +
  scale_y_continuous(breaks = seq(0, 100, 25), labels = function(b) paste0(b, "%")) +
  scale_x_discrete(expand = expansion(add = c(0.52, 1.05))) +
  coord_flip(ylim = c(0, 100), expand = FALSE, clip = "off") +
  labs(x = NULL, y = "Share of the site's motif dependence (median per site)") +
  theme_pub() + no_cat_axis("y") +
  theme(legend.position = "none",
        axis.text.y = element_text(colour = rev(unname(COL5)), hjust = 0, size = 14.5),
        plot.margin = margin(6, 16, 4, 6))
save_fig(g4, "r_fb_share", 10.4, 4.3)
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
  ZN0 <- c("AP-1 core", "Flank, within 50 bp", "Beyond 50 bp")
  TS <- read.csv(file.path(ROOT, "Analysis/fb_scan/tables/top_tile_per_site.csv"))
  TS$cls <- ifelse(TS$stable_dynamic == "Stable", "Stable",
            ifelse(TS$age_class == "Age Up", "Age Up",
            ifelse(TS$age_class == "Age Down", "Age Down",
            ifelse(TS$ex_class == "Exercise Up", "Exercise Up",
            ifelse(TS$ex_class == "Exercise Down", "Exercise Down", NA)))))
  TS <- TS[!is.na(TS$cls), ]
  Z <- do.call(rbind, lapply(CLS5, function(g) {
    v <- TS$top_zone[TS$cls == g]
    data.frame(class = g, n_sites = length(v), zone = ZN0, pct = 100 * sapply(ZN0, function(z) mean(v == z)))
  }))
  Z$class <- factor(Z$class, levels = rev(CLS5))
  ZN <- ZN0
  Z$zone <- factor(Z$zone, levels = ZN)
  Z <- Z[order(Z$class, Z$zone), ]
  core <- Z[Z$zone == "AP-1 core", ]
  chance <- 100 * 2 / 60
  p_core <- if (!is.null(TT)) TT$p[TT$comparison == "top tile is the AP-1 core, Stable vs Dynamic"][1] else NA

  g5 <- ggplot(Z, aes(class, pct, fill = zone)) +
    geom_col(width = 0.56, colour = "white", linewidth = 1.0, position = position_stack(reverse = TRUE)) +
    annotate("segment", x = 0.48, xend = 5.52, y = chance, yend = chance, linewidth = 0.45,
             linetype = "22", colour = INK) +
    geom_text(data = core, aes(class, pct, label = sprintf("%.0f%%", pct)), inherit.aes = FALSE,
              hjust = 1.35, family = "Arial", fontface = "bold", size = 5.8, colour = "#5B4300") +
    annotate("text", x = 5.74, y = chance, label = "chance: 2 of 60 tiles", hjust = 0, vjust = 0.5,
             family = "Arial", size = 4.2, colour = MUTED) +
    annotate("text", x = 5.74, y = 100, label = sprintf("Stable vs Dynamic, %s (n.s.)", fmtP(p_core)),
             hjust = 1, vjust = 0.5, family = "Arial", size = 4.3, colour = INK) +
    scale_fill_manual(values = c("AP-1 core" = GOLD, "Flank, within 50 bp" = GOLD_PALE,
                                 "Beyond 50 bp" = "#EDECE8")) +
    scale_y_continuous(breaks = seq(0, 100, 25), labels = function(b) paste0(b, "%")) +
    scale_x_discrete(expand = expansion(add = c(0.52, 1.05))) +
    coord_flip(ylim = c(0, 100), expand = FALSE, clip = "off") +
    labs(x = NULL, y = "Sites whose most important 10 bp tile falls in that zone") +
    theme_pub() + no_cat_axis("y") +
    theme(legend.position = "bottom", legend.justification = "left",
          axis.text.y = element_text(colour = rev(unname(COL5)), hjust = 0, size = 14.5),
          plot.margin = margin(6, 16, 4, 6))
  save_fig(g5, "r_fb_top_tile", 10.8, 4.4)
  print(core[, c("class", "n_sites", "pct")], row.names = FALSE)
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
# the sliding scan supersedes the 60-tile one where it exists: same quantity, sampled every 2 bp instead of every 10
sf <- file.path(ROOT, "Analysis/fb_scan/tables/slide_density_profile.csv")
pf <- file.path(ROOT, "Analysis/fb_scan/tables/importance_density_profile.csv")
cf <- file.path(ROOT, "Analysis/fb_scan/tables/importance_density_concentration.csv")
mf <- file.path(ROOT, "Analysis/fb_scan/tables/importance_density_matched.csv")
SLIDE <- file.exists(sf)
if (SLIDE) {
  pf <- sf
  cf <- file.path(ROOT, "Analysis/fb_scan/tables/slide_density_concentration.csv")
  mf <- cf
}
if (file.exists(pf)) {
  PR <- read.csv(pf)
  if (!SLIDE) PR <- PR[PR$grouping == "stable_dynamic" & PR$version == "corrected", ]
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
    { if (SLIDE && any(W$sig))
        annotate("rect", xmin = min(W$pos[W$sig]), xmax = max(W$pos[W$sig]), ymin = -dlim, ymax = dlim,
                 fill = INK, alpha = 0.07) } +
    { if (SLIDE) geom_line(linewidth = 0.9, colour = "#55544F")
      else geom_segment(aes(xend = pos, yend = 0, colour = sig), linewidth = 1.9, show.legend = FALSE) } +
    { if (SLIDE) geom_point(data = W[W$sig, ], size = 0.9, colour = INK)
      else geom_point(aes(colour = sig), size = 1.9, show.legend = FALSE) } +
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

# ==========================================================================================================
# 7. three other readings of the same density, for picking one. Nothing here changes the measurement; what
#    changes is what is plotted. Each panel prints the effect size so a transformed axis can be checked.
# ==========================================================================================================
cuf <- file.path(ROOT, "Analysis/fb_scan/tables/density_cumulative.csv")
dif <- file.path(ROOT, "Analysis/fb_scan/tables/density_difference.csv")
c5f <- file.path(ROOT, "Analysis/fb_scan/tables/density_concentration_5class.csv")
if (file.exists(cuf)) {
  CU <- read.csv(cuf); for (v in c("median", "lo", "hi")) CU[[v]] <- 100 * CU[[v]]

  # ---- variation 1: cumulative. Per-position differences of a point each add up to a visible gap ----------
  PO <- CU[CU$grouping == "pooled", ]
  PO$cls <- factor(PO$class, levels = c("Stable", "Dynamic"))
  w <- reshape(PO[, c("radius_bp", "class", "median")], idvar = "radius_bp",
               timevar = "class", direction = "wide")
  names(w) <- sub("median\\.", "", names(w))
  w$gap <- w$Stable - w$Dynamic
  gx <- w$radius_bp[which.max(w$gap)]; gy <- max(w$gap)
  half <- function(cl) approx(w[[cl]], w$radius_bp, xout = 50)$y
  hS <- half("Stable"); hD <- half("Dynamic")
  nn <- tapply(PO$n_sites, PO$cls, max)

  g7 <- ggplot(PO, aes(radius_bp, median, colour = cls, fill = cls)) +
    geom_ribbon(aes(ymin = lo, ymax = hi), alpha = 0.22, linewidth = 0, show.legend = FALSE) +
    geom_line(linewidth = 1.7, show.legend = FALSE) +
    annotate("segment", x = gx, xend = gx, y = w$Dynamic[w$radius_bp == gx], yend = w$Stable[w$radius_bp == gx],
             linewidth = 0.9, colour = INK) +
    annotate("point", x = c(gx, gx), y = c(w$Dynamic[w$radius_bp == gx], w$Stable[w$radius_bp == gx]),
             size = 2.1, colour = INK) +
    annotate("text", x = gx + 12, y = mean(c(w$Dynamic[w$radius_bp == gx], w$Stable[w$radius_bp == gx])),
             label = sprintf("%.1f points apart\nat ±%d bp", gy, gx), hjust = 0, vjust = 0.5,
             family = "Arial", size = 4.6, colour = INK, lineheight = 0.95) +
    annotate("segment", x = 0, xend = max(hS, hD), y = 50, yend = 50, linewidth = 0.4,
             linetype = "22", colour = FAINT) +
    annotate("text", x = 8, y = 101, label = sprintf("Stable  (n = %s)", nn[["Stable"]]), hjust = 0,
             vjust = 1, family = "Arial", fontface = "bold", size = 5, colour = BLUE) +
    annotate("text", x = 8, y = 92, label = sprintf("Dynamic  (n = %s)", nn[["Dynamic"]]), hjust = 0,
             vjust = 1, family = "Arial", fontface = "bold", size = 5, colour = ROSE) +
    scale_colour_manual(values = c(Stable = BLUE, Dynamic = ROSE)) +
    scale_fill_manual(values = c(Stable = BLUE, Dynamic = ROSE)) +
    scale_x_continuous(breaks = seq(0, 300, 50)) +
    scale_y_continuous(breaks = seq(0, 100, 25), labels = function(b) paste0(b, "%")) +
    coord_cartesian(xlim = c(0, 300), ylim = c(0, 104), expand = FALSE, clip = "off") +
    labs(x = "Distance from the AP-1 motif (± bp)",
         y = "Share of the site's importance\nwithin that distance") +
    theme_pub() + theme(plot.margin = margin(14, 16, 4, 6))
  g7 <- g7 + annotate("text", x = 300, y = 6,
                      label = sprintf("half of it lies within ±%.0f bp (Stable) and ±%.0f bp (Dynamic)", hS, hD),
                      hjust = 1, family = "Arial", size = 4.2, colour = MUTED)
  save_fig(g7, "r_fb_density_cumulative", 9.4, 4.4)
  cat(sprintf("     cumulative: max gap %.1f pp at ±%d bp; half-radius %.0f vs %.0f bp\n", gy, gx, hS, hD))

  # ---- variation 2: the difference on its own axis -------------------------------------------------------
  DF <- read.csv(dif)
  DF$sig <- DF$fdr < 0.05
  g8 <- ggplot(DF, aes(pos, diff_pp)) +
    annotate("rect", xmin = -10, xmax = 10, ymin = -Inf, ymax = Inf, fill = GOLD_PALE) +
    geom_hline(yintercept = 0, linewidth = 0.4, colour = FAINT) +
    geom_ribbon(aes(ymin = lo, ymax = hi), alpha = 0.20, linewidth = 0, fill = BLUE) +
    geom_line(linewidth = 1.4, colour = BLUE) +
    geom_point(data = DF[DF$sig, ], size = 2.6, colour = INK) +
    annotate("text", x = 0, y = max(DF$hi) * 1.04, label = "AP-1 core", vjust = -0.2, family = "Arial",
             fontface = "bold", size = 4.5, colour = "#9A7400") +
    annotate("text", x = 300, y = min(DF$lo), label = "filled points: FDR < 0.05", hjust = 1, vjust = 0,
             family = "Arial", size = 4.2, colour = MUTED) +
    scale_x_continuous(breaks = seq(-300, 300, 100)) +
    coord_cartesian(xlim = c(-300, 300), ylim = c(min(DF$lo) * 1.1, max(DF$hi) * 1.28),
                    expand = FALSE, clip = "off") +
    labs(x = "Position relative to the AP-1 motif (bp)",
         y = "Stable − Dynamic\n(percentage points of a site's importance)") +
    theme_pub() + theme(plot.margin = margin(16, 14, 4, 6))
  save_fig(g8, "r_fb_density_diff", 9.4, 4.2)

  # ---- variation 3: the five classes ---------------------------------------------------------------------
  CL <- CU[CU$grouping == "class", ]
  CL$cls <- factor(CL$class, levels = CLS5)
  C5 <- read.csv(c5f)
  C5$cls <- factor(C5$class, levels = rev(CLS5))
  for (v in c("share_within_25bp", "lo", "hi")) C5[[v]] <- 100 * C5[[v]]
  C5$lab <- sapply(C5$fdr_vs_stable, function(p) if (is.na(p)) "reference" else sub("^P", "FDR", fmtP(p)))
  g9 <- ggplot(C5, aes(cls, share_within_25bp, colour = cls)) +
    geom_errorbar(aes(ymin = lo, ymax = hi), width = 0, linewidth = 1.5) +
    geom_point(size = 4.2) +
    geom_text(aes(label = sprintf("%.1f%%", share_within_25bp)), hjust = -0.45, vjust = -0.75,
              family = "Arial", fontface = "bold", size = 4.8, show.legend = FALSE) +
    geom_text(aes(y = 20.6, label = lab), hjust = 0, vjust = -0.9, family = "Arial", size = 4.2,
              colour = MUTED) +
    geom_hline(yintercept = C5$share_within_25bp[C5$class == "Stable"], linewidth = 0.4,
               linetype = "22", colour = BLUE) +
    scale_colour_manual(values = COL5) +
    scale_y_continuous(breaks = seq(20, 40, 5), labels = function(b) paste0(b, "%")) +
    scale_x_discrete(expand = expansion(add = 0.7),
                     labels = function(x) sprintf("%s\n(n = %d)", x, C5$n_sites[match(x, C5$class)])) +
    coord_flip(xlim = c(0.4, 5.8), ylim = c(20, 40), expand = FALSE, clip = "off") +
    labs(x = NULL, y = "Share of the site's importance within ±25 bp of the motif") +
    theme_pub() + no_cat_axis("y") +
    theme(legend.position = "none",
          axis.text.y = element_text(colour = rev(unname(COL5)), hjust = 0, size = 14),
          plot.margin = margin(10, 20, 4, 6))
  save_fig(g9, "r_fb_density_classes", 9.4, 4.4)
  print(C5[, c("class", "n_sites", "share_within_25bp", "lo", "hi", "fdr_vs_stable")], row.names = FALSE)
} else {
  cat("  skipped the density variations: run Scripts/59_fb_density_variations.py on Katana first\n")
}
