const pptxgen = require("pptxgenjs");
const path = require("path");
const { applyTheme } = require("/Users/z5162865/.claude/skills/synced/2361eaa7-acd4-463c-8e26-dee27087b623_1ef5f431-f760-4047-8679-39ac25a5dd75/pptx/scripts/apply_theme.js");

const ROOT = "/Users/z5162865/WongLab/AlphaGenome_Ap1";
const F = (p) => path.join(ROOT, p);

const THEME = {
  name: "AP1 Mono",
  headFontFace: "Cambria",
  bodyFontFace: "Calibri",
  colors: {
    dk1: "111111", lt1: "FFFFFF", dk2: "595959", lt2: "F2F2F2",
    accent1: "C2551C", accent2: "1F6FB2", accent3: "111111", accent4: "8C8C8C",
    accent5: "D9D9D9", accent6: "595959", hlink: "1F6FB2", folHlink: "595959",
  },
};

const pres = new pptxgen();
pres.layout = "LAYOUT_WIDE";            // 13.3 x 7.5
pres.theme = { headFontFace: THEME.headFontFace, bodyFontFace: THEME.bodyFontFace };
pres.author = "Jack";
pres.title = "AlphaGenome perturbation of AP-1 sites";
const C = pres.SchemeColor;

// ------------------------------------------------------------------ layouts
pres.defineSlideMaster({
  title: "TITLE", background: { color: "FFFFFF" },
  objects: [
    { placeholder: { options: { name: "title", type: "title", x: 0.9, y: 2.3, w: 11.5, h: 1.5, fontSize: 40, bold: true, color: C.text1, align: "left", valign: "bottom" }, text: " " } },
    { placeholder: { options: { name: "body", type: "body", x: 0.9, y: 4.05, w: 11.5, h: 1.5, fontSize: 17, color: C.text2, align: "left", valign: "top" }, text: " " } },
  ],
});

pres.defineSlideMaster({
  title: "TITLE_FIGURE", background: { color: "FFFFFF" },
  objects: [
    { placeholder: { options: { name: "title", type: "title", x: 0.55, y: 0.33, w: 12.2, h: 0.75, fontSize: 28, bold: true, color: C.text1, align: "left", valign: "middle" }, text: " " } },
    { placeholder: { options: { name: "body", type: "body", x: 0.55, y: 6.5, w: 12.2, h: 0.72, fontSize: 15, color: C.text1, align: "left", valign: "top" }, text: " " } },
  ],
  slideNumber: { x: 12.75, y: 7.05, fontSize: 10, color: "A6A6A6" },
});

pres.defineSlideMaster({
  title: "TITLE_BODY", background: { color: "FFFFFF" },
  objects: [
    { placeholder: { options: { name: "title", type: "title", x: 0.55, y: 0.33, w: 12.2, h: 0.75, fontSize: 28, bold: true, color: C.text1, align: "left", valign: "middle" }, text: " " } },
    { placeholder: { options: { name: "body", type: "body", x: 0.55, y: 1.35, w: 12.2, h: 5.6, fontSize: 17, color: C.text1, align: "left", valign: "top" }, text: " " } },
  ],
  slideNumber: { x: 12.75, y: 7.05, fontSize: 10, color: "A6A6A6" },
});

// ------------------------------------------------------------------ helpers
const sizeOf = require("image-size").default || require("image-size");
// fit an image inside a box, preserving its real aspect ratio — stretching is the commonest slide defect
function FIG(slide, file, box) {
  const { width, height } = sizeOf(require("fs").readFileSync(F(file)));
  const r = width / height;
  let w = box.w, h = w / r;
  if (h > box.h) { h = box.h; w = h * r; }
  slide.addImage({ path: F(file), x: box.x + (box.w - w) / 2, y: box.y + (box.h - h) / 2, w, h });
}
const BOX = { x: 0.55, y: 1.2, w: 12.2, h: 5.1 };
// a slide that opens with a method band: the band sits under the title, the figure takes what is left
const BAND = { x: 0.55, y: 1.05, w: 7.9, h: 1.1 };
const BOXB = { x: 0.55, y: 2.22, w: 12.2, h: 4.0 };

// Lay out a class figure from its manifest: every panel is its own image you can move or replace, and every
// word on the slide is a real text box, so wording changes never mean regenerating a figure.
function legendRow(slide, items, x, y) {
  let lx = x;
  items.forEach(([lab, col], i) => {
    slide.addShape(pres.ShapeType.ellipse, { x: lx, y: y + 0.045, w: 0.12, h: 0.12, fill: { color: col }, line: { color: col, width: 0 }, objectName: `lg-dot-${y}-${i}` });
    const w = 0.22 + lab.length * 0.072;
    slide.addText(lab, { isTextBox: true, x: lx + 0.17, y, w, h: 0.22, fontSize: 10.5, color: "111111", margin: 0, objectName: `lg-txt-${y}-${i}` });
    lx += 0.17 + w + 0.12;
  });
}

function classFigure(slide, name, box, opts = {}) {
  const man = JSON.parse(require("fs").readFileSync(F(`Analysis/slide_figures/${name}.json`), "utf8"));
  const nCol = man.cols.length, nRow = man.panels[0].length;
  const stack = opts.stack ?? false;                       // one measure per row instead of per column
  const gut = opts.gutter ?? 0.12, head = stack ? 0.0 : 0.82, foot = opts.foot ?? 0.92, gap = 0.08;
  const cw = stack ? box.w - gut : (box.w - gut) / nCol;
  const ch = stack ? (box.h - foot) / nCol : (box.h - head - foot) / nRow;
  man.cols.forEach((c, ci) => {
    const x = box.x + gut;
    if (stack) {
      const y0 = box.y + ci * ch;
      slide.addText(c.title, { isTextBox: true, x, y: y0, w: cw, h: 0.26, fontSize: 11.5, bold: true, color: "111111", margin: 0, objectName: `${name}-title-${ci}` });
      FIG(slide, `Analysis/slide_figures/panels/${man.panels[ci][0]}`, { x, y: y0 + 0.28, w: cw, h: ch - 0.34 });
      return;
    }
    slide.addText(c.title, { isTextBox: true, x: x + ci * cw, y: box.y, w: cw - gap, h: 0.48, fontSize: 12.5, bold: true, color: "111111", valign: "top", margin: 0, objectName: `${name}-title-${ci}` });
    slide.addText(c.subtitle, { isTextBox: true, x: x + ci * cw, y: box.y + 0.47, w: cw - gap, h: 0.34, fontSize: 9.5, color: "595959", valign: "top", margin: 0, objectName: `${name}-sub-${ci}` });
    man.panels[ci].forEach((f, ri) => {
      FIG(slide, `Analysis/slide_figures/panels/${f}`, { x: x + ci * cw, y: box.y + head + ri * ch, w: cw - gap, h: ch - gap });
    });
  });
  if (nRow > 1) man.rows.forEach((r, ri) => {
    slide.addText(r, { isTextBox: true, x: box.x - 0.52, y: box.y + head + ri * ch + ch / 2 - 0.75, w: 1.5, h: 1.5,
      fontSize: 10, color: "111111", align: "center", valign: "middle", rotate: 270, margin: 0, objectName: `${name}-row-${ri}` });
  });
  const lx0 = opts.legendX ?? box.x;
  let lx = lx0;
  const lw = opts.legendW ?? box.w;
  const ly = opts.legendY ?? (box.y + box.h - foot + 0.10);
  man.legend.forEach(([lab, col], i) => {
    slide.addShape(pres.ShapeType.ellipse, { x: lx, y: ly + 0.045, w: 0.12, h: 0.12, fill: { color: col.replace("#", "") }, line: { color: col.replace("#", ""), width: 0 }, objectName: `${name}-dot-${i}` });
    const w = 0.22 + lab.length * 0.072;
    slide.addText(lab, { isTextBox: true, x: lx + 0.17, y: ly, w, h: 0.22, fontSize: 10.5, color: "111111", margin: 0, objectName: `${name}-leg-${i}` });
    lx += 0.17 + w + 0.12;
  });
  if (opts.showKey !== false) slide.addText(man.key, { isTextBox: true, x: lx0, y: ly + 0.30, w: lw, h: 0.55, fontSize: 9, color: "595959", margin: 0, objectName: `${name}-key` });
}

function statRow(slide, items, y, h = 1.25) {
  const n = items.length, gap = 0.3, w = (12.2 - gap * (n - 1)) / n;
  items.forEach((it, i) => {
    const x = 0.55 + i * (w + gap);
    slide.addShape(pres.ShapeType.roundRect, { x, y, w, h, rectRadius: 0.06, fill: { color: "F2F2F2" }, line: { color: "E2E2E2", width: 1 }, objectName: `stat-${i}` });
    slide.addText(it.v, { isTextBox: true, x: x + 0.12, y: y + 0.1, w: w - 0.24, h: 0.6, fontSize: 30, bold: true, color: it.c || "111111", align: "center", margin: 0 });
    slide.addText(it.l, { isTextBox: true, x: x + 0.12, y: y + 0.72, w: w - 0.24, h: h - 0.8, fontSize: 12, color: "595959", align: "center", margin: 0 });
  });
}

function takeaway(slide, text) {
  slide.addText([{ text: "Takeaway   ", options: { bold: true, color: "C2551C" } }, { text, options: { color: "111111" } }],
    { placeholder: "body" });
}

// ================================================================== 1 title
pres.addSection({ title: "Opening" });
let s = pres.addSlide({ masterName: "TITLE", sectionTitle: "Opening" });
s.addText("What makes an AP-1 site respond to ageing?", { placeholder: "title" });
s.addText([
  { text: "In-silico mutagenesis of every cardiomyocyte AP-1 element with AlphaGenome", options: { breakLine: true } },
  { text: "3,678 sites · 364 mutations each · validated against snATAC, RNA and a fibroblast MPRA", options: { color: "595959", fontSize: 15 } },
], { placeholder: "body" });
s.addShape(pres.ShapeType.line, { x: 0.98, y: 3.82, w: 2.2, h: 0, line: { color: "111111", width: 2.5 }, objectName: "rule" });
s.addNotes("Mouse cardiomyocytes, young sedentary / aged sedentary / aged exercise. Sites defined by TOBIAS footprints of Fos motifs.");

// ================================================================== 2 design
pres.addSection({ title: "Design" });
s = pres.addSlide({ masterName: "TITLE_BODY", sectionTitle: "Design" });
s.addText("Two ways to call a site responsive, and they disagree", { placeholder: "title" });
s.addText([
  { text: "Threshold (binary call)", options: { bold: true, fontSize: 19, breakLine: true } },
  { text: "The TOBIAS footprint crosses the calling threshold between conditions. Age Up 396, Age Down 179, Exercise Up 186, Exercise Down 276. Stable = called in all three conditions, 2,907 sites.", options: { breakLine: true, paraSpaceAfter: 14 } },
  { text: "Δscore (graded)", options: { bold: true, fontSize: 19, breakLine: true } },
  { text: "The top fifth of the change in binding score, whether or not the call ever crossed the line. 736 sites each way.", options: { breakLine: true, paraSpaceAfter: 14 } },
  { text: "Only 22% of Δscore-Up sites ever cross the threshold, and 40% of threshold-Stable sites show a large score change.", options: { color: "595959" } },
], { isTextBox: true, x: 0.55, y: 1.35, w: 6.3, h: 4.0, fontSize: 17, color: "111111", margin: 0 });
statRow(s, [
  { v: "ρ 0.57", l: "Δscore vs the measured snATAC change with age", c: "C2551C" },
  { v: "ρ 0.31", l: "threshold class vs the same measured change", c: "1F6FB2" },
], 5.55, 1.35);
s.addText("Every comparison below matches each responsive site to a Stable site on predicted accessibility, peak width, TSS distance and measured accessibility.",
  { isTextBox: true, x: 7.1, y: 1.5, w: 5.65, h: 1.5, fontSize: 15, color: "111111", fill: { color: "F2F2F2" }, margin: 12, objectName: "matching-note" });
s.addText([
  { text: "Which to believe?", options: { bold: true, fontSize: 17, breakLine: true } },
  { text: "Δscore tracks the measured snATAC change about twice as well, so it is the better description of what actually changed. The threshold call is what the field reports, so both are carried through and every result below is labelled with the one it comes from.", options: { fontSize: 15, color: "595959" } },
], { isTextBox: true, x: 7.1, y: 3.25, w: 5.65, h: 2.1, margin: 0 });
s.addNotes("The two definitions pick out largely different sites. Report every result labelled by which definition it comes from.");

// ================================================================== 3 method schematic
s = pres.addSlide({ masterName: "TITLE_FIGURE", sectionTitle: "Design" });
s.addText("What one site looks like", { placeholder: "title" });
FIG(s, "Analysis/slide_figures/s0_element.png", { x: 0.55, y: 1.25, w: 12.2, h: 3.85 });
statRow(s, [
  { v: "3,678", l: "cardiomyocyte AP-1 sites" },
  { v: "364", l: "edits of each one" },
  { v: "772 / 2,907", l: "Dynamic against Stable sites", c: "C2551C" },
], 5.3, 0.95);
takeaway(s, "every edit is read as the change in predicted accessibility over the surrounding 170 bp");
s.addNotes("The 131 kb window is what AlphaGenome takes as input; we only ever edit the middle 600 bp. Each tile is 3 replicate random replacements, averaged. Effects are log2(edited / native).");

// ================================================================== 4 the null
// ================================================================== 3b the method at one real element
s = pres.addSlide({ masterName: "TITLE_FIGURE", sectionTitle: "Design" });
s.addText("The same thing at one real element", { placeholder: "title" });
s.addText("Every 10 bp window from \u2212300 to +300 bp is replaced by random bases in turn, three times, and the bar is what the prediction loses over the surrounding 170 bp. Gold is the called AP-1 motif. The lane below is every other footprinted motif TOBIAS calls in this element, which is what the windows are read against.",
  { isTextBox: true, x: 0.8, y: 1.12, w: 11.9, h: 0.52, fontSize: 11.5, color: "595959", margin: 0, objectName: "ml-intro" });
FIG(s, "Analysis/slide_figures/panels/s11_method_locus.png", { x: 0.55, y: 1.9, w: 12.2, h: 3.5 });
s.addText("chr2:167,698,004 \u00b7 the motif window costs \u22120.56 against a median of \u22120.66 across all 3,678 sites, and the deepest window is a neighbour \u2014 as it is at 69% of sites",
  { isTextBox: true, x: 0.8, y: 5.55, w: 11.9, h: 0.3, fontSize: 10.5, color: "111111", margin: 0, objectName: "ml-cap" });
takeaway(s, "the motif is one window among sixty, and usually not the deepest \u2014 everything that follows is this, repeated over 3,678 sites");
s.addNotes("A demonstration, not a result, but chosen to be typical rather than striking: the motif window, the total dependence and the deepest-neighbour pattern all sit near the median. An earlier version used a site in the top 8% for motif depth, which made the motif look far more dominant than it is. The footprint lane is every named family TOBIAS calls here, with overlapping calls of one family merged and counted; only the two catch-all bins (Other, Zinc finger other) are dropped. The same scan is also read over the 11 bp motif and over the whole peak, and the titration and single-base experiments run on the same 3,678 sites.");

pres.addSection({ title: "Findings" });
s = pres.addSlide({ masterName: "TITLE_FIGURE", sectionTitle: "Findings" });
s.addText("Responsive sites lose less to every edit of their element", { placeholder: "title" });
s.addText("The three edits, and the control", { isTextBox: true, x: 0.55, y: 1.08, w: 5.4, h: 0.28, fontSize: 12.5, bold: true, color: "111111", margin: 0, objectName: "ed-t" });
FIG(s, "Analysis/slide_figures/panels/d_edits.png", { x: 0.55, y: 1.38, w: 5.4, h: 1.3 });
classFigure(s, "s2_matched_null", { x: 0.55, y: 2.82, w: 12.2, h: 3.36 }, { foot: 1.0 });
takeaway(s, "every class loses less than Stable sites do at every scale, and the 300 bp control shows it is not a global offset");
s.addNotes("Top row is the uncontrolled comparison: if the accessibility difference is itself biology, it lives there. Bottom row matches each class to its own Stable controls on predicted and measured accessibility, peak width and TSS distance. The 300 bp control is the thing to check: Age Up differs slightly even there (-0.152 vs -0.111 matched), but that offset is 0.04 against a 0.26 effect at 170 bp and runs the opposite way, so it does not explain the result.");

// ================================================================== 4 partners
s = pres.addSlide({ masterName: "TITLE_FIGURE", sectionTitle: "Findings" });
s.addText("Both classes lean on the same partners", { placeholder: "title" });
FIG(s, "Analysis/slide_figures/d_tile.png", { x: 0.55, y: 1.05, w: 7.9, h: 1.1 });
FIG(s, "Analysis/slide_figures/panels/s8_partners_all.png", { x: 3.0, y: 2.3, w: 7.3, h: 3.35 });
legendRow(s, [["Stable", "B3B2AE"], ["Age Up", "C2551C"], ["Age Down", "1F6FB2"], ["Exercise Up", "2E7D32"], ["Exercise Down", "7B3FA0"]], 1.4, 5.72);
s.addText("six families shown; tests run over all 20 class-by-family comparisons, smallest FDR 0.02",
  { isTextBox: true, x: 1.4, y: 6.02, w: 11.1, h: 0.3, fontSize: 9.5, color: "595959", margin: 0, objectName: "pa-key" });
takeaway(s, "the ranking is the validation \u2014 MEF2, GATA, NFI, TEAD, none of it told to the model; the class split is the question, and no class leans on a different partner");
s.addNotes("2.49 million footprinted cardiomyocyte TFBS mapped onto the scan. Dynamic is split into its four directional classes because a gain and a loss could cancel in the pooled comparison. Tests run over all 88 matched class-by-family comparisons, not just the six families drawn. CTCF was reported as a lead earlier and is now retracted twice over: the original tile-level significance was pseudo-replication, and with the full scan and measured-accessibility matching the point estimate itself is zero (Dynamic -0.551 vs matched Stable -0.548, p = 0.98). The nearest thing to a real signal is MEF2 in Age Up, -0.53 against -0.84 in matched controls, p = 0.03 uncorrected, and in Exercise Down, p = 0.06 - the same two classes that are short of a close MEF2 on the previous slide. Neither survives correction on its own; they are worth quoting only because the proximity test points the same way.");

// ================================================================== 5 spacing
s = pres.addSlide({ masterName: "TITLE_FIGURE", sectionTitle: "Findings" });
s.addText("MEF2 sits right beside the AP-1 motif, and it matters", { placeholder: "title" });
s.addText("How often a partner sits there", { isTextBox: true, x: 0.8, y: 1.15, w: 5.4, h: 0.3, fontSize: 13.5, bold: true, color: "111111", margin: 0, objectName: "sp-t1" });
s.addText("every footprint of that family, important or not", { isTextBox: true, x: 0.8, y: 1.46, w: 5.4, h: 0.3, fontSize: 10.5, color: "595959", margin: 0, objectName: "sp-s1" });
FIG(s, "Analysis/slide_figures/panels/s7_occurrence.png", { x: 0.8, y: 1.85, w: 5.4, h: 3.5 });
s.addText("How much it matters there", { isTextBox: true, x: 6.9, y: 1.15, w: 5.6, h: 0.3, fontSize: 13.5, bold: true, color: "111111", margin: 0, objectName: "sp-t2" });
s.addText("against footprint-free tiles the same distance away", { isTextBox: true, x: 6.9, y: 1.46, w: 5.6, h: 0.3, fontSize: 10.5, color: "595959", margin: 0, objectName: "sp-s2" });
FIG(s, "Analysis/slide_figures/panels/s7_importance.png", { x: 6.9, y: 1.85, w: 5.6, h: 3.5 });
legendRow(s, [["MEF2", "C2551C"], ["KLF / SP", "9A9993"]], 0.8, 5.45);
s.addText("shaded band on the left: 11\u201325 bp, where a neighbour abuts the motif.   MEF2 is 2.67\u00d7 enriched there against 1.47\u00d7 for KLF/SP (FDR 0.004).",
  { isTextBox: true, x: 0.8, y: 5.78, w: 11.7, h: 0.4, fontSize: 9.5, color: "595959", margin: 0, objectName: "sp-key" });
takeaway(s, "MEF2 is 2.67\u00d7 enriched at 11\u201325 bp and is worth \u22120.9 log2 there; KLF/SP is the commonest neighbour and is worth nothing anywhere");
s.addNotes("Left panel counts every MEF2 footprint, important or not: it is occurrence, not importance. Middle panel is importance, measured against footprint-free tiles at the same distance, so the general fall-off with distance is removed. Right panel: 5.6% of Age Up sites carry a MEF2 footprint within 25 bp against 13.7% of their matched Stable controls (Fisher p 0.004, FDR 0.03); Exercise Down is similar but weaker. No 10-11 bp repeat survives correction, so no fixed protein-protein register.");

// ================================================================== 5b MEF2 by class
s = pres.addSlide({ masterName: "TITLE_FIGURE", sectionTitle: "Findings" });
s.addText("Age-gained sites are short of that MEF2", { placeholder: "title" });
classFigure(s, "s9_mef2_by_class", { x: 2.6, y: 1.25, w: 8.3, h: 4.9 });
takeaway(s, "Age Up 5.8% and Exercise Down 5.4% carry a close MEF2, against 8.5% of Stable sites \u2014 and the gap is sharper when only the MEF2 the model uses is counted");
s.addNotes("Matched tests are exact McNemar on the paired sites, raw tests are Fisher, BH across all 20. Surviving: Age Up on both definitions (FDR 0.021 leaned-on, 0.045 footprint-present) and Exercise Down on the leaned-on definition (FDR 0.026); Exercise Down footprint-present does not (FDR 0.12). Age Down and Exercise Up show nothing either way. An earlier quote of p = 0.003 and p = 0.0005 for these came from Fisher applied to the matched groups as if unpaired; McNemar is the right test and gives p = 0.0021 and p = 0.0013.");

// ================================================================== 5c presence vs rank
s = pres.addSlide({ masterName: "TITLE_FIGURE", sectionTitle: "Findings" });
s.addText("Is the partner missing, or just worth less?", { placeholder: "title" });
s.addText("Rank the 60 windows within each element, take the three most damaging, and ask whether any covers that family\u2019s footprint. A within-element question, so a class with less total dependence is not penalised for it.",
  { isTextBox: true, x: 0.8, y: 1.1, w: 11.9, h: 0.46, fontSize: 11.5, color: "595959", margin: 0, objectName: "or-intro" });
s.addText("MEF2", { isTextBox: true, x: 1.5, y: 1.72, w: 5.0, h: 0.28, fontSize: 13.5, bold: true, color: "111111", margin: 0, objectName: "or-t1" });
FIG(s, "Analysis/slide_figures/panels/s10_mef2_all.png", { x: 1.5, y: 2.06, w: 5.0, h: 3.3 });
s.addText("A second AP-1 footprint", { isTextBox: true, x: 7.0, y: 1.72, w: 5.0, h: 0.28, fontSize: 13.5, bold: true, color: "111111", margin: 0, objectName: "or-t2" });
FIG(s, "Analysis/slide_figures/panels/s10_ap1.png", { x: 7.0, y: 2.06, w: 5.0, h: 3.3 });
legendRow(s, [["does the element have one at all", "C2551C"], ["when it has one, is it ranked top", "4F4F4F"]], 1.5, 5.5);
s.addText("log2 odds ratio against all Stable sites with a 95% CI; dotted line = no difference, so \u22121 is half the odds and +1 is double.   "
  + "MEF2, does it have one: Age Up OR 0.45, Exercise Down 0.56, Dynamic 0.58, all FDR \u2264 0.005.   "
  + "AP-1, when it has one: Dynamic OR 1.45, FDR 0.036.",
  { isTextBox: true, x: 0.8, y: 5.88, w: 11.9, h: 0.5, fontSize: 9.5, color: "595959", margin: 0, objectName: "or-key" });
takeaway(s, "the MEF2 deficit is about having one, not about how it is ranked \u2014 and responsive elements lean more on a second AP-1 site when they have one");
s.addNotes("Presence vs rank is the point. MEF2, does the element have one: Dynamic OR 0.58, Age Up 0.45, Exercise Down 0.56, all surviving BH. Given a MEF2 is present, every class ranks it in the top three about equally often (0.87-2.03, nothing significant). AP-1 given present: Dynamic 1.45x, FDR 0.036, with the direction repeating in all four directional classes (1.43-1.62). CTCF is the reason to keep the matched analysis in the report even though it is off the slides: OR 3.0-6.5 here with FDR < 0.001, and exactly 1.0 once accessibility is matched.");

// ================================================================== 6 example
s = pres.addSlide({ masterName: "TITLE_FIGURE", sectionTitle: "Findings" });
s.addText("The same contrast at one pair of promoters", { placeholder: "title" });
FIG(s, "Analysis/example_track/example_scan_track.png", { x: 0.55, y: 1.15, w: 12.2, h: 5.2 });
takeaway(s, "the stable site has MEF2 abutting its AP-1 motif at 19 bp; the age-gained site has no MEF2 anywhere in the element");
s.addNotes("Both are promoters a few hundred bp from a TSS, matched on predicted accessibility (0.96 vs 1.12) and peak width (1359 vs 1329 bp). At the Sqstm1 promoter AP-1 is bound in all three conditions, the motif window is worth -1.00 and 81% of the loss sits within 85 bp. At the Psmd11 promoter AP-1 is absent in young and gained in both aged states, the motif is still the strongest window but only -0.58, and the loss is smeared over 114 bp against 83. Psmd11 rises with age, matching the AP-1 gain; the exercise direction does not match, so do not claim it.");

// ================================================================== 7 diffuse
s = pres.addSlide({ masterName: "TITLE_FIGURE", sectionTitle: "Findings" });
s.addText("Responsive sites spread their dependence out", { placeholder: "title" });
s.addText("Where the dependence sits across the element", { isTextBox: true, x: 0.55, y: 1.14, w: 7.6, h: 0.3, fontSize: 13.5, bold: true, color: "111111", margin: 0, objectName: "prof-t1" });
s.addText("mean loss at each of the 60 windows, smoothed over 3 \u2014 summed across windows, exactly the quantity on the right. Gold = the \u00b115 bp core.", { isTextBox: true, x: 0.55, y: 1.45, w: 7.6, h: 0.3, fontSize: 10.5, color: "595959", margin: 0, objectName: "prof-s1" });
FIG(s, "Analysis/slide_figures/panels/s12_profile_absolute.png", { x: 0.55, y: 1.76, w: 7.6, h: 1.95 });
s.addText("...and the same as a share of each site\u2019s own total", { isTextBox: true, x: 0.55, y: 3.78, w: 7.6, h: 0.3, fontSize: 13.5, bold: true, color: "111111", margin: 0, objectName: "prof-t2" });
s.addText("depth divided out, so only the shape is compared", { isTextBox: true, x: 0.55, y: 4.07, w: 7.6, h: 0.3, fontSize: 10.5, color: "595959", margin: 0, objectName: "prof-s2" });
FIG(s, "Analysis/slide_figures/panels/s12_profile_shape.png", { x: 0.55, y: 4.32, w: 7.6, h: 1.85 });
classFigure(s, "s3_dependence", { x: 8.35, y: 1.14, w: 4.4, h: 5.2 }, { gutter: 0.02, foot: 0.0, stack: true, legendX: 0.55, legendW: 12.2, legendY: 6.26, showKey: false });
takeaway(s, "Stable sites peak highest at the core; Age Up, Exercise Down and Dynamic carry more of their dependence out in the flanks \u2014 Age Down simply has less of everything");
s.addNotes("The profiles are means, not medians, so the area under a curve equals the per-site totals in the panels on the right; with medians the two panels are different functions and do not correspond. Note where the difference lives: at the single deepest window Stable and Age Up are almost identical (0.68 vs 0.63 mean), and the 18% gap in loss within 85 bp comes from the shoulders, not the peak. Spread is the standard deviation of window position when each of the 60 windows is weighted by its share of that site's total loss, measured about the centre of mass of the loss rather than about the motif. For scale: all the loss in one window gives 0 bp, confined to the central 90 bp gives 26, to the central 170 bp gives 49, and spread evenly over the whole 600 bp gives 173. Observed values run 29 to 216 with a median of 87, so a typical element draws on a region well over 100 bp wide. Caveat worth having ready: spread is partly confounded with weakness. Across all sites it correlates with the Gini concentration at rho -0.71, but with the NUMBER of windows that matter at only -0.22, and the broadest fifth of sites have FEWER windows that matter (5 vs 7), a shallower motif window (-0.44 vs -0.99) and less total loss (6.4 vs 9.3) than the tightest fifth. So a broad profile here means no stretch dominates, not that many stretches contribute. The honest sentence is that responsive elements have no dominant anchor -- their dependence is both weaker and more evenly scattered -- which is the same thing the near-versus-far panels say.");

// ================================================================== 8 titration
s = pres.addSlide({ masterName: "TITLE_FIGURE", sectionTitle: "Findings" });
s.addText("Position is a dial, and a moved motif never replaces the original", { placeholder: "title" });
s.addText("Where the added copy works best", { isTextBox: true, x: 0.55, y: 1.12, w: 3.1, h: 0.3, fontSize: 12.5, bold: true, color: "111111", margin: 0, objectName: "copy-title" });
FIG(s, "Analysis/slide_figures/panels/d_addcopy.png", { x: 0.55, y: 1.45, w: 3.1, h: 1.95 });
s.addText("...and where it lands when moved", { isTextBox: true, x: 0.55, y: 3.52, w: 3.1, h: 0.3, fontSize: 12.5, bold: true, color: "111111", margin: 0, objectName: "move-title" });
FIG(s, "Analysis/slide_figures/panels/d_relocate.png", { x: 0.55, y: 3.85, w: 3.1, h: 1.95 });
s.addText("both curves drawn once per class; colours as in the legend", { isTextBox: true, x: 0.55, y: 5.85, w: 3.2, h: 0.5, fontSize: 9.5, color: "595959", margin: 0, objectName: "curve-caption" });
classFigure(s, "s4_titration", { x: 4.3, y: 1.25, w: 8.45, h: 4.9 });
takeaway(s, "a second copy is worth most at \u00b120\u201340 bp and a moved motif recovers only a third of what destroying it cost");
s.addNotes("The 1-base and 7-base core results are not repeated here: they say the same thing as the motif shuffle on the earlier slide. One wrong base costs -0.45 of the -0.56 total, the dose slope tracks core quality read from the genome at rho 0.30, and orientation is irrelevant as expected for a near-palindrome.");

// ================================================================== 9 age-gained
s = pres.addSlide({ masterName: "TITLE_FIGURE", sectionTitle: "Findings" });
s.addText("Age-gained sites read AP-1 dosage more sharply", { placeholder: "title" });
s.addText("Age-gained sites respond more to AP-1 dosage, both ways", { isTextBox: true, x: 0.8, y: 1.15, w: 8.0, h: 0.32, fontSize: 13.5, bold: true, color: "111111", margin: 0, objectName: "ag-t1" });
s.addText("\u0394score definition: the top fifth of the change in TOBIAS binding score with age, against matched Unchanged sites", { isTextBox: true, x: 0.8, y: 1.47, w: 8.0, h: 0.32, fontSize: 10.5, color: "595959", margin: 0, objectName: "ag-s1" });
FIG(s, "Analysis/slide_figures/panels/s5_age_gained.png", { x: 0.8, y: 1.9, w: 8.0, h: 3.3 });
legendRow(s, [["Matched Unchanged", "B3B2AE"], ["Age-gained (\u0394score)", "C2551C"]], 0.8, 5.3);
s.addText([
  { text: "Mass action   ", options: { bold: true, fontSize: 17, color: "C2551C" } },
  { text: "if ageing raises AP-1 activity, occupancy rises first at the highest-affinity motifs \u2014 exactly the sites that gain. Losing binding is not the reverse: age-lost sites show no sequence signature at all.", options: { fontSize: 15 } },
], { isTextBox: true, x: 0.55, y: 5.5, w: 12.2, h: 1.0, fill: { color: "F2F2F2" }, margin: 12, objectName: "massaction" });

s.addNotes("The titration experiment was designed before this hypothesis existed, so this is a prediction rather than a fit. Caveat: the motif-sensitivity result does not survive matching on measured accessibility; the core-quality result does.");

// ================================================================== 10 MPRA
s = pres.addSlide({ masterName: "TITLE_FIGURE", sectionTitle: "Findings" });
s.addText("An independent reporter assay agrees", { placeholder: "title" });
s.addText("Same answer: the motif is a small part", { isTextBox: true, x: 0.8, y: 2.55, w: 5.4, h: 0.3, fontSize: 13, bold: true, color: "111111", margin: 0, objectName: "mp-t1" });
FIG(s, "Analysis/slide_figures/panels/s6_share.png", { x: 0.8, y: 2.9, w: 5.4, h: 2.9 });
s.addText("But they disagree on which elements depend more", { isTextBox: true, x: 6.9, y: 2.55, w: 5.6, h: 0.3, fontSize: 13, bold: true, color: "111111", margin: 0, objectName: "mp-t2" });
FIG(s, "Analysis/slide_figures/panels/s6_direction.png", { x: 6.9, y: 2.9, w: 5.6, h: 2.9 });
s.addText("below 1\u00d7: responsive elements lose less than Stable ones   \u00b7   above 1\u00d7: they lose more   \u00b7   both p < 0.001",
  { isTextBox: true, x: 0.8, y: 5.85, w: 11.7, h: 0.3, fontSize: 9.5, color: "595959", margin: 0, objectName: "mp-key" });
statRow(s, [
  { v: "7.5%", l: "AlphaGenome: the motif's share of local dependence", c: "1F6FB2" },
  { v: "9.3%", l: "fibroblast MPRA: the same quantity, measured in the lab", c: "C2551C" },
], 1.2, 1.15);
takeaway(s, "two assays, different cell types and readouts, both put the AP-1 motif under a tenth of its 170 bp surroundings");
s.addNotes("The comparable quantity is a ratio within each element, so assay scale cancels. The class direction disagrees: in the MPRA responsive elements depend MORE on their sequence, which separates intrinsic potential from what sets accessibility in the genome.");

// ================================================================== 11 so what
pres.addSection({ title: "Conclusion" });
s = pres.addSlide({ masterName: "TITLE_BODY", sectionTitle: "Conclusion" });
s.addText("What this changes", { placeholder: "title" });
s.addText([
  { text: "A responsive AP-1 site is a marginal one, not a differently-built one", options: { bold: true, fontSize: 19, breakLine: true } },
  { text: "Weaker AP-1 motif, fewer AP-1 motifs around it, and its dependence spread thinly across the element instead of anchored on the motif. The same partners matter, in the same order, in every class.", options: { breakLine: true, paraSpaceAfter: 15 } },
  { text: "Two things survive matching, and only two", options: { bold: true, fontSize: 19, breakLine: true } },
  { text: "Scale: the 11 bp motif and the whole peak are worth the same to both classes, but replacing the 170 bp element costs a Dynamic site 0.24 less — what you expect when the dependence is spread wider. Partner: age-gained sites lack a MEF2 abutting the motif, 5.6% against 13.7% of matched controls.", options: { breakLine: true, paraSpaceAfter: 15 } },
  { text: "That is why it fluctuates", options: { bold: true, fontSize: 19, breakLine: true } },
  { text: "A weak motif with no close MEF2 to hold it sits near the occupancy threshold, so its binding tracks how much AP-1 protein is available — which is what ageing changes.", options: { breakLine: true, paraSpaceAfter: 15 } },
  { text: "Gain and loss are different events", options: { bold: true, fontSize: 19, breakLine: true } },
  { text: "Age-gained sites have better cores and read AP-1 dosage more sharply. Age-lost sites show no sequence signature of their own.", options: {} },
], { placeholder: "body", w: 12.2, h: 5.1 });
s.addNotes("The mutation-response profile is a null and a well-controlled one: given the same chromatin, the two classes lose the same amount to the same kind of edit. What is not null is the scale at which that is true and the one missing partner.");

// ================================================================== 12 next
s = pres.addSlide({ masterName: "TITLE_BODY", sectionTitle: "Conclusion" });
s.addText("Where this goes next", { placeholder: "title" });
s.addText([
  { text: "Fine-tuned condition-aware head", options: { bold: true, fontSize: 19, breakLine: true } },
  { text: "A frozen-trunk probe on our own cardiomyocyte and fibroblast ATAC already reaches profile r = 0.65 against 0.59 for the built-in heart head, training 0.01% of the model. Next: 1 bp insertion targets, then LoRA.", options: { breakLine: true, paraSpaceAfter: 14 } },
  { text: "Does AP-1 activity actually rise with age?", options: { bold: true, fontSize: 19, breakLine: true } },
  { text: "The mass-action reading rests on it. One check against the RNA.", options: { breakLine: true, paraSpaceAfter: 14 } },
  { text: "Fibroblast scan", options: { bold: true, fontSize: 19, breakLine: true } },
  { text: "If MEF2 also tops the partner ranking in fibroblasts, the ranking is a model artefact rather than cardiac biology. This is the sharpest available control.", options: {} },
], { placeholder: "body", w: 12.2, h: 4.6 });
s.addNotes("Also open: whether the MEF2 deficit at Dynamic sites is detection or genuine differential occupancy of equally-present motifs.");

(async () => {
  const out = F("Analysis/AP1_AlphaGenome_results.pptx");
  await pres.writeFile({ fileName: out });
  await applyTheme(out, THEME);
  console.log("wrote", out);
})();
