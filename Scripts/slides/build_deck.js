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
  const gut = opts.gutter ?? 0.72, head = 0.82, foot = opts.foot ?? 0.92, gap = 0.08;
  const cw = (box.w - gut) / nCol, ch = (box.h - head - foot) / nRow;
  man.cols.forEach((c, ci) => {
    const x = box.x + gut + ci * cw;
    slide.addText(c.title, { isTextBox: true, x, y: box.y, w: cw - gap, h: 0.48, fontSize: 12.5, bold: true, color: "111111", valign: "top", margin: 0, objectName: `${name}-title-${ci}` });
    slide.addText(c.subtitle, { isTextBox: true, x, y: box.y + 0.47, w: cw - gap, h: 0.34, fontSize: 9.5, color: "595959", valign: "top", margin: 0, objectName: `${name}-sub-${ci}` });
    man.panels[ci].forEach((f, ri) => {
      FIG(slide, `Analysis/slide_figures/panels/${f}`, { x, y: box.y + head + ri * ch, w: cw - gap, h: ch - gap });
    });
  });
  man.rows.forEach((r, ri) => {
    slide.addText(r, { isTextBox: true, x: box.x - 0.52, y: box.y + head + ri * ch + ch / 2 - 0.75, w: 1.5, h: 1.5,
      fontSize: 10, color: "111111", align: "center", valign: "middle", rotate: 270, margin: 0, objectName: `${name}-row-${ri}` });
  });
  let lx = box.x + gut;
  const ly = box.y + box.h - foot + 0.10;
  man.legend.forEach(([lab, col], i) => {
    slide.addShape(pres.ShapeType.ellipse, { x: lx, y: ly + 0.045, w: 0.12, h: 0.12, fill: { color: col.replace("#", "") }, line: { color: col.replace("#", ""), width: 0 }, objectName: `${name}-dot-${i}` });
    const w = 0.22 + lab.length * 0.072;
    slide.addText(lab, { isTextBox: true, x: lx + 0.17, y: ly, w, h: 0.22, fontSize: 10.5, color: "111111", margin: 0, objectName: `${name}-leg-${i}` });
    lx += 0.17 + w + 0.12;
  });
  slide.addText(man.key, { isTextBox: true, x: box.x, y: ly + 0.30, w: box.w, h: 0.55, fontSize: 9, color: "595959", margin: 0, objectName: `${name}-key` });
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
pres.addSection({ title: "Findings" });
s = pres.addSlide({ masterName: "TITLE_FIGURE", sectionTitle: "Findings" });
s.addText("Most of the difference is chromatin, not sequence", { placeholder: "title" });
classFigure(s, "s2_matched_null", { x: 1.05, y: 1.1, w: 11.7, h: 5.1 });
takeaway(s, "every class behaves the same way: the motif and the whole peak close up under matching, the 170 bp element does not \u2014 except Age Down, which closes up everywhere");
s.addNotes("Top row is the uncontrolled comparison: if the accessibility difference is itself biology, it lives there. Bottom row matches each class to its own Stable controls on predicted and measured accessibility, peak width and TSS distance. The 300 bp control is the thing to check: Age Up differs slightly even there (-0.152 vs -0.111 matched), but that offset is 0.04 against a 0.26 effect at 170 bp and runs the opposite way, so it does not explain the result.");

// ================================================================== 4 partners
s = pres.addSlide({ masterName: "TITLE_FIGURE", sectionTitle: "Findings" });
s.addText("Both classes lean on the same partners", { placeholder: "title" });
FIG(s, "Analysis/slide_figures/d_tile.png", { x: 0.55, y: 1.05, w: 7.9, h: 1.1 });
s.addText("All sites, nothing controlled", { isTextBox: true, x: 0.8, y: 2.25, w: 5.4, h: 0.3, fontSize: 13, bold: true, color: "111111", margin: 0, objectName: "pa-t1" });
FIG(s, "Analysis/slide_figures/panels/s8_partners_all.png", { x: 0.8, y: 2.6, w: 5.4, h: 3.0 });
s.addText("Each class against its own matched Stable controls", { isTextBox: true, x: 6.9, y: 2.25, w: 5.6, h: 0.3, fontSize: 13, bold: true, color: "111111", margin: 0, objectName: "pa-t2" });
FIG(s, "Analysis/slide_figures/panels/s8_partners_matched.png", { x: 6.9, y: 2.6, w: 5.6, h: 3.0 });
legendRow(s, [["Stable", "B3B2AE"], ["Age Up", "C2551C"], ["Age Down", "1F6FB2"], ["Exercise Up", "2E7D32"], ["Exercise Down", "7B3FA0"]], 0.8, 5.68);
s.addText("six families shown, tests run over all 40 class-by-family comparisons; smallest FDR 0.04, and only before matching",
  { isTextBox: true, x: 0.8, y: 5.98, w: 11.7, h: 0.3, fontSize: 9.5, color: "595959", margin: 0, objectName: "pa-key" });
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
classFigure(s, "s9_mef2_by_class", { x: 2.2, y: 1.1, w: 9.4, h: 5.1 });
takeaway(s, "Age Up 3.2% against 10.5% of its matched controls, Exercise Down 1.2% against 9.8% \u2014 and the deficit is sharper when only the MEF2 the model uses is counted");
s.addNotes("Matched tests are exact McNemar on the paired sites, raw tests are Fisher, BH across all 20. Surviving: Age Up on both definitions (FDR 0.021 leaned-on, 0.045 footprint-present) and Exercise Down on the leaned-on definition (FDR 0.026); Exercise Down footprint-present does not (FDR 0.12). Age Down and Exercise Up show nothing either way. An earlier quote of p = 0.003 and p = 0.0005 for these came from Fisher applied to the matched groups as if unpaired; McNemar is the right test and gives p = 0.0021 and p = 0.0013.");

// ================================================================== 6 example
s = pres.addSlide({ masterName: "TITLE_FIGURE", sectionTitle: "Findings" });
s.addText("The same contrast at one pair of promoters", { placeholder: "title" });
FIG(s, "Analysis/example_track/example_scan_track.png", { x: 0.55, y: 1.15, w: 12.2, h: 5.2 });
takeaway(s, "the stable site has MEF2 abutting its AP-1 motif at 19 bp; the age-gained site has no MEF2 anywhere in the element");
s.addNotes("Both are promoters a few hundred bp from a TSS, matched on predicted accessibility (0.96 vs 1.12) and peak width (1359 vs 1329 bp). At the Sqstm1 promoter AP-1 is bound in all three conditions, the motif window is worth -1.00 and 81% of the loss sits within 85 bp. At the Psmd11 promoter AP-1 is absent in young and gained in both aged states, the motif is still the strongest window but only -0.58, and the loss is smeared over 114 bp against 83. Psmd11 rises with age, matching the AP-1 gain; the exercise direction does not match, so do not claim it.");

// ================================================================== 7 diffuse
s = pres.addSlide({ masterName: "TITLE_FIGURE", sectionTitle: "Findings" });
s.addText("Responsive sites spread their dependence out", { placeholder: "title" });
FIG(s, "Analysis/slide_figures/panels/d_profile.png", { x: 0.55, y: 1.5, w: 2.5, h: 2.1 });
s.addText("one site\u2019s sliding-window profile: each bar is a 10 bp window replaced, shaded band = within \u00b185 bp of the motif",
  { isTextBox: true, x: 0.55, y: 3.65, w: 2.5, h: 0.9, fontSize: 9.5, color: "595959", margin: 0, objectName: "profile-caption" });
classFigure(s, "s3_dependence", { x: 3.75, y: 1.1, w: 9.0, h: 5.1 }, { gutter: 0.62 });
takeaway(s, "the same shift in every class and on all three measures \u2014 strongest in Age Up and Exercise Down, absent only in Age Down");
s.addNotes("10 of 16 scan features separate Dynamic from matched Stable. Share of dependence within 85 bp 0.63 vs 0.74, spread 94 vs 86 bp, largest single-base effect 0.74 vs 0.90. Age Up is the most diffuse class of all.");

// ================================================================== 8 titration
s = pres.addSlide({ masterName: "TITLE_FIGURE", sectionTitle: "Findings" });
s.addText("Motif strength is a switch; position is a dial", { placeholder: "title" });
s.addText("The dose curve, class by class", { isTextBox: true, x: 0.55, y: 1.12, w: 2.6, h: 0.3, fontSize: 12.5, bold: true, color: "111111", margin: 0, objectName: "dose-title" });
FIG(s, "Analysis/slide_figures/panels/d_dose.png", { x: 0.55, y: 1.45, w: 2.5, h: 1.85 });
s.addText("Where a second copy helps", { isTextBox: true, x: 0.55, y: 3.42, w: 2.6, h: 0.3, fontSize: 12.5, bold: true, color: "111111", margin: 0, objectName: "copy-title" });
FIG(s, "Analysis/slide_figures/panels/d_addcopy.png", { x: 0.55, y: 3.75, w: 2.5, h: 1.85 });
s.addText("both curves drawn once per class; colours as in the legend", { isTextBox: true, x: 0.55, y: 5.65, w: 2.7, h: 0.5, fontSize: 9.5, color: "595959", margin: 0, objectName: "curve-caption" });
classFigure(s, "s4_titration", { x: 3.75, y: 1.1, w: 9.0, h: 5.1 }, { gutter: 0.62 });
takeaway(s, "one wrong base costs \u22120.45 of the \u22120.56 total and a moved motif recovers only a third \u2014 but every class titrates the same way");
s.addNotes("The dose-response slope tracks core quality measured independently from the genome sequence, rho 0.30. Orientation is irrelevant, as expected for a near-palindrome.");

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
