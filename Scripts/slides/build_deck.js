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
FIG(s, "Analysis/slide_figures/d_block.png", BAND);
FIG(s, "Analysis/slide_figures/s2_matched_null.png", BOXB);
takeaway(s, "the motif and the whole peak behave the same once context is matched; a residual 170 bp difference survives, and the next slide shows what it is");
s.addNotes("Raw, Dynamic sites lose less to every wide edit. Matching now holds measured snATAC accessibility constant as well as the predicted value, peak width and TSS distance; that removes most of the gap but not all of it. What survives is confined to the 170 bp scale: the motif and the whole peak are worth the same to both classes.");

// ================================================================== 4 partners
s = pres.addSlide({ masterName: "TITLE_FIGURE", sectionTitle: "Findings" });
s.addText("Both classes lean on the same partners", { placeholder: "title" });
FIG(s, "Analysis/slide_figures/d_tile.png", BAND);
FIG(s, "Analysis/slide_figures/s8_partners.png", { x: 0.55, y: 2.15, w: 12.2, h: 4.15 });
takeaway(s, "the ranking is the validation \u2014 MEF2, GATA, NFI, TEAD, none of it told to the model; the class split is the question, and no class leans on a different partner");
s.addNotes("2.49 million footprinted cardiomyocyte TFBS mapped onto the scan. Dynamic is split into its four directional classes because a gain and a loss could cancel in the pooled comparison. Tests run over all 88 matched class-by-family comparisons, not just the six families drawn. CTCF was reported as a lead earlier and is now retracted twice over: the original tile-level significance was pseudo-replication, and with the full scan and measured-accessibility matching the point estimate itself is zero (Dynamic -0.551 vs matched Stable -0.548, p = 0.98). The nearest thing to a real signal is MEF2 in Age Up, -0.53 against -0.84 in matched controls, p = 0.03 uncorrected, and in Exercise Down, p = 0.06 - the same two classes that are short of a close MEF2 on the previous slide. Neither survives correction on its own; they are worth quoting only because the proximity test points the same way.");

// ================================================================== 5 spacing
s = pres.addSlide({ masterName: "TITLE_FIGURE", sectionTitle: "Findings" });
s.addText("Age-gained sites are missing their closest MEF2", { placeholder: "title" });
FIG(s, "Analysis/slide_figures/s7_spacing.png", { x: 0.55, y: 1.25, w: 12.2, h: 4.9 });
takeaway(s, "MEF2 sits 2.67\u00d7 more often right beside the motif than chance, and it is the age-gained sites that lack one (FDR 0.03)");
s.addNotes("Left panel counts every MEF2 footprint, important or not: it is occurrence, not importance. Middle panel is importance, measured against footprint-free tiles at the same distance, so the general fall-off with distance is removed. Right panel: 5.6% of Age Up sites carry a MEF2 footprint within 25 bp against 13.7% of their matched Stable controls (Fisher p 0.004, FDR 0.03); Exercise Down is similar but weaker. No 10-11 bp repeat survives correction, so no fixed protein-protein register.");

// ================================================================== 6 example
s = pres.addSlide({ masterName: "TITLE_FIGURE", sectionTitle: "Findings" });
s.addText("The same contrast at one pair of promoters", { placeholder: "title" });
FIG(s, "Analysis/example_track/example_scan_track.png", { x: 0.55, y: 1.15, w: 12.2, h: 5.2 });
takeaway(s, "the stable site has MEF2 abutting its AP-1 motif at 19 bp; the age-gained site has no MEF2 anywhere in the element");
s.addNotes("Both are promoters a few hundred bp from a TSS, matched on predicted accessibility (0.96 vs 1.12) and peak width (1359 vs 1329 bp). At the Sqstm1 promoter AP-1 is bound in all three conditions, the motif window is worth -1.00 and 81% of the loss sits within 85 bp. At the Psmd11 promoter AP-1 is absent in young and gained in both aged states, the motif is still the strongest window but only -0.58, and the loss is smeared over 114 bp against 83. Psmd11 rises with age, matching the AP-1 gain; the exercise direction does not match, so do not claim it.");

// ================================================================== 7 diffuse
s = pres.addSlide({ masterName: "TITLE_FIGURE", sectionTitle: "Findings" });
s.addText("Responsive sites spread their dependence out", { placeholder: "title" });
FIG(s, "Analysis/slide_figures/d_base.png", BAND);
FIG(s, "Analysis/slide_figures/s3_dependence.png", BOXB);
takeaway(s, "this is the part that survives matching: a responsive site has no single anchor holding it open");
s.addNotes("10 of 16 scan features separate Dynamic from matched Stable. Share of dependence within 85 bp 0.63 vs 0.74, spread 94 vs 86 bp, largest single-base effect 0.74 vs 0.90. Age Up is the most diffuse class of all.");

// ================================================================== 8 titration
s = pres.addSlide({ masterName: "TITLE_FIGURE", sectionTitle: "Findings" });
s.addText("Motif strength is a switch; position is a dial", { placeholder: "title" });
FIG(s, "Analysis/slide_figures/d_titr.png", BAND);
FIG(s, "Analysis/slide_figures/s4_titration.png", { x: 0.55, y: 2.15, w: 12.2, h: 4.15 });
takeaway(s, "one wrong base costs \u22120.45 of the \u22120.56 total and a moved motif recovers only a third \u2014 but every class titrates the same way");
s.addNotes("The dose-response slope tracks core quality measured independently from the genome sequence, rho 0.30. Orientation is irrelevant, as expected for a near-palindrome.");

// ================================================================== 9 age-gained
s = pres.addSlide({ masterName: "TITLE_FIGURE", sectionTitle: "Findings" });
s.addText("Age-gained sites read AP-1 dosage more sharply", { placeholder: "title" });
FIG(s, "Analysis/slide_figures/s5_age_gained.png", { x: 0.55, y: 1.2, w: 12.2, h: 4.15 });
s.addText([
  { text: "Mass action   ", options: { bold: true, fontSize: 17, color: "C2551C" } },
  { text: "if ageing raises AP-1 activity, occupancy rises first at the highest-affinity motifs \u2014 exactly the sites that gain. Losing binding is not the reverse: age-lost sites show no sequence signature at all.", options: { fontSize: 15 } },
], { isTextBox: true, x: 0.55, y: 5.5, w: 12.2, h: 1.0, fill: { color: "F2F2F2" }, margin: 12, objectName: "massaction" });

s.addNotes("The titration experiment was designed before this hypothesis existed, so this is a prediction rather than a fit. Caveat: the motif-sensitivity result does not survive matching on measured accessibility; the core-quality result does.");

// ================================================================== 10 MPRA
s = pres.addSlide({ masterName: "TITLE_FIGURE", sectionTitle: "Findings" });
s.addText("An independent reporter assay agrees", { placeholder: "title" });
FIG(s, "Analysis/slide_figures/s6_mpra.png", { x: 0.55, y: 2.45, w: 12.2, h: 3.85 });
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
