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
  { text: "The TOBIAS footprint is absent in young and called in aged. 396 Age Up sites, 179 Age Down. Stable = called in all three conditions.", options: { breakLine: true, paraSpaceAfter: 14 } },
  { text: "Δscore (graded)", options: { bold: true, fontSize: 19, breakLine: true } },
  { text: "The top fifth of the change in binding score, whether or not the call ever crossed the line. 736 sites each way.", options: { breakLine: true, paraSpaceAfter: 14 } },
  { text: "Only 22% of Δscore-Up sites ever cross the threshold, and 40% of threshold-Stable sites show a large score change.", options: { color: "595959" } },
], { isTextBox: true, x: 0.55, y: 1.35, w: 6.3, h: 4.0, fontSize: 17, color: "111111", margin: 0 });
statRow(s, [
  { v: "0.57", l: "Δscore vs measured snATAC change (age)", c: "C2551C" },
  { v: "0.31", l: "threshold class vs measured change", c: "1F6FB2" },
], 5.55, 1.35);
s.addText("Every comparison below matches each responsive site to a Stable site on predicted accessibility, peak width, TSS distance and measured accessibility.",
  { isTextBox: true, x: 7.1, y: 1.5, w: 5.65, h: 1.5, fontSize: 15, color: "111111", fill: { color: "F2F2F2" }, margin: 12, objectName: "matching-note" });
s.addText("The Δscore definition tracks the measured data about twice as well.", { isTextBox: true, x: 7.1, y: 3.25, w: 5.65, h: 0.9, fontSize: 15, italic: true, color: "595959", margin: 0 });
s.addNotes("The two definitions pick out largely different sites. Report every result labelled by which definition it comes from.");

// ================================================================== 3 method schematic
s = pres.addSlide({ masterName: "TITLE_FIGURE", sectionTitle: "Design" });
s.addText("What we did to every site", { placeholder: "title" });
FIG(s, "Analysis/slide_figures/s1_method_schematic.png", BOX);
takeaway(s, "3,678 sites \u00d7 364 edits each, every edit read as the change in predicted accessibility over the surrounding 170 bp");
s.addNotes("The 131 kb window is what AlphaGenome takes as input; we only ever edit the middle 600 bp. Each tile is 3 replicate random replacements, averaged. Effects are log2(edited / native).");

// ================================================================== 4 the null
pres.addSection({ title: "Findings" });
s = pres.addSlide({ masterName: "TITLE_FIGURE", sectionTitle: "Findings" });
s.addText("Most of the difference is chromatin, not sequence", { placeholder: "title" });
FIG(s, "Analysis/slide_figures/s2_matched_null.png", BOX);
takeaway(s, "the motif and the whole peak behave the same once context is matched; a residual 170 bp difference survives, and the next slide shows what it is");
s.addNotes("Raw, Dynamic sites lose less to every wide edit. After matching on predicted accessibility, peak width and TSS distance, nothing survives. A classifier on the whole response profile is at chance, AUC 0.49.");

// ================================================================== 4 partners
s = pres.addSlide({ masterName: "TITLE_FIGURE", sectionTitle: "Findings" });
s.addText("Both classes lean on the same partners", { placeholder: "title" });
FIG(s, "Analysis/slide_figures/s8_partners.png", { x: 0.55, y: 1.2, w: 9.0, h: 5.0 });
s.addText([
  { text: "2.49 million", options: { bold: true, fontSize: 22, color: "C2551C", breakLine: true } },
  { text: "footprinted cardiomyocyte TFBS mapped onto the scan", options: { fontSize: 14, color: "595959", breakLine: true, paraSpaceAfter: 14 } },
  { text: "The ranking is the validation: MEF2, GATA, NFI and TEAD come out on top, and the model was never told any of them.", options: { fontSize: 14.5, breakLine: true, paraSpaceAfter: 12 } },
  { text: "The class split is the question, and the answer is no. Not one family matters more to a Dynamic site than to its matched Stable partner.", options: { fontSize: 14.5, breakLine: true, paraSpaceAfter: 12 } },
  { text: "Whatever makes a site responsive, it is not a different cast of neighbours.", options: { fontSize: 14.5, color: "595959", italic: true } },
], { isTextBox: true, x: 9.75, y: 1.5, w: 3.0, h: 4.4, margin: 0 });
takeaway(s, "the cast of partners is identical in both classes — what differs is how the dependence is spread across them");
s.addNotes("Importance is the loss from replacing a 10 bp tile over that factor's footprint, measured against the site's own footprint-free tiles, holding under five adjustments. The class comparison is at site level inside the 442 matched pairs; nothing survives BH correction. AP-1's own footprints are the one uncorrected p < 0.05, and Dynamic sites lean slightly MORE on them.");

// ================================================================== 5 spacing
s = pres.addSlide({ masterName: "TITLE_FIGURE", sectionTitle: "Findings" });
s.addText("The partner that matters is the one that sits closest", { placeholder: "title" });
FIG(s, "Analysis/slide_figures/s7_spacing.png", BOX);
takeaway(s, "MEF2 is 2.67\u00d7 enriched immediately beside the motif (FDR 0.004); KLF/SP, the commonest factor, is not");
s.addNotes("MEF2 is enriched 2.67x at 11-25 bp against 1.47x for KLF/SP, modal spacing 14 bp between centres. No 10-11 bp repeat survives correction, so no fixed protein-protein register.");

// ================================================================== 6 example
s = pres.addSlide({ masterName: "TITLE_FIGURE", sectionTitle: "Findings" });
s.addText("One locus, two architectures", { placeholder: "title" });
FIG(s, "Analysis/example_track/example_scan_track.png", { x: 0.55, y: 1.15, w: 12.2, h: 5.2 });
takeaway(s, "at the age-gained site every cardiac cofactor is bound in all three conditions — only AP-1 switches");
s.addNotes("Top: a Stable site, 85% of its dependence on the AP-1 motif. Bottom: an age-gained site where the AP-1 motif is worth -0.18 and the strongest other window -1.48. Tnfaip8 rises with age, matching the AP-1 gain.");

// ================================================================== 7 diffuse
s = pres.addSlide({ masterName: "TITLE_FIGURE", sectionTitle: "Findings" });
s.addText("Responsive sites spread their dependence out", { placeholder: "title" });
FIG(s, "Analysis/slide_figures/s3_dependence.png", BOX);
takeaway(s, "this is the part that survives matching: a responsive site has no single anchor holding it open");
s.addNotes("10 of 16 scan features separate Dynamic from matched Stable. Share of dependence within 85 bp 0.63 vs 0.74, spread 94 vs 86 bp, largest single-base effect 0.74 vs 0.90. Age Up is the most diffuse class of all.");

// ================================================================== 8 titration
s = pres.addSlide({ masterName: "TITLE_FIGURE", sectionTitle: "Findings" });
s.addText("Motif strength is a switch; position is a dial", { placeholder: "title" });
FIG(s, "Analysis/slide_figures/s4_titration.png", { x: 0.55, y: 1.2, w: 12.2, h: 3.75 });
statRow(s, [
  { v: "−0.45", l: "cost of changing one core base" },
  { v: "−0.56", l: "cost of changing all seven" },
  { v: "−0.29", l: "destroyed, re-inserted nearby" },
  { v: "±20–40 bp", l: "where a motif works best", c: "C2551C" },
], 5.15, 1.15);
takeaway(s, "a motif moved elsewhere recovers only a third of what destroying it cost \u2014 position matters, not just presence");
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
  { text: "Weaker AP-1 motif, fewer AP-1 motifs around it, sparser and less MEF2-rich neighbourhood, dependence spread thinly across the element. Nothing about how it responds to mutation differs.", options: { breakLine: true, paraSpaceAfter: 16 } },
  { text: "That is why it fluctuates", options: { bold: true, fontSize: 19, breakLine: true } },
  { text: "A weak motif in a poorly-supported enhancer sits near the occupancy threshold, so its binding tracks how much AP-1 protein is available — which is what ageing changes.", options: { breakLine: true, paraSpaceAfter: 16 } },
  { text: "Gain and loss are different events", options: { bold: true, fontSize: 19, breakLine: true } },
  { text: "Age-gained sites have better cores and read dosage more sharply. Age-lost sites are MEF2-rich, high-scoring sites drifting down, with no sequence signature of their own.", options: {} },
], { placeholder: "body", w: 12.2, h: 4.6 });
s.addNotes("The negative result is well controlled and is itself the point: responsiveness is not encoded in local sequence architecture.");

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
