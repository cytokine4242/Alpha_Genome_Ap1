// Lab-meeting deck for the fibroblast AlphaGenome results.
// One result per slide, one plot per result, and the biological takeaway written under it; the statistics and the
// caveats live in the speaker notes so the slide itself stays readable from the back of the room.
//
//   npm install pptxgenjs image-size            # once, in this directory
//   python3 make_fb_simple_figures.py           # the plots
//   NODE_PATH="$PWD/node_modules" node build_fb_deck.js
//
// Writes Analysis/AP1_fibroblast_lab_meeting.pptx.

const pptxgen = require("pptxgenjs");
const path = require("path");
const fs = require("fs");
const { applyTheme } = require("/Users/z5162865/.claude/skills/synced/2361eaa7-acd4-463c-8e26-dee27087b623_1ef5f431-f760-4047-8679-39ac25a5dd75/pptx/scripts/apply_theme.js");

const ROOT = "/Users/z5162865/WongLab/AlphaGenome_Ap1";
const F = (p) => path.join(ROOT, p);
const OUTFILE = F("Analysis/AP1_fibroblast_lab_meeting.pptx");

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
pres.layout = "LAYOUT_WIDE";                       // 13.3 x 7.5
pres.theme = { headFontFace: THEME.headFontFace, bodyFontFace: THEME.bodyFontFace };
pres.author = "Jack";
pres.title = "AP-1 sites in cardiac fibroblasts: what AlphaGenome says";
const C = pres.SchemeColor;

// ------------------------------------------------------------------ layouts
pres.defineSlideMaster({
  title: "TITLE", background: { color: "FFFFFF" },
  objects: [
    { placeholder: { options: { name: "title", type: "title", x: 0.9, y: 2.3, w: 11.5, h: 1.6, fontSize: 38, bold: true, color: C.text1, align: "left", valign: "bottom" }, text: " " } },
    { placeholder: { options: { name: "body", type: "body", x: 0.9, y: 4.1, w: 11.5, h: 1.6, fontSize: 17, color: C.text2, align: "left", valign: "top" }, text: " " } },
  ],
});

// a result slide: title, the plot, and one sentence of biology under it
pres.defineSlideMaster({
  title: "RESULT", background: { color: "FFFFFF" },
  objects: [
    { placeholder: { options: { name: "title", type: "title", x: 0.55, y: 0.3, w: 12.2, h: 0.8, fontSize: 28, bold: true, color: C.text1, align: "left", valign: "middle" }, text: " " } },
  ],
  slideNumber: { x: 12.75, y: 7.05, fontSize: 10, color: "A6A6A6" },
});

pres.defineSlideMaster({
  title: "TITLE_BODY", background: { color: "FFFFFF" },
  objects: [
    { placeholder: { options: { name: "title", type: "title", x: 0.55, y: 0.3, w: 12.2, h: 0.8, fontSize: 28, bold: true, color: C.text1, align: "left", valign: "middle" }, text: " " } },
    { placeholder: { options: { name: "body", type: "body", x: 0.55, y: 1.35, w: 12.2, h: 5.5, fontSize: 17, color: C.text1, align: "left", valign: "top" }, text: " " } },
  ],
  slideNumber: { x: 12.75, y: 7.05, fontSize: 10, color: "A6A6A6" },
});

// ------------------------------------------------------------------ helpers
const sizeOf = require("image-size").default || require("image-size");
// fit an image inside a box, preserving its real aspect ratio - stretching is the commonest slide defect
function FIG(slide, file, box) {
  const { width, height } = sizeOf(fs.readFileSync(F(file)));
  const r = width / height;
  let w = box.w, h = w / r;
  if (h > box.h) { h = box.h; w = h * r; }
  slide.addImage({ path: F(file), x: box.x + (box.w - w) / 2, y: box.y + (box.h - h) / 2, w, h });
}

// title / figure / takeaway, with the takeaway in an accent rule so the eye lands on it last
function resultSlide({ title, figure, method, takeaway, caveat, notes }) {
  const s = pres.addSlide({ masterName: "RESULT" });
  s.addText(title, { placeholder: "title" });
  s.addText(method, { isTextBox: true, x: 0.55, y: 1.06, w: 12.2, h: 0.34, fontSize: 12.5, italic: true, color: "595959", margin: 0 });
  FIG(s, figure, { x: 1.1, y: 1.46, w: 11.1, h: caveat ? 3.9 : 4.22 });
  const y0 = caveat ? 5.52 : 5.80;
  s.addShape(pres.ShapeType.rect, { x: 0.55, y: y0, w: 0.055, h: caveat ? 0.8 : 0.95, fill: { color: "C2551C" }, line: { width: 0 } });
  s.addText(takeaway, { isTextBox: true, x: 0.82, y: y0, w: 11.9, h: caveat ? 0.8 : 0.95, fontSize: 17, color: "111111", valign: "middle", margin: 0 });
  if (caveat) s.addText(caveat, { isTextBox: true, x: 0.82, y: 6.40, w: 11.9, h: 0.66, fontSize: 13, color: "595959", valign: "top", margin: 0 });
  s.addNotes(notes);
  return s;
}

// pptxgenjs runs array items together unless each one asks for its own paragraph
const bullets = (items) => items.map((t) => (typeof t === "string"
  ? { text: t, options: { bullet: { code: "2022" }, breakLine: true, paraSpaceAfter: 9, fontSize: 17 } }
  : { text: t.text, options: { bullet: t.sub ? { code: "25AA" } : { code: "2022" }, indentLevel: t.sub ? 1 : 0,
      breakLine: true, paraSpaceBefore: t.sub ? 0 : 11, paraSpaceAfter: t.sub ? 4 : 3,
      fontSize: t.sub ? 14.5 : 17, color: t.sub ? "595959" : "111111", bold: !!t.bold } }));

// ================================================================== 1. title
{
  const s = pres.addSlide({ masterName: "TITLE" });
  s.addText("AP-1 sites in cardiac fibroblasts", { placeholder: "title" });
  s.addText([
    { text: "What AlphaGenome says separates Stable from Dynamic sites", options: { fontSize: 19, color: "111111", breakLine: true, paraSpaceAfter: 12 } },
    { text: "Fine-tuned 1 bp ATAC head (held-out Spearman 0.68, against 0.27 for the built-in head)  ·  5,018 fibroblast AP-1 sites  ·  all effects are predicted accessibility, not measured", options: { fontSize: 13.5, color: "595959", breakLine: true } },
  ], { placeholder: "body" });
  s.addNotes("Everything here is in silico: AlphaGenome's predicted ATAC signal after we edit the sequence. "
    + "The head was fine-tuned on our own 1 bp Tn5 insertion data with Borzoi fold-3 region-based splits, so the numbers quoted are on held-out regions. "
    + "Site classes are the TOBIAS categories as given - Stable = bound in all three conditions; Dynamic = anything else; Age Up / Age Down are the directional age classes.");
}

// ================================================================== 2. the three experiments
{
  const s = pres.addSlide({ masterName: "TITLE_BODY" });
  s.addText("Three perturbations, all in the native 131 kb window", { placeholder: "title" });
  s.addText(bullets([
    { text: "Transplant — take the element out and drop it into 20 closed, non-open backgrounds", bold: true },
    { text: "scored against the same sequence shuffled, so the host contributes nothing. 50 donors per class.", sub: true },
    { text: "Swap — put every element into every other site's location, at two widths", bold: true },
    { text: "62,250 swaps; separates what the element IS from WHERE it sits.", sub: true },
    { text: "Shuffle, motif by motif — scramble the site's own AP-1 motif, and each partner motif, one at a time", bold: true },
    { text: "each partner is scored against a footprint-free window of the same width at the mirrored distance, so the comparison is inside one site and needs no accessibility matching.", sub: true },
    { text: "Why that last control matters", bold: true },
    { text: "Dynamic sites are about half as accessible as Stable ones. Any comparison that does not hold that constant will find a difference whether or not one exists - two results below are exactly that trap.", sub: true },
  ]), { placeholder: "body" });
  s.addNotes("TOBIAS is used only to say where the partner footprints are and which family each belongs to. "
    + "No footprint score enters any calculation; every number on these slides is an AlphaGenome effect.");
}

// ================================================================== 3. transplant
resultSlide({
  title: "A Stable element opens closed chromatin on its own",
  method: "Each element moved into 20 closed backgrounds and scored against the same sequence shuffled; one point per donor, whole-peak width. Black, median with a bootstrap 95% CI; Mann–Whitney.",
  figure: "Analysis/slide_figures/r_fb_transplant.png",
  takeaway: "Stable AP-1 elements carry enough sequence to build their own accessible site. Dynamic ones do not — "
    + "they need the neighbourhood they came from. That is \"context, not sequence autonomy\", shown by moving the sequence.",
  notes: "Per-donor medians over 20 hosts: Stable +4.65 log2 (25x), Dynamic +2.35 (5.1x), p = 4e-5.\n\n"
    + "By class: Age Down +2.62, Age Up +2.66, Exercise Down +2.20, Exercise Up +2.14, Stable +4.65.\n\n"
    + "IMPORTANT CAVEAT, be ready for this question: the donors are not matched on accessibility - Stable donors come from loci "
    + "about 1.8x more accessible (measured ATAC: 32.2 vs 18.0 CPM/kb). Matching on MEASURED fibroblast ATAC (model-independent), "
    + "the residual Stable advantage is +0.84 +- 0.30 log2, i.e. about 1.8x, p = 0.005 (covariance adjustment, 224 donors). "
    + "The 47 matched pairs point the same way (+0.99) but are underpowered (p = 0.13). "
    + "So: direction solid, magnitude about a third to a half of the raw gap. Do not defend 25x vs 5x as a sequence effect.\n\n"
    + "Also worth knowing: 0 of 50 Stable donors lack a measured ATAC peak, against 5-8 per Dynamic class.",
  caveat: "Stable donors also come from more accessible loci. Holding that constant (measured ATAC, 224 donors) a Stable element still gains 1.8× more than a Dynamic one, P = 0.005 — the direction holds, but about a third of the gap above is accessibility rather than sequence.",
});

// ================================================================== 4. scale
resultSlide({
  title: "The element's identity lives at peak scale, not in 170 bp",
  method: "Every element swapped into every other site's location: the swapped site's accessibility regressed on the donor element and on the recipient location. 62,250 swaps; 95% CI from a cluster bootstrap over the 250 donors.",
  figure: "Analysis/slide_figures/r_fb_scale.png",
  takeaway: "At 170 bp an element mostly takes on the level of wherever you put it; at full peak width it imposes its own. "
    + "The information that opens chromatin is spread over the whole peak — which is why a 170 bp assay cannot see it.",
  notes: "62,250 swaps. 170 bp: beta(element) 0.24 vs beta(location) 0.61, R2 = 0.20. Whole peak: 0.71 vs 0.20, R2 = 0.46. "
    + "Transfer fraction 0.29 at 170 bp, 0.83 at peak.\n\n"
    + "Link to the MPRA, if asked: the MPRA library is 170 bp, and its class difference there is tiny (flank shuffle 0.128 vs 0.016 log2, d = 0.14). "
    + "Both assays agree that 170 bp is a poor window for this question. They are NOT the same measurement - the MPRA is episomal activity, "
    + "this is chromatin accessibility in genomic context.\n\n"
    + "AVOID on a slide: our in-situ 170 bp shuffle has Stable sites losing more (-2.51 vs -1.96, p = 4e-40), which is the opposite "
    + "direction to the MPRA flank result. Different assay, but it reads as a contradiction. Keep it off unless someone asks.",
});

// ================================================================== 5. specificity
resultSlide({
  title: "Only the AP-1-losing sites let go — and only of AP-1",
  method: "The site's own AP-1 motif and each partner motif shuffled separately inside the same window; partners net of a width- and distance-matched footprint-free control. Violin, all 2,352 sites; point, median with a bootstrap 95% CI; Mann–Whitney.",
  figure: "Analysis/slide_figures/r_fb_specificity.png",
  takeaway: "Sites that lose AP-1 binding with age were never strongly AP-1-driven: they depend 36% less on their own motif "
    + "while using every other motif exactly as much. Sites that gain AP-1 with age are indistinguishable from Stable.",
  notes: "2,352 sites. Own AP-1 motif: Stable -1.21 (n=1179), AP-1 gained -1.19 (n=605, p = 0.93 vs Stable), AP-1 lost -0.78 (n=209, p = 2.3e-5 vs Stable).\n\n"
    + "Partner motifs, net of control: Stable -0.20, gained -0.25, lost -0.22 - no class difference (lost vs Stable p = 0.66). "
    + "So the deficit is AP-1-specific, not general weakness. No single partner family takes over: nothing survives FDR in the per-family split.\n\n"
    + "Agreement: pooled Stable vs Dynamic at the motif is a null, which matches the manuscript's core-scramble result (no class difference). "
    + "Age-gained being on par with Stable is the sequence-level version of 'age-gained reporter activity on par with stable'.\n\n"
    + "Caveat: the own-AP-1 bar has no matched control, because one at distance 0 would overlap the motif. That inflates AP-1 relative "
    + "to the partners in every class equally, so the class comparison is unaffected.",
});

// ================================================================== 6. AP-1 share
resultSlide({
  title: "Stable sites are more AP-1-centric; Dynamic sites spread it around",
  method: "Per site: |own AP-1| / (|own AP-1| + |partner motifs|). Both terms are AlphaGenome shuffle effects measured inside the same window, so accessibility cannot drive the comparison. Violin with box (median, IQR); Mann–Whitney.",
  figure: "Analysis/slide_figures/r_fb_share.png",
  takeaway: "A Stable site puts 80% of its motif dependence on AP-1 itself; a Dynamic site 74%, with the rest spread across partner motifs. "
    + "Modest, but it is the AP-1-centric vs distributed picture, measured rather than inferred.",
  notes: "Stable 0.804 (n=1176) vs Dynamic 0.745 (n=1156), p = 7.5e-11. Every Dynamic class sits below Stable: "
    + "AP-1 gained 0.754 (p = 6e-5), AP-1 lost 0.736 (p = 6e-6), other Dynamic 0.727 (p = 1.5e-8). Kruskal across the four, p = 1.7e-10.\n\n"
    + "Detection check, in case it is asked: footprint calling needs signal, so Stable sites carry more detected partners "
    + "(5.89 vs 5.01 per site, p = 2e-82). Restricting to the 1,808 sites that hit the cap of 6 detected partners, the gap is unchanged: "
    + "0.806 vs 0.736, p = 2.1e-13. Not a detection artefact.\n\n"
    + "Note this is a DIFFERENT statement from the previous slide. There, only the AP-1-losing class drops in absolute AP-1 dependence. "
    + "Here, all Dynamic classes are modestly less AP-1-centric in relative terms.\n\n"
    + "The equivalent tile-scan measure (motif share vs a 10-50 bp collar) does NOT survive accessibility matching (p = 0.11); "
    + "the collar is mostly noise at low-accessibility sites. This version is better powered because it only scores footprinted motifs.",
});

// ================================================================== 7. what did not hold
{
  const s = pres.addSlide({ masterName: "TITLE_BODY" });
  s.addText("Two differences that turned out to be accessibility", { placeholder: "title" });
  s.addText(bullets([
    { text: "A second AP-1 copy looked worth more to Dynamic sites", bold: true },
    { text: "+0.72 vs +0.60 log2, p < 1e-4 — matched on accessibility, p = 0.93. Gone.", sub: true },
    { text: "Upgrading an imperfect core to a perfect TGAGTCA looked worth more to Dynamic sites", bold: true },
    { text: "+0.66 vs +0.52, p < 1e-4 — matched, p = 0.94. Gone.", sub: true },
    { text: "What did survive matching", bold: true },
    { text: "Destroying the motif costs Dynamic sites less (+0.16, p = 4e-4); the dose-response is shallower; relocating the motif helps them less.", sub: true },
    { text: "Two clean negatives worth having", bold: true },
    { text: "Reverse-complementing the whole motif costs nothing in either class (-0.03 log2, p = 0.89) — expected for a palindromic TRE, and a sanity check on the head.", sub: true },
    { text: "The spacing preference for a second copy is identical: best at ±20–40 bp, decaying by ±150 bp, with no class difference at any distance (all FDR > 0.58).", sub: true },
  ]), { placeholder: "body" });
  s.addNotes("From the fibroblast titration, 2,359 sites, 860 accessibility-matched pairs (caliper 0.1 log10). "
    + "This slide is the argument for why the matched analysis is worth doing - two of the five largest raw class differences "
    + "in the titration were entirely explained by the accessibility gap.\n\n"
    + "Age Down is the outlier again, consistently with the other arms: one base of the core changed -0.62 vs Stable -0.83; motif gone -0.72 vs -1.07.");
}

// ================================================================== 8. where this sits
{
  const s = pres.addSlide({ masterName: "TITLE_BODY" });
  s.addText("How this sits with the MPRA work", { placeholder: "title" });
  s.addText(bullets([
    { text: "Agrees", bold: true },
    { text: "Dynamic elements have weaker endogenous chromatin, and depend on context rather than being sequence-autonomous.", sub: true },
    { text: "Shuffling the AP-1 core alone does not separate pooled Stable from pooled Dynamic.", sub: true },
    { text: "Age-gained elements behave like Stable ones — here, at their own AP-1 motif.", sub: true },
    { text: "C/EBP is the strongest partner motif at fibroblast AP-1 sites (net −0.71 log2, ahead of MEF2 and NFI).", sub: true },
    { text: "Adds", bold: true },
    { text: "The class difference is a whole-peak property, not a 170 bp one.", sub: true },
    { text: "Within Dynamic, the AP-1-losing sites are the ones that are sequence-different; the gaining ones are not.", sub: true },
    { text: "Keep off the slides", bold: true },
    { text: "Our in-situ 170 bp shuffle costs Stable sites more (−2.51 vs −1.96, p = 4e-40), the opposite of the MPRA flank result. Different assay — episomal activity vs chromatin accessibility — but it reads as a contradiction.", sub: true },
  ]), { placeholder: "body" });
  s.addNotes("If the flank question comes up: the MPRA destroys the flanks of a 170 bp fragment on a plasmid and reads transcription; "
    + "we move an intact ~0.8 kb element into closed chromatin and read accessibility. There is no chromatin to open in the MPRA. "
    + "The honest shared conclusion is that neither assay finds much class information in 170 bp.\n\n"
    + "Everything here is predicted accessibility from a model fine-tuned on our own fibroblast ATAC, not a measurement, and not the "
    + "MPRA-trained probe used in the manuscript.");
}

pres.writeFile({ fileName: OUTFILE }).then(() => {
  applyTheme(OUTFILE, THEME);
  console.log("wrote", OUTFILE, `(${pres.slides.length} slides)`);
});
