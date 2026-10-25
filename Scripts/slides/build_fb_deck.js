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
const DIA = require("./diagram_shapes.js");
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
  s.addShape(pres.ShapeType.rect, { x: 0.55, y: y0, w: 0.055, h: caveat ? 0.8 : 0.95, fill: { color: "3F3F3F" }, line: { width: 0 } });
  s.addText(takeaway, { isTextBox: true, x: 0.82, y: y0, w: 11.9, h: caveat ? 0.8 : 0.95, fontSize: 17, color: "111111", valign: "middle", margin: 0 });
  if (caveat) s.addText(caveat, { isTextBox: true, x: 0.82, y: 6.40, w: 11.9, h: 0.66, fontSize: 13, color: "595959", valign: "top", margin: 0 });
  s.addNotes(notes);
  return s;
}

// a method slide. The diagram is drawn as native PowerPoint objects - text boxes, lines, rounded rectangles and
// freeform curves - so every part of it can be moved, recoloured or deleted in PowerPoint without regenerating.
function diagramSlide({ title, draw, method, takeaway, notes }) {
  const s = pres.addSlide({ masterName: "RESULT" });
  s.addText(title, { placeholder: "title" });
  if (method) s.addText(method, { isTextBox: true, x: 0.55, y: 1.06, w: 12.2, h: 0.34, fontSize: 12.5, italic: true, color: "595959", margin: 0 });
  const top = method ? 1.46 : 1.26;
  const dh = ((draw.designHeight || 400) * 12.2) / 1100;        // centre the drawing in the space it has
  if (dh > 6.45 - top) console.warn(`  ! ${title}: drawing is ${dh.toFixed(2)} in tall, overruns the takeaway`);
  draw(pres, s, { x: 0.55, y: top + Math.max(0, (6.45 - top - dh) / 2), w: 12.2 });
  s.addShape(pres.ShapeType.rect, { x: 0.55, y: 6.55, w: 0.055, h: 0.62, fill: { color: "3F3F3F" }, line: { width: 0 } });
  s.addText(takeaway, { isTextBox: true, x: 0.82, y: 6.55, w: 11.9, h: 0.62, fontSize: 15.5, color: "111111", valign: "middle", margin: 0 });
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

// ================================================================== 2. what the model is
diagramSlide({
  title: "What AlphaGenome is, and what it is doing here",
  draw: DIA.d0Model,
  takeaway: "A network that reads 131 kb of raw sequence and returns predicted assay signal at every base. It is "
    + "given no experimental data at prediction time — which is exactly why editing the sequence and predicting "
    + "again counts as an experiment.",
  notes: "Say this before any result, because every number in the deck depends on it.\n\n"
    + "WHAT IT IS: a sequence-to-function model (Google DeepMind, 2025). Convolutional layers pick up local "
    + "motif-scale features, a transformer lets every position see the whole 131,072 bp window, and a stack of "
    + "decoders upsample back to base resolution, where a head per track turns it into predicted signal for "
    + "thousands of assay x cell-type combinations.\n\n"
    + "ON THE WINDOW: AlphaGenome takes a variable context and will go up to 1 Mb. We run it at 131,072 bp, which is "
    + "wide enough to carry the regulatory neighbourhood of a site and cheap enough to repeat tens of thousands of "
    + "times. Every perturbation in the deck uses that same width, so nothing here turns on the choice.\n\n"
    + "WHAT IT IS GIVEN: one-hot DNA. No ATAC, no expression, no annotation, no cell-type label beyond our choice of "
    + "which output track to read. If someone asks whether the model 'already knows' a site is dynamic - it cannot; "
    + "it never sees the labels, the conditions, or the TOBIAS calls.\n\n"
    + "WHAT WE READ: the fibroblast ATAC head, fine-tuned on our own Tn5 data. Held-out Spearman 0.68 against 0.27 "
    + "for the built-in head, which is the reason we fine-tuned rather than using it off the shelf. The FB trackset "
    + "is the mean of the four conditions.\n\n"
    + "WHY IT LICENSES THE REST: because the output is a deterministic function of the input sequence, two "
    + "predictions that differ only in an edit differ only because of the edit. That is the whole logic of the next "
    + "slide, and of every perturbation in the deck.\n\n"
    + "THE LIMIT, say it here rather than be asked: these are predictions, not measurements. AlphaGenome is good at "
    + "ranking sequences on accessibility; it is not a reporter assay. Anything surprising is a hypothesis for the "
    + "bench, not a result from it.",
});

// ================================================================== 3. how an effect is measured
diagramSlide({
  title: "Everything here is one prediction minus another",
  draw: DIA.d1Pipeline,
  takeaway: "No edit is ever compared against another site. The native and the edited sequence are the same 131 kb "
    + "window through the same weights, so whatever changes in the prediction came from the edit.",
  notes: "Worth saying out loud at the start, because it is what makes the rest defensible: the comparison is always "
    + "within one window, and the model never sees a label.\n\n"
    + "The head is ours - fine-tuned on our own 1 bp fibroblast Tn5 insertions with Borzoi fold-3 region-based splits, "
    + "so 0.68 is on held-out regions. The built-in AlphaGenome ATAC head manages 0.27 on the same data, which is why "
    + "none of this was run on it.\n\n"
    + "FB is the mean of the four conditions (young/old x sedentary/exercise), not a single track.\n\n"
    + "TOBIAS is used only to say WHERE a footprint is and which family it belongs to. No footprint score enters any "
    + "calculation on these slides.",
});

// ================================================================== 3. transplant, the design
diagramSlide({
  title: "Transplant: can the element build a site on its own?",
  draw: DIA.d2Transplant,
  takeaway: "The control is the same letters shuffled, dropped into the same background — so whatever the host "
    + "contributes, it contributes to both and cancels in the ratio.",
  notes: "Hosts are genomic windows with no called ATAC peak. 20 per element, 50 donor elements per class, scored "
    + "at the whole-peak width.\n\n"
    + "Why the shuffled control and not an empty host: an empty host still has its own baseline signal, and that "
    + "baseline differs between hosts. Scoring element against shuffled inside the SAME host removes it.\n\n"
    + "What this cannot do: the donors are not matched on accessibility. That is handled on the next slide and the "
    + "caveat line there carries the matched number.",
});

// ================================================================== 4. transplant
resultSlide({
  title: "A Stable element opens closed chromatin on its own",
  method: "Each element moved into 20 closed backgrounds and scored against the same sequence shuffled; one point per donor, whole-peak width. Box, median and IQR, whiskers 1.5 × IQR; 50 donors per class — Up = AP-1 gained, Down = AP-1 lost; Mann–Whitney, Stable against the four Dynamic classes pooled.",
  figure: "Analysis/slide_figures/r_fb_transplant.png",
  takeaway: "Stable AP-1 elements carry enough sequence to build their own accessible site. Dynamic ones do not — "
    + "they need the neighbourhood they came from. That is \"context, not sequence autonomy\", shown by moving the sequence.",
  notes: "Per-donor medians over 20 hosts: Stable +4.65 log2 (25x), Dynamic pooled +2.35 (5.1x), p = 1.5e-5.\n\n"
    + "By class, as plotted: Stable 25x, Age Up 6.3x, Age Down 6.1x, Exercise Down 4.6x, Exercise Up 4.4x. "
    + "Every Dynamic class sits well below Stable and they are close to each other, which is why pooling them is fair here.\n\n"
    + "IMPORTANT CAVEAT, be ready for this question: the donors are not matched on accessibility - Stable donors come from loci "
    + "about 1.8x more accessible (measured ATAC: 32.2 vs 18.0 CPM/kb). Matching on MEASURED fibroblast ATAC (model-independent), "
    + "the residual Stable advantage is +0.84 +- 0.30 log2, i.e. about 1.8x, p = 0.005 (covariance adjustment, 224 donors). "
    + "The 47 matched pairs point the same way (+0.99) but are underpowered (p = 0.13). "
    + "So: direction solid, magnitude about a third to a half of the raw gap. Do not defend 25x vs 5x as a sequence effect.\n\n"
    + "Also worth knowing: 0 of 50 Stable donors lack a measured ATAC peak, against 5-8 per Dynamic class.",
  caveat: "Stable donors also come from more accessible loci. Holding that constant (measured ATAC, 224 donors) a Stable element still gains 1.8× more than a Dynamic one, P = 0.005 — the direction holds, but about a third of the gap above is accessibility rather than sequence.",
});

// ================================================================== 5. swap, the design
diagramSlide({
  title: "Swap: is the level set by the element, or by the place?",
  draw: DIA.d3Swap,
  takeaway: "Write every element into every other element's location and regress the result on both. "
    + "Whichever term carries more weight is the one that sets a site's accessibility.",
  notes: "62,250 swaps, self-swaps excluded, at two widths: the 170 bp MPRA-sized window and the whole peak.\n\n"
    + "Each donor appears in about 250 of those rows, so the rows are nowhere near independent. The intervals on the "
    + "next slide come from a cluster bootstrap over the 250 donors; an OLS standard error on 62,250 rows would be "
    + "far too small and would make both coefficients look certain.",
});

// ================================================================== 6. element or place
resultSlide({
  title: "The element's identity lives at peak scale, not in 170 bp",
  method: "Every element written into every other site's location. The accessibility that results is regressed on two things: how accessible the donor element's own site was, and how accessible the recipient location was. Both weights are positive — they are drawn on opposite sides of zero only so the pair can be compared, and the longer bar is the stronger influence. P tests the two weights against each other. 62,250 swaps; 95% CI from a cluster bootstrap over the 250 donors.",
  figure: "Analysis/slide_figures/r_fb_scale.png",
  takeaway: "At 170 bp an element mostly takes on the level of wherever you put it; at full peak width it imposes its own. "
    + "The information that opens chromatin is spread over the whole peak — which is why a 170 bp assay cannot see it.",
  caveat: "Read it as a transplant. Take a strong element and put it somewhere weak: at 170 bp it comes out weak, because the location carries most of the weight (0.61 against 0.24). At full peak width it stays strong, because the element now carries most of it (0.71 against 0.20). Nothing about the element changed — only how much of it we moved.",
  notes: "62,250 swaps. 170 bp: beta(element) 0.24 vs beta(location) 0.61, R2 = 0.20. Whole peak: 0.71 vs 0.20, R2 = 0.46. "
    + "Transfer fraction 0.29 at 170 bp, 0.83 at peak.\n\n"
    + "Link to the MPRA, if asked: the MPRA library is 170 bp, and its class difference there is tiny (flank shuffle 0.128 vs 0.016 log2, d = 0.14). "
    + "Both assays agree that 170 bp is a poor window for this question. They are NOT the same measurement - the MPRA is episomal activity, "
    + "this is chromatin accessibility in genomic context.\n\n"
    + "AVOID on a slide: our in-situ 170 bp shuffle has Stable sites losing more (-2.51 vs -1.96, p = 4e-40), which is the opposite "
    + "direction to the MPRA flank result. Different assay, but it reads as a contradiction. Keep it off unless someone asks.",
});

// ================================================================== 7. motif shuffle, the design
diagramSlide({
  title: "Shuffling one motif at a time, inside the site",
  draw: DIA.d4Shuffle,
  takeaway: "Each partner motif is scored against a footprint-free window of the same width at the mirrored "
    + "distance, so the two differ only in whether a motif was there. Everything stays inside one site.",
  notes: "This is the control that makes the next two slides immune to the accessibility gap. Dynamic sites are "
    + "about half as accessible as Stable ones, so any measure that compares a Dynamic number against a Stable "
    + "number is exposed to that confound. Here both terms come from the same window of the same site.\n\n"
    + "TOBIAS says only WHERE a footprint is and which family it belongs to. No footprint score enters the "
    + "calculation; the effect is always AlphaGenome's.\n\n"
    + "One asymmetry to admit if asked: the site's own AP-1 motif has no matched control, because a control at "
    + "distance 0 would overlap the motif. That inflates AP-1 relative to the partners in every class equally, so "
    + "the comparison BETWEEN classes is unaffected.",
});

// ================================================================== 8. motif specificity
resultSlide({
  title: "The sites that lose AP-1 let go of it — whether with age or exercise",
  method: "The site's own AP-1 motif and each partner motif shuffled separately inside the same window; partners net of a width- and distance-matched footprint-free control. Box, median and IQR, whiskers 5th–95th centile, notch a 95% interval on the median. Classes as given: Up = AP-1 gained, Down = AP-1 lost.",
  figure: "Analysis/slide_figures/r_fb_specificity.png",
  takeaway: "It is the direction that matters, not the stimulus. Both losing classes depend less on their own AP-1 motif "
    + "— Age Down by 36%, Exercise Down by 21% — while both gaining classes are indistinguishable from Stable. "
    + "No class differs at its partner motifs.",
  notes: "2,340 sites. Own AP-1 motif, median log2, with BH-adjusted P against Stable:\n"
    + "  Stable -1.21 (n=1179, reference) · Age Up -1.19 (n=605, FDR 0.93) · Age Down -0.78 (n=209, FDR 9.1e-5) · "
    + "Exercise Up -1.14 (n=233, FDR 0.49) · Exercise Down -0.96 (n=114, FDR 0.043)\n\n"
    + "Partner motifs, net of control: -0.20 / -0.25 / -0.22 / -0.23 / -0.33. Nothing survives BH (smallest adjusted 0.20), "
    + "so the deficit is AP-1-specific in both losing classes. "
    + "So the deficit is AP-1-specific, not general weakness. No single partner family takes over: nothing survives FDR in the per-family split.\n\n"
    + "Agreement: pooled Stable vs Dynamic at the motif is a null, which matches the manuscript's core-scramble result (no class difference). "
    + "Age-gained being on par with Stable is the sequence-level version of 'age-gained reporter activity on par with stable'.\n\n"
    + "Note on the split: a site that is dynamic with both age and exercise is counted under its age class, so the two "
    + "exercise groups are the sites that are NOT also age-dynamic. Splitting them is what separated the old 'other "
    + "Dynamic' group (p = 0.055, n = 347) into Exercise Up, which looks like Stable, and Exercise Down, which does not.\n\n"
    + "Caveat: the own-AP-1 bar has no matched control, because one at distance 0 would overlap the motif. That inflates AP-1 relative "
    + "to the partners in every class equally, so the class comparison is unaffected.",
});

// ================================================================== 9. the scan, the design
diagramSlide({
  title: "Scanning the whole window, tile by tile and base by base",
  draw: DIA.d5Scan,
  takeaway: "Rather than choosing a region to destroy, replace every 10 bp tile across 600 bp in turn — and then "
    + "every single base across the middle 61 — and ask where the cost lands.",
  notes: "60 tiles of 10 bp from -300 to +290, one prediction each; then 61 positions x 3 alternative bases.\n\n"
    + "The analysis step that matters is on the next slide: each site's profile has its OWN noise floor subtracted "
    + "(the median |effect| beyond +-200 bp) before being normalised to its own total. Without that, Dynamic sites "
    + "look more 'distributed' simply because they are less accessible and therefore noisier, and noise is flat "
    + "across position.",
});

// ================================================================== 10. where the dependence sits
resultSlide({
  title: "Where the dependence sits: on the motif, and more so if Stable",
  method: "Each 10 bp tile replaced in turn; each site's own noise floor removed, then normalised to its own total, so the curve is where the dependence sits, not how much there is. Median per site, band a bootstrap 95% CI; Mann–Whitney per position, BH-FDR across the 60. Cropped to ±150 bp.",
  figure: "Analysis/slide_figures/r_fb_density.png",
  takeaway: "A spike on the AP-1 core sitting on a shallow shoulder that reaches ±100 bp — and the model was never told "
    + "where the motif was. The spike is taller for Stable sites, and the two tiles covering the core are the only positions "
    + "in the whole window where the classes differ at all.",
  caveat: "Holding accessibility constant the gap halves but holds: 33.1% vs 31.3% of importance within ±25 bp, P = 0.006 over 860 matched pairs. The older ratio version of this measure did not survive that check (P = 0.11); this one does, because the noise floor is removed per site first.",
  notes: "NUMBERS TO ANNOTATE WITH, strongest first:\n"
    + "  - the peak tile carries 10.8% (Stable) / 9.9% (Dynamic) of a site's importance, against 1.7% for an average tile - about 6x enrichment\n"
    + "  - within +-25 bp (5 of 60 tiles, 8% of the sequence): Stable 33.8% vs Dynamic 30.7%, p = 2.0e-6\n"
    + "  - accessibility-matched, 860 pairs: 33.1% vs 31.3%, paired difference +1.5 points, Wilcoxon p = 0.006\n"
    + "  - only 2 of 60 positions differ at FDR < 0.05, and they are the two tiles over the core (FDR = 0.001). Every other position FDR > 0.8\n"
    + "  - in 56% of Stable and 52% of Dynamic sites the single most important tile IS a core tile, against a chance of 3.3% - 16x. "
    + "That one does NOT differ by class (Fisher p = 0.063; matched McNemar p = 0.21), so quote it as the shared fact, not as a difference. "
    + "Figure for it is on disk as r_fb_top_tile.png if anyone wants to see it.\n\n"
    + "2,359 sites, 60 tiles of 10 bp from -300 to +290, fine-tuned FB head, readout atac_170.\n\n"
    + "Why the floor subtraction matters, if asked: Dynamic sites are about half as accessible, so their log2 ratios are noisier, "
    + "and noise is flat across position. Without removing each site's own floor (median |effect| beyond +-200 bp: 0.047 Dynamic, "
    + "0.037 Stable) a Dynamic profile looks more distributed whether or not it is. That is exactly the artefact this corrects for.\n\n"
    + "Agreement: this is the AlphaGenome version of 'attention is more distributed at dynamic enhancers'. Same direction, "
    + "different model and different readout.\n\n"
    + "Honest size: the matched difference is 1.5 percentage points on a base of about a third. Real, repeatedly so, but small - "
    + "do not oversell it as two different architectures.",
});

// ================================================================== 10a-c. the same result, read three other ways
// Three alternative readings of the density result, in the deck so they can be compared side by side. Keep one and delete the
// rest - the measurement is identical in all four, only the plotted quantity changes. ALT marks them in the deck.
const ALT = "ALTERNATIVE VIEW OF THE PREVIOUS SLIDE — same measurement, keep one.  ";

resultSlide({
  title: "Stable sites keep half their dependence within 41 bp, Dynamic 48",
  method: ALT + "Cumulative share of a site's importance within ±x bp, against x. Median per site, band a bootstrap 95% CI. The value at ±20 bp is the 33.8% / 30.7% quoted on the previous slide, so the two readings cannot disagree.",
  figure: "Analysis/slide_figures/r_fb_density_cumulative.png",
  takeaway: "Per-position differences of a point each add up: by ±60 bp the two curves are 5.3 points apart. Half of a "
    + "Stable site's dependence sits within ±41 bp of the motif, against ±48 bp for a Dynamic one — the same result as "
    + "the previous slide, said as a distance rather than a percentage.",
  caveat: "Nothing is rescaled: the axis runs 0 to 100% and both curves end at 1. The half-radii are read off a profile sampled every 10 bp, so quote them as about 40 against about 50 bp until the sliding scan lands.",
  notes: "Stable vs Dynamic, cumulative share within ±x bp:\n"
    + "  ±20 bp  33.8% / 30.7%   (the previous slide's number)\n"
    + "  ±40 bp  49.0% / 44.5%\n"
    + "  ±60 bp  62.8% / 57.5%   <- widest gap, 5.3 points\n"
    + "  ±100 bp 80.0% / 75.7%\n"
    + "Half-maximum radius by linear interpolation: Stable 41.4 bp, Dynamic 48.1 bp.\n\n"
    + "Why this view: the two profiles on the previous slide differ by about a percentage point per position, which is "
    + "invisible when both are drawn at full height. Integrating makes a real difference legible without inflating it - "
    + "the quantity at every radius is still a share of the site's own total.\n\n"
    + "The honest size is unchanged from the previous slide: about 1.5 points once accessibility is matched, on a base "
    + "of a third. Do not let the separation in this plot be heard as two different architectures.",
});

resultSlide({
  title: "The classes differ on the motif and nowhere else",
  method: ALT + "Stable minus Dynamic at each position, on its own axis. Band from resampling sites within each class; filled points are FDR < 0.05, Mann–Whitney with BH across the 60 tiles.",
  figure: "Analysis/slide_figures/r_fb_density_diff.png",
  takeaway: "One curve instead of two overlapping ones. It is flat across the entire 600 bp except for a single spike of "
    + "+0.94 points on the tile covering the AP-1 core — the difference between the classes is not spread along the "
    + "element, it is on the motif itself.",
  notes: "Two of 60 positions reach FDR < 0.05 and they are the two core tiles: -10 bp (+0.94 points, FDR 0.0008) and "
    + "0 bp (+0.84 points, FDR 0.0007). The next largest difference anywhere is +0.31 points at -20 bp, FDR 0.81; every "
    + "other position FDR > 0.9.\n\n"
    + "Use this one if the point you want to land is WHERE the classes differ rather than BY HOW MUCH. It is the most "
    + "conservative of the three - the whole class difference is two points above a flat line, and the plot says so.\n\n"
    + "A difference plot hides the shared signal: both classes still put a large spike on the motif, and that spike is "
    + "the main result of the previous slide. Say it out loud if you show this instead.",
});

resultSlide({
  title: "The two classes that lose AP-1 spread their dependence furthest",
  method: ALT + "Share of a site's importance within ±25 bp of the motif, by class. Point, median per site; bar, bootstrap 95% CI; Mann–Whitney against Stable, BH across the four comparisons.",
  figure: "Analysis/slide_figures/r_fb_density_classes.png",
  takeaway: "All four Dynamic classes are less motif-centred than Stable, and they order the same way the motif-shuffle "
    + "slide did: the two AP-1-losing classes (28.3%, 28.6%) sit below the two gaining ones (31.3%, 31.5%), which sit "
    + "below Stable at 33.8%. Two independent measures, the same split.",
  caveat: "Exercise Down is 115 sites and its interval (24.2–32.2%) overlaps the gaining classes, so the ordering within Dynamic is suggestive, not established. Every class is below Stable, and that part is solid.",
  notes: "Share within ±25 bp, median per site, with BH-adjusted P against Stable:\n"
    + "  Stable         33.8%  (n=1179, reference)\n"
    + "  Age Up         31.3%  (n=613,  FDR 0.0038)\n"
    + "  Age Down       28.3%  (n=212,  FDR 0.00045)\n"
    + "  Exercise Up    31.5%  (n=240,  FDR 0.0087)\n"
    + "  Exercise Down  28.6%  (n=115,  FDR 0.0084)\n\n"
    + "This is the only one of the three alternatives that shows something the pooled comparison cannot. It is also the "
    + "one that connects to the motif-specificity slide: there, the LOSING classes let go of the motif in absolute terms; here, the same two "
    + "classes are the ones whose dependence is least concentrated on it. Same partition, different measurement.\n\n"
    + "What it is NOT: this is not accessibility. The matched analysis on the previous slide (860 pairs, p = 0.006) holds "
    + "for the pooled contrast; the per-class split has not been matched, because the smaller classes do not support it. "
    + "If pushed, that is the honest answer - the pooled difference survives matching, the ordering within Dynamic is "
    + "unmatched and n=115 at the bottom.\n\n"
    + "A site dynamic with both age and exercise is counted under its age class, as everywhere else in the deck.",
});

// ================================================================== 11. how the dependence is divided
resultSlide({
  title: "Dynamic sites hand more of the load to the partner motifs",
  method: "Per site: |own AP-1| / (|own AP-1| + |partner motifs|), the share of its motif dependence that sits inside the AP-1 motif rather than outside it. Both terms are AlphaGenome shuffle effects measured inside the same window, so accessibility cannot drive the comparison. Median per site; error bar, bootstrap 95% CI; Mann–Whitney.",
  figure: "Analysis/slide_figures/r_fb_share.png",
  takeaway: "Outside the motif is where every class differs: a Stable site leaves 20% of its motif dependence to "
    + "partners, each of the four Dynamic classes 25–27% — about a third more — and all four sit on the same side, "
    + "gaining and losing alike. AP-1 still leads everywhere; it just carries less of the load.",
  notes: "Stable 0.804 (n=1179) vs Dynamic pooled 0.745 (n=1161), p = 7.2e-11. By class: Age Up 0.755 (p = 7.9e-5), "
    + "Age Down 0.734 (p = 4.1e-6), Exercise Up 0.726 (p = 1.9e-5), Exercise Down 0.732 (p = 8.7e-6).\n\n"
    + "Worth contrasting with the previous slide out loud: there, only the LOSING classes moved. Here all four move, "
    + "including both gaining classes. So the absolute dependence on AP-1 and the share of the total it carries are "
    + "two different statements, and only the second separates every Dynamic class from Stable.\n\n"
    + "The 40%: partner share 19.6% for Stable against 24.5 / 26.6 / 27.3% for the Dynamic classes, so 27.3/19.6 = 1.39x. "
    + "It is the same number as 80 vs 73, read from the other end - say it whichever way, but do not say both as if they were two results.\n\n"
    + "Detection check, in case it is asked: footprint calling needs signal, so Stable sites carry more detected partners "
    + "(5.89 vs 5.01 per site, p = 2e-82). Restricting to the 1,808 sites that hit the cap of 6 detected partners, the gap is unchanged: "
    + "0.806 vs 0.736, p = 2.1e-13. Not a detection artefact.\n\n"
    + "Note this is a DIFFERENT statement from the motif-specificity slide. There, only the AP-1-losing class drops in absolute AP-1 dependence. "
    + "Here, all Dynamic classes are modestly less AP-1-centric in relative terms.\n\n"
    + "The equivalent tile-scan measure (motif share vs a 10-50 bp collar) does NOT survive accessibility matching (p = 0.11); "
    + "the collar is mostly noise at low-accessibility sites. This version is better powered because it only scores footprinted motifs.",
});

// ================================================================== 12. the AP-1 ladder
diagramSlide({
  title: "The AP-1 ladder, and the spacing experiments",
  draw: DIA.d6Titration,
  takeaway: "Rather than only destroying the motif, weaken it one base at a time, strengthen it, flip it, duplicate "
    + "it or move it — each against a control that changes the same amount of sequence but no motif.",
  notes: "This is where the next slide's numbers come from. 2,359 sites; 860 accessibility-matched pairs at caliper "
    + "0.1 log10.\n\n"
    + "The ladder (1 to 7 bases of the core) is what gives a dose-response rather than a single on/off contrast - it "
    + "is how we can say the Dynamic curve is SHALLOWER rather than just lower.\n\n"
    + "Two controls worth naming if asked:\n"
    + "  - the second copy is paired with a scrambled insert of the same width at the same distance, so the comparison "
    + "is 'a motif there' against 'that much sequence there'\n"
    + "  - moving the motif is paired with moving it to a spot that has no footprint\n\n"
    + "The reverse complement is the sanity check on the head: a TRE is a near-palindrome, only 5 of the 11 bases "
    + "actually move, and the predicted cost is -0.03 log2 (p = 0.89). If that had come out large, nothing else here "
    + "would be trustworthy.",
});

// ================================================================== 13. what did not hold
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

// ================================================================== 14. where this sits
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
    { text: "Within Dynamic, the AP-1-LOSING sites are the ones that are sequence-different, with age and with exercise alike; the gaining ones look like Stable. Direction, not stimulus.", sub: true },
    { text: "Unprompted, AlphaGenome makes the AP-1 core the most important 10 bp in over half of all sites — equally in both classes. The classes differ in how much of the dependence sits there, not in whether AP-1 leads.", sub: true },
    { text: "Dynamic sites spread their dependence slightly wider, and it survives accessibility matching (33.1% vs 31.3% within ±25 bp, P = 0.006) — the sequence-model version of “attention is more distributed at dynamic enhancers”.", sub: true },
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
