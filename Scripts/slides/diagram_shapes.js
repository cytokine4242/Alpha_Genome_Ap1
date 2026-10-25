// The method diagrams, drawn as native PowerPoint objects rather than placed as images.
//
// Every element is editable in PowerPoint: text is a real text box, the DNA is a line, the motif and the scrambled
// stretches are rounded rectangles, and the predicted ATAC tracks are freeform shapes (custGeom), so they can be
// recoloured, moved, restyled or nudged without coming back to the generator.
//
// Layout is written in design units on an 1100-wide canvas, the same numbers the earlier SVG version used; U()
// converts to inches. Font sizes are in points and therefore do NOT scale with the drawing - that is deliberate,
// so the lettering stays legible whatever width the diagram is placed at.

const INK = "15181C";
const MUTED = "5A6471";
const HAIR = "AFB5BD";
const SIG = "55637A";       // predicted ATAC signal
const GOLD = "D9A32B";      // the AP-1 motif
const GOLD_D = "8A6400";
const EDIT = "E2622F";      // sequence we changed
const EDIT_D = "A2380F";
const PART = "4E8E8A";      // a partner motif

// ------------------------------------------------------------------ frame
// Maps the 1100-unit design canvas onto a rectangle of the slide.
function frame(pres, slide, { x, y, w, vbw = 1100 }) {
  const k = w / vbw;
  const F = {
    pres, slide, k, x0: x, y0: y,
    X: (u) => x + u * k,
    Y: (u) => y + u * k,
    S: (u) => u * k,
  };

  // --- text: always a real text box, anchored on its own baseline-ish centre
  F.txt = (s, u, v, o = {}) => {
    const bw = o.bw || 260, bh = o.bh || 20;
    const anchor = o.anchor || "start";
    const left = anchor === "middle" ? F.X(u) - F.S(bw) / 2 : anchor === "end" ? F.X(u) - F.S(bw) : F.X(u);
    slide.addText(s, {
      x: left, y: F.Y(v) - F.S(bh) / 2, w: F.S(bw), h: F.S(bh),
      fontFace: o.face || "Arial", fontSize: o.size || 12, bold: !!o.bold,
      color: o.color || INK, align: anchor === "middle" ? "center" : anchor === "end" ? "right" : "left",
      valign: "middle", margin: 0, isTextBox: true, wrap: o.wrap !== false,
      charSpacing: o.ls || 0, lineSpacingMultiple: o.lsm || 0.95,
    });
  };
  F.head = (s, u, v, o = {}) => F.txt(s, u, v, { size: 13, bold: true, bw: 420, ...o });
  F.cap = (s, u, v, o = {}) => F.txt(s, u, v, { size: 11, color: MUTED, ...o });

  // --- line
  F.line = (u1, v1, u2, v2, o = {}) => {
    const x1 = F.X(Math.min(u1, u2)), y1 = F.Y(Math.min(v1, v2));
    slide.addShape(pres.ShapeType.line, {
      x: x1, y: y1, w: F.S(Math.abs(u2 - u1)), h: F.S(Math.abs(v2 - v1)),
      flipH: u2 < u1, flipV: v2 < v1,
      line: {
        color: o.color || HAIR, width: o.w || 1,
        dashType: o.dash || "solid",
        ...(o.arrow ? { endArrowType: "triangle" } : {}),
        ...(o.transparency !== undefined ? { transparency: o.transparency } : {}),
      },
    });
  };

  // --- rounded block sitting on the DNA (the motif, or a scrambled stretch)
  F.block = (u, v, w_, o = {}) => {
    const h_ = o.h || 11;
    slide.addShape(pres.ShapeType.roundRect, {
      x: F.X(u) - F.S(w_) / 2, y: F.Y(v) - F.S(h_) / 2, w: F.S(w_), h: F.S(h_),
      rectRadius: 0.02,
      fill: o.fill === false ? { type: "none" } : { color: o.color || GOLD, transparency: o.t === undefined ? 20 : o.t },
      line: { color: o.line || o.color || GOLD, width: o.lw || 1, dashType: o.dash || "solid" },
    });
  };

  // --- a scrambled stretch: a band with a zig-zag through it, so it needs no label
  F.scramble = (u, v, w_, o = {}) => {
    const col = o.color || EDIT, h_ = o.h || 13;
    F.block(u, v, w_, { color: col, t: 86, line: col, lw: 0.75, h: h_ });
    const n = Math.max(4, Math.round(w_ / 5)), pts = [];
    for (let i = 0; i <= n; i++) pts.push({ x: F.S((w_ * i) / n), y: F.S(i % 2 ? 9 : 1) });
    slide.addShape(pres.ShapeType.custGeom, {
      x: F.X(u) - F.S(w_) / 2, y: F.Y(v) - F.S(5), w: F.S(w_), h: F.S(10),
      points: pts, fill: { type: "none" }, line: { color: col, width: 1.4 },
    });
  };

  // --- the predicted ATAC track, as one freeform shape
  F.track = (u, v, w_, peaks, o = {}) => {
    const col = o.color || SIG, hmax = o.hmax || 62, N = 96;
    const noise = (t, s) => Math.sin(t * 11 + s) * 0.5 + Math.sin(t * 29 + s * 2.1) * 0.3 + Math.sin(t * 53 + s * 3.7) * 0.2;
    const pts = [{ x: 0, y: F.S(hmax) }];
    for (let i = 0; i <= N; i++) {
      const t = i / N;
      let hv = 2.0 + 1.6 * Math.abs(noise(t * 6, o.seed || 1));
      for (const p of peaks) hv += p.a * Math.exp(-Math.pow((t - p.c) * w_, 2) / (2 * p.w * p.w));
      pts.push({ x: F.S(t * w_), y: F.S(hmax - Math.min(hv, hmax)) });
    }
    pts.push({ x: F.S(w_), y: F.S(hmax) }, { close: true });
    slide.addShape(pres.ShapeType.custGeom, {
      x: F.X(u), y: F.Y(v - hmax), w: F.S(w_), h: F.S(hmax), points: pts,
      fill: { color: col, transparency: o.t === undefined ? 86 : o.t },
      line: { color: col, width: 1.3 },
    });
  };

  // --- the DNA itself
  F.dna = (u1, u2, v) => F.line(u1, v, u2, v, { color: INK, w: 1.6, transparency: 45 });

  // --- a measurement bracket: the quantity a panel produces
  F.delta = (u, v1, v2, label, o = {}) => {
    const col = o.color || INK;
    F.line(u, v1, u, v2, { color: col, w: 1.2 });
    F.line(u, v1, u + 7, v1, { color: col, w: 1.2 });
    F.line(u, v2, u + 7, v2, { color: col, w: 1.2 });
    F.txt(label, u + 13, (v1 + v2) / 2, { size: 12, bold: true, color: col, bw: 150 });
  };
  return F;
}

// ------------------------------------------------------------------ shared pieces
// lift is how far above the DNA the label sits: enough to clear whatever the track is doing there
function motifMark(F, u, v, { label = null, w = 26, t = 20, color = GOLD, labelColor = GOLD_D, lift = 24 } = {}) {
  F.block(u, v, w, { color, t, lw: 1.1 });
  if (label) F.txt(label, u, v - lift, { anchor: "middle", size: 11.5, bold: true, color: labelColor, bw: 120 });
}

// ================================================================== d1  the pipeline
function d1Pipeline(pres, slide, box) {
  const F = frame(pres, slide, box);
  const LX = 136, LW = 232, PX = 616, PW = 292;
  const yA = 108, yB = 230;

  F.head("NATIVE", 14, yA - 34, { bw: 220 });
  F.cap("as it is", 14, yA - 16, { bw: 220 });
  F.dna(LX, LX + LW, yA);
  motifMark(F, LX + LW / 2, yA, { label: "AP-1" });

  F.head("EDITED", 14, yB - 34, { bw: 220 });
  F.cap("one thing changed", 14, yB - 16, { bw: 240 });
  F.dna(LX, LX + LW, yB);
  F.scramble(LX + LW / 2, yB, 30);

  slide.addShape(pres.ShapeType.roundRect, {
    x: F.X(432), y: F.Y(yA - 56), w: F.S(132), h: F.S(yB - yA + 92), rectRadius: 0.06,
    fill: { color: INK, transparency: 95 }, line: { color: INK, width: 1, transparency: 68 },
  });
  F.txt("AlphaGenome", 498, yA + 48, { anchor: "middle", size: 13, bold: true, bw: 160 });
  F.cap("1 bp ATAC head,", 498, yA + 72, { anchor: "middle", bw: 160 });
  F.cap("fine-tuned on our", 498, yA + 87, { anchor: "middle", bw: 160 });
  F.cap("own fibroblast Tn5", 498, yA + 102, { anchor: "middle", bw: 170 });

  F.line(LX + LW + 14, yA, 426, yA, { color: INK, w: 1.4, arrow: true });
  F.line(LX + LW + 14, yB, 426, yB, { color: EDIT, w: 1.4, arrow: true });
  F.line(572, yA, PX - 12, yA, { color: INK, w: 1.4, arrow: true });
  F.line(572, yB, PX - 12, yB, { color: EDIT, w: 1.4, arrow: true });

  F.cap("predicted accessibility", PX, yA - 62, { bw: 260 });
  F.track(PX, yA, PW, [{ c: 0.5, w: 26, a: 48 }], { seed: 1.3 });
  F.dna(PX, PX + PW, yA);
  F.track(PX, yB, PW, [{ c: 0.5, w: 26, a: 17 }], { seed: 1.3, color: EDIT, t: 88 });
  F.dna(PX, PX + PW, yB);

  F.delta(PX + PW + 32, yA - 48, yB - 17, "log₂ ratio");
  F.cap("one number", PX + PW + 45, yB + 6, { bw: 170 });
  F.cap("per site", PX + PW + 45, yB + 21, { bw: 170 });

  F.line(14, 288, 1086, 288, { dash: "dash" });
  F.cap("both predictions come from the same 131 kb window, so only the edit can move them", 14, 310, { bw: 900 });
  return 334;
}

// ================================================================== d2  transplant
function d2Transplant(pres, slide, box) {
  const F = frame(pres, slide, box);
  const W = 252, C = [44, 392, 740];
  const yT = 128, yB = 258;

  F.head("1   A REAL SITE", C[0], 40);
  F.track(C[0], yT, W, [{ c: 0.5, w: 32, a: 52 }], { seed: 2.2 });
  F.dna(C[0], C[0] + W, yT);
  motifMark(F, C[0] + W / 2, yT, { label: "AP-1", lift: 70 });
  F.block(C[0] + W / 2, yT + 17, 176, { fill: false, line: INK, lw: 1, dash: "dash", h: 17 });
  F.cap("lift out its whole peak", C[0] + W / 2, yT + 42, { anchor: "middle", bw: 240 });

  F.head("2   AN EMPTY PLACE", C[1], 40);
  F.track(C[1], yT, W, [], { seed: 4.1 });
  F.dna(C[1], C[1] + W, yT);
  F.cap("a closed background", C[1] + W / 2, yT + 24, { anchor: "middle", bw: 240 });
  F.cap("with no called peak", C[1] + W / 2, yT + 40, { anchor: "middle", bw: 240 });

  F.head("3   WHAT IT BUILDS THERE", C[2], 40);
  F.track(C[2], yT, W, [{ c: 0.5, w: 28, a: 46 }], { seed: 2.6 });
  F.dna(C[2], C[2] + W, yT);
  motifMark(F, C[2] + W / 2, yT);
  F.txt("the element", C[2], yT + 22, { size: 12, bold: true, bw: 220 });

  F.track(C[2], yB, W, [{ c: 0.5, w: 28, a: 6 }], { seed: 2.6, color: EDIT, t: 88 });
  F.dna(C[2], C[2] + W, yB);
  F.scramble(C[2] + W / 2, yB, 30);
  F.txt("the same letters, shuffled", C[2], yB + 22, { size: 12, color: EDIT_D, bw: 260 });
  F.cap("the host is in both, so it cancels", C[2], yB + 44, { bw: 300 });

  F.delta(C[2] + W + 24, yT - 46, yB - 8, "effect");
  F.line(C[0] + W + 18, yT - 16, C[1] - 18, yT - 16, { color: MUTED, w: 1.2, arrow: true, transparency: 35 });
  F.line(C[1] + W + 18, yT - 16, C[2] - 18, yT - 16, { color: MUTED, w: 1.2, arrow: true, transparency: 35 });

  F.line(14, 318, 1086, 318, { dash: "dash" });
  F.cap("20 backgrounds per element   ·   50 elements per class   ·   whole-peak width", 14, 340, { bw: 800 });
  return 362;
}

// ================================================================== d3  swap
function d3Swap(pres, slide, box) {
  const F = frame(pres, slide, box);
  const AX = 44, AW = 232, BX = 648, BW = 268;
  const yD = 112, yR = 244;

  F.head("DONOR", AX, 40);
  F.cap("a strong element", AX, 58, { bw: 240 });
  F.track(AX, yD, AW, [{ c: 0.5, w: 28, a: 50 }], { seed: 3.1 });
  F.dna(AX, AX + AW, yD);
  motifMark(F, AX + AW / 2, yD);      // by this slide the gold block needs no label

  F.head("RECIPIENT", AX, yR - 44);
  F.cap("a weak location", AX, yR - 26, { bw: 240 });
  F.track(AX, yR, AW, [{ c: 0.5, w: 28, a: 15 }], { seed: 5.4 });
  F.dna(AX, AX + AW, yR);
  motifMark(F, AX + AW / 2, yR, { t: 55 });

  F.line(AX + AW + 46, yD + 16, AX + AW + 46, yR - 26, { color: EDIT, w: 1.6, arrow: true });
  F.txt("write the donor's sequence", AX + AW + 66, (yD + yR) / 2 - 9, { size: 12, bold: true, color: EDIT_D, bw: 300 });
  F.txt("into the recipient's place", AX + AW + 66, (yD + yR) / 2 + 9, { size: 12, bold: true, color: EDIT_D, bw: 300 });

  F.line(602, 30, 602, 296, { dash: "dash" });
  F.head("WHICH HEIGHT DOES IT TAKE?", BX, 40, { bw: 420 });

  F.track(BX, yD + 24, BW, [{ c: 0.5, w: 26, a: 50 }], { seed: 3.1, t: 92 });
  F.dna(BX, BX + BW, yD + 24);
  motifMark(F, BX + BW / 2, yD + 24);
  F.txt("the donor's", BX + BW + 18, yD + 16, { size: 12, bold: true, color: SIG, bw: 220 });
  F.cap("the sequence wins", BX + BW + 18, yD + 33, { bw: 220 });

  F.track(BX, yR + 16, BW, [{ c: 0.5, w: 26, a: 15 }], { seed: 5.4, color: EDIT, t: 90 });
  F.dna(BX, BX + BW, yR + 16);
  motifMark(F, BX + BW / 2, yR + 16);
  F.txt("the recipient's", BX + BW + 18, yR + 8, { size: 12, bold: true, color: EDIT_D, bw: 220 });
  F.cap("the place wins", BX + BW + 18, yR + 25, { bw: 220 });

  F.txt("?", BX + BW / 2, (yD + yR) / 2 + 34, { anchor: "middle", size: 30, bold: true, color: HAIR, bw: 70 });

  F.line(14, 310, 1086, 310, { dash: "dash" });
  F.cap("62,250 swaps   ·   two widths: 170 bp and the whole peak   ·   "
    + "intervals from a cluster bootstrap over the 250 donors", 14, 332, { bw: 1000 });
  return 354;
}

// ================================================================== d4  motif shuffle
function d4Shuffle(pres, slide, box) {
  const F = frame(pres, slide, box);
  const X0 = 228, W = 552, CX = X0 + W / 2;
  const PA = CX - 134, PB = CX + 134, CTRLX = CX - 56;
  const y0 = 110;

  F.head("ONE SITE", 14, 40, { bw: 220 });
  F.cap("its own motif,", 14, 58, { bw: 220 });
  F.cap("and its partners", 14, 74, { bw: 220 });

  F.track(X0, y0, W, [{ c: 0.5, w: 26, a: 46 }, { c: (PA - X0) / W, w: 18, a: 12 }, { c: (PB - X0) / W, w: 18, a: 10 }], { seed: 6.3 });
  F.dna(X0, X0 + W, y0);
  motifMark(F, CX, y0, { label: "AP-1", lift: 64 });
  motifMark(F, PA, y0, { label: "partner", w: 22, color: PART, labelColor: PART });
  motifMark(F, PB, y0, { label: "partner", w: 22, color: PART, labelColor: PART });
  F.line(14, 150, 1086, 150, { dash: "dash" });

  const rows = [
    ["scramble AP-1", "ap1", "big drop"],
    ["scramble a partner", "partner", "small drop"],
    ["scramble an empty window", "control", "no drop"],
  ];
  rows.forEach(([lab, kind, out], i) => {
    const y = 204 + i * 62, ctrl = kind === "control";
    F.txt(lab, 14, y, { size: 12, bold: !ctrl, color: ctrl ? MUTED : INK, bw: 210 });
    const peaks = {
      ap1: [{ c: 0.5, w: 26, a: 14 }, { c: (PA - X0) / W, w: 18, a: 12 }, { c: (PB - X0) / W, w: 18, a: 10 }],
      partner: [{ c: 0.5, w: 26, a: 46 }, { c: (PA - X0) / W, w: 18, a: 2 }, { c: (PB - X0) / W, w: 18, a: 10 }],
      control: [{ c: 0.5, w: 26, a: 46 }, { c: (PA - X0) / W, w: 18, a: 12 }, { c: (PB - X0) / W, w: 18, a: 10 }],
    }[kind];
    F.track(X0, y, W, peaks, { seed: 6.3, color: ctrl ? SIG : EDIT, t: 88 });
    F.dna(X0, X0 + W, y);
    if (kind === "ap1") {
      F.scramble(CX, y, 28);
      motifMark(F, PA, y, { w: 22, color: PART }); motifMark(F, PB, y, { w: 22, color: PART });
    } else if (kind === "partner") {
      F.scramble(PA, y, 24);
      motifMark(F, CX, y); motifMark(F, PB, y, { w: 22, color: PART });
    } else {
      F.scramble(CTRLX, y, 24, { color: MUTED });
      F.block(CTRLX, y, 32, { fill: false, line: MUTED, lw: 1.1, dash: "dash", h: 17 });
      motifMark(F, CX, y);
      motifMark(F, PA, y, { w: 22, color: PART }); motifMark(F, PB, y, { w: 22, color: PART });
    }
    F.txt(out, X0 + W + 20, y, { size: 12, bold: !ctrl, color: ctrl ? MUTED : EDIT_D, bw: 150 });
  });

  const bx = X0 + W + 128;
  F.line(bx, 266, bx + 10, 266, { color: INK, w: 1.2 });
  F.line(bx + 10, 266, bx + 10, 328, { color: INK, w: 1.2 });
  F.line(bx, 328, bx + 10, 328, { color: INK, w: 1.2 });
  F.cap("a partner's effect", bx + 18, 290, { bw: 180 });
  F.txt("is the difference", bx + 18, 306, { size: 12, bold: true, bw: 180 });

  F.line(14, 362, 1086, 362, { dash: "dash" });
  F.cap("the empty window is the same width, at the mirrored distance, with no footprint in it "
    + "— so the comparison never leaves the site", 14, 384, { bw: 1000 });
  return 406;
}

// ================================================================== d5  scan
function d5Scan(pres, slide, box) {
  const F = frame(pres, slide, box);
  const X0 = 162, W = 760, CX = X0 + W / 2, y0 = 118;

  F.head("WALK A 10 bp WINDOW ACROSS THE LOCUS", 14, 40, { bw: 520 });
  F.track(X0, y0, W, [{ c: 0.5, w: 26, a: 48 }], { seed: 7.7 });
  F.dna(X0, X0 + W, y0);
  motifMark(F, CX, y0, { label: "AP-1", lift: 66 });

  const NT = 30, tw = W / NT;                      // thirty tiles stand in for sixty, one of them live
  for (let k = 0; k < NT; k++) {
    const live = k === 11;
    slide.addShape(pres.ShapeType.roundRect, {
      x: F.X(X0 + k * tw + 0.8), y: F.Y(y0 + 11), w: F.S(tw - 1.6), h: F.S(17), rectRadius: 0.01,
      fill: live ? { color: EDIT, transparency: 25 } : { type: "none" },
      line: { color: EDIT, width: live ? 1.1 : 0.75, transparency: live ? 0 : 62 },
    });
  }
  F.line(X0, y0 + 40, X0 + W, y0 + 40, {});
  F.cap("−300 bp", X0, y0 + 58, { bw: 140 });
  F.cap("+290 bp", X0 + W, y0 + 58, { anchor: "end", bw: 140 });
  F.cap("60 tiles, one prediction each", CX, y0 + 58, { anchor: "middle", bw: 300 });

  F.line(14, 212, 1086, 212, { dash: "dash" });
  F.head("THEN ZOOM IN AND CHANGE EVERY BASE", 14, 252, { bw: 520 });

  const seq = "A C G T T C C T G A G T C A C C G T".split(" ");
  const hi = 9;
  const runs = seq.map((ch, i) => ({
    text: ch + (i < seq.length - 1 ? "  " : ""),
    options: { color: i === hi ? EDIT_D : INK, bold: i === hi },
  }));
  slide.addText(runs, {
    x: F.X(CX - 230), y: F.Y(300) - F.S(16), w: F.S(460), h: F.S(32),
    fontFace: "Arial", fontSize: 13, align: "center", valign: "middle", margin: 0, isTextBox: true,
    charSpacing: 1.2,
  });
  F.block(CX - 230 + (460 * (hi + 0.5)) / seq.length, 300, 24, { color: EDIT, t: 86, lw: 1.1, h: 26 });
  F.cap("61 positions", CX - 244, 300, { anchor: "end", bw: 170 });
  F.line(CX - 230, 274, X0 + 290, y0 + 30, { dash: "dash", transparency: 30 });
  F.line(CX + 230, 274, X0 + 302, y0 + 30, { dash: "dash", transparency: 30 });

  const hx = CX - 230 + (460 * (hi + 0.5)) / seq.length;
  F.line(hx, 316, hx, 330, { color: EDIT, w: 1.2 });
  F.line(hx - 38, 330, hx + 38, 330, { color: EDIT, w: 1.1 });
  ["A", "C", "G"].forEach((alt, j) => {
    const x = hx + (j - 1) * 30;
    F.line(x, 330, x, 338, { color: EDIT, w: 1.1 });
    F.block(x, 350, 24, { color: EDIT, t: 88, lw: 1.1, h: 24 });
    F.txt(alt, x, 350, { anchor: "middle", size: 12.5, bold: true, color: EDIT_D, bw: 60 });
  });
  F.cap("each of the other three bases, one prediction each", hx + 62, 350, { bw: 420 });

  F.line(14, 392, 1086, 392, { dash: "dash" });
  F.cap("each site's own noise floor is removed before its profile is normalised, because a less accessible site "
    + "is a noisier one", 14, 414, { bw: 1000 });
  return 436;
}

// ================================================================== d6  titration
function d6Titration(pres, slide, box) {
  const F = frame(pres, slide, box);
  const y = 126, xs = [124, 336, 548, 760];
  const core = "TGAGTCA";
  const states = [
    ["TGAGTCA", "perfect core", 20],
    ["TGACTCA", "1 base changed", 42],
    ["TCACTTA", "3 bases changed", 64],
    ["CATCAGG", "core gone", null],
  ];

  F.head("WEAKEN THE MOTIF, ONE BASE AT A TIME", 14, 40, { bw: 520 });
  F.line(xs[0] + 92, y - 54, xs[3] - 92, y - 54, { color: MUTED, w: 1.1, arrow: true, transparency: 55 });
  F.cap("weaker grip", (xs[0] + xs[3]) / 2, y - 70, { anchor: "middle", bw: 200 });

  states.forEach(([alt, lab, t], i) => {
    const x = xs[i];
    F.dna(x - 76, x + 76, y);
    if (t === null) F.scramble(x, y, 48);
    else F.block(x, y, 48, { color: GOLD, t, lw: 1.1 });
    const runs = alt.split("").map((ch, j) => ({
      text: ch + (j < 6 ? " " : ""),
      options: { color: ch !== core[j] ? EDIT_D : INK, bold: ch !== core[j] },
    }));
    slide.addText(runs, {
      x: F.X(x - 80), y: F.Y(y + 28) - F.S(14), w: F.S(160), h: F.S(28),
      fontFace: "Arial", fontSize: 12.5, align: "center", valign: "middle", margin: 0, isTextBox: true,
      charSpacing: 1.5,
    });
    F.cap(lab, x, y + 54, { anchor: "middle", bw: 220 });
  });
  F.cap("also run", 912, y - 18, { bw: 200 });
  ["·  upgrade the core", "·  half site", "·  reverse complement"].forEach((t, i) =>
    F.txt(t, 912, y + 4 + i * 19, { size: 12, bw: 230 }));

  F.line(14, 214, 1086, 214, { dash: "dash" });
  F.head("OR ADD A SECOND COPY AT A SET DISTANCE", 14, 250, { bw: 560 });

  const BX = 236, BW = 520, yy = 296;
  F.dna(BX, BX + BW, yy);
  motifMark(F, BX + 150, yy);
  motifMark(F, BX + 350, yy);
  F.line(BX + 150, yy + 24, BX + 350, yy + 24, { color: MUTED, w: 1.1 });
  F.line(BX + 150, yy + 19, BX + 150, yy + 29, { color: MUTED, w: 1.1 });
  F.line(BX + 350, yy + 19, BX + 350, yy + 29, { color: MUTED, w: 1.1 });
  F.txt("d", BX + 250, yy + 18, { anchor: "middle", size: 12, bold: true, color: MUTED, bw: 60 });
  F.txt("a second AP-1", BX + BW + 22, yy, { size: 12, bold: true, color: GOLD_D, bw: 240 });

  F.dna(BX, BX + BW, yy + 74);
  motifMark(F, BX + 150, yy + 74);
  F.scramble(BX + 350, yy + 74, 26, { color: MUTED });
  F.block(BX + 350, yy + 74, 34, { fill: false, line: MUTED, lw: 1.1, dash: "dash", h: 18 });
  F.cap("control", 14, yy + 74, { bw: 150 });
  F.cap("or the same width of scramble", BX + BW + 22, yy + 74, { bw: 300 });
  F.cap("d = 20 to 150 bp, either side", BX + 250, yy + 108, { anchor: "middle", bw: 320 });
  return 424;
}

// ================================================================== d0  what AlphaGenome is
// Three columns: what the model is given, what it is, what it returns. The point the whole deck rests on is the
// footer - the prediction is a function of sequence alone, which is what makes an edit an experiment.
function d0Model(pres, slide, box) {
  const F = frame(pres, slide, box);
  const AX = 14, AW = 286;            // 1  the input
  const BX = 352, BW = 244;           // 2  the network
  const CX = 650, CW = 258;           // 3  the predictions

  // ---- 1. a window of sequence, and nothing else -----------------------------------------------------------
  F.head("1   WHAT GOES IN", AX, 36);
  const dnaY = 148, aMid = AX + AW / 2;
  F.cap("131,072 bp of mouse genome", aMid, 108, { anchor: "middle", bw: 300 });
  F.cap("one AP-1 site, held in the middle", aMid, 124, { anchor: "middle", bw: 300 });
  F.dna(AX + 6, AX + AW - 6, dnaY);
  F.line(AX + 6, dnaY - 8, AX + 6, dnaY + 8, { color: INK, w: 1.3, transparency: 45 });
  F.line(AX + AW - 6, dnaY - 8, AX + AW - 6, dnaY + 8, { color: INK, w: 1.3, transparency: 45 });
  F.block(aMid, dnaY, 16, { h: 13, color: GOLD, t: 15, lw: 1.1 });

  // a zoom on that gold block, so "sequence" is something on the slide rather than a word
  const zx = AX + 26, zw = 234, zy = 238, zTop = 216;
  F.line(aMid - 9, dnaY + 5, zx, zTop, { color: HAIR, w: 0.9, dash: "dash" });
  F.line(aMid + 9, dnaY + 5, zx + zw, zTop, { color: HAIR, w: 0.9, dash: "dash" });
  F.block(zx + zw / 2, zy, zw, { h: 44, color: "FFFFFF", t: 0, line: HAIR, lw: 1 });
  F.block(zx + zw / 2, zy, 88, { h: 24, color: GOLD, t: 70, line: GOLD, lw: 0.6 });
  F.txt("ccaTGAGTCAgtt", zx + zw / 2, zy, { anchor: "middle", size: 13, bold: true, bw: 230,
        face: "Consolas", ls: 1.4 });
  F.cap("AP-1 core", zx + zw / 2, 276, { anchor: "middle", bw: 160, color: GOLD_D, size: 11 });
  F.cap("one-hot A C G T, and nothing else:", AX, 302, { bw: 300 });
  F.cap("no ATAC, no cell type, no annotation", AX, 318, { bw: 300 });

  // ---- 2. the network --------------------------------------------------------------------------------------
  F.head("2   THE MODEL", BX, 36);
  slide.addShape(pres.ShapeType.roundRect, {
    x: F.X(BX), y: F.Y(60), w: F.S(BW), h: F.S(238), rectRadius: 0.05,
    fill: { color: INK, transparency: 95 }, line: { color: INK, width: 1, transparency: 68 },
  });
  F.txt("AlphaGenome", BX + BW / 2, 86, { anchor: "middle", size: 14, bold: true, bw: 230 });
  F.cap("sequence in, signal out", BX + BW / 2, 106, { anchor: "middle", bw: 240 });

  const bands = [
    [152, "convolutions", "motifs, and how they sit together"],
    [208, "transformer", "all 131 kb attended to at once"],
    [264, "decoders and heads", "back up to 1 bp, one per track"],
  ];
  for (const [y, t1, t2] of bands) {
    F.block(BX + BW / 2, y, BW - 44, { h: 44, color: SIG, t: 90, line: SIG, lw: 0.9 });
    F.txt(t1, BX + BW / 2, y - 10, { anchor: "middle", size: 12, bold: true, bw: 210 });
    F.cap(t2, BX + BW / 2, y + 9, { anchor: "middle", bw: 210, size: 10.5 });
  }
  F.line(BX + BW / 2, 175, BX + BW / 2, 186, { color: SIG, w: 1.2, arrow: true });
  F.line(BX + BW / 2, 231, BX + BW / 2, 242, { color: SIG, w: 1.2, arrow: true });
  F.cap("trained on thousands of genomic tracks", BX, 318, { bw: 300 });

  // ---- 3. the predictions ----------------------------------------------------------------------------------
  F.head("3   WHAT COMES OUT", CX, 36);
  const LX = CX + CW + 14;
  const rows = [
    [112, "ATAC  ·  fibroblast", "the one we read", SIG, [{ c: 0.42, w: 24, a: 34 }, { c: 0.68, w: 15, a: 17 }]],
    [174, "ATAC  ·  cardiomyocyte", null, HAIR, [{ c: 0.44, w: 22, a: 21 }, { c: 0.70, w: 14, a: 12 }]],
    [236, "RNA-seq", "exon by exon", HAIR, [{ c: 0.22, w: 7, a: 25 }, { c: 0.37, w: 5, a: 30 }, { c: 0.63, w: 9, a: 20 }]],
    [298, "ChIP-seq  ·  CAGE", "and more", HAIR, [{ c: 0.33, w: 12, a: 22 }, { c: 0.58, w: 18, a: 16 }]],
  ];
  rows.forEach(([y, lab, sub, col, peaks], i) => {
    F.track(CX, y, CW, peaks, { seed: 1.4 + i * 1.7, hmax: 42, color: col, t: i === 0 ? 84 : 91 });
    F.dna(CX, CX + CW, y);
    F.txt(lab, LX, y - 27, { size: 11.5, bold: i === 0, bw: 164, color: i === 0 ? INK : MUTED });
    F.cap(sub || "", LX, y - 10, { bw: 164, size: 10.5, color: i === 0 ? GOLD_D : HAIR });
  });
  F.line(CX - 8, 70, CX - 8, 112, { color: GOLD_D, w: 3 });
  F.cap("1 bp resolution, right across the window", CX, 318, { bw: 320 });
  F.cap("plus contact maps and splicing", CX, 334, { bw: 320 });

  // ---- the flow, and the sentence the deck rests on ---------------------------------------------------------
  F.line(AX + AW + 8, 180, BX - 8, 180, { color: INK, w: 1.4, arrow: true });
  F.line(BX + BW + 8, 180, CX - 16, 180, { color: INK, w: 1.4, arrow: true });

  F.line(14, 358, 1086, 358, { dash: "dash" });
  F.txt("The prediction is a function of the sequence alone: change a base, predict again, and the difference "
      + "is caused by the change.", 14, 380, { size: 12, bold: true, bw: 1060 });
  return 400;
}

// the height each drawing occupies, in design units, so the caller can centre it
d0Model.designHeight = 400;
d1Pipeline.designHeight = 334;
d2Transplant.designHeight = 362;
d3Swap.designHeight = 354;
d4Shuffle.designHeight = 406;
d5Scan.designHeight = 436;
d6Titration.designHeight = 424;

module.exports = {
  d0Model, d1Pipeline, d2Transplant, d3Swap, d4Shuffle, d5Scan, d6Titration,
  COLORS: { INK, MUTED, HAIR, SIG, GOLD, GOLD_D, EDIT, EDIT_D, PART },
};
