# Slide deck

`build_deck.js` writes `Analysis/AP1_AlphaGenome_results.pptx`: a 12-slide black-on-white summary of
the AlphaGenome AP-1 perturbation results, with every figure pulled from `Analysis/` so the deck
rebuilds whenever the analyses are rerun.

    npm install pptxgenjs                     # once, in this directory
    NODE_PATH="$PWD/node_modules" node build_deck.js

The second command needs `NODE_PATH` because the skill's `apply_theme.js` (which writes the theme
colours into the packed file) resolves `pptxgenjs` from its own location, not from this one.

Speaker notes carry the statistics and caveats behind each slide; `markitdown` on the .pptx prints
all the text, and LibreOffice converts it to PDF for checking.
