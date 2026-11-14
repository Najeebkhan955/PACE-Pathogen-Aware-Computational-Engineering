# PACE manuscript figure archive

This directory contains an **extraction script** for the figures embedded in the author's revised PACE manuscript. The script extracts the original image bytes from the provided Word file, preserves image content without resampling, and produces a SHA-256 manifest and caption mapping.

## Recreate the figure files

```bash
python -m pip install python-docx
python figures/manuscript/extract_pace_figures.py Manuscript.docx --out extracted_figures
```

Expected manuscript figures: **PACE_Figure_01.png through PACE_Figure_09.png**. Figure 9 is the matched-cohort recognition ROC-AUC comparison. The exact image bytes depend on the provided manuscript version; keep the original manuscript alongside the generated archive for provenance.

## Scientific provenance

These are **previously embedded research figures**, not independently regenerated computational plots. Extraction does not establish that the plotted data were reproduced. To regenerate individual graphs from raw data, use the source-level scripts and relevant S01–S12 research datasets, checking parameters and versions carefully.

The nine original PNG figures and their manifest were also packaged separately for the author. **Do not treat this folder as containing those PNGs until the binary files are uploaded.**

No existing repository files were replaced.
