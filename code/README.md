# PACE code and reproducibility tools

- [`original_PACE/`](original_PACE/) — 14 unmodified author-supplied Python scripts and original `RUN_ORDER.md`.
- [`../S12/DATA/`](../S12/DATA/) — scripts and notebooks used in subsequent independent analyses and validations.
- [`verify_checksums.py`](verify_checksums.py) — validates scientific data/code against `MASTER_MANIFEST_SHA256.csv`.
- [`verify_repository_manifest.py`](verify_repository_manifest.py) — validates additional published repository documentation.
- [`reassemble_archives.py`](reassemble_archives.py) — verifies and rejoins losslessly split ESM-2 and GA/NSGA-II result archives.

```bash
python code/verify_checksums.py
python code/verify_repository_manifest.py
python code/reassemble_archives.py --verify-only
python code/reassemble_archives.py
```

Source/model reproducibility is conditional on specific data and compatible dependencies. Original fitted model binaries and exact historical GA search state were not supplied; these scripts and supplementary predictions do not claim otherwise.
