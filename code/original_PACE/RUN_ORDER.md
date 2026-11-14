# PACE pipeline — reproduction order

Environment: Python 3.12; numpy, pandas, scikit-learn, scipy, biopython, shap, matplotlib, openpyxl, Pillow.

Download the three public archives first:
  https://codeload.github.com/bcgsc/AMPlify/tar.gz/refs/heads/master
  https://codeload.github.com/tlawrence3/amPEPpy/tar.gz/refs/heads/master
  https://codeload.github.com/zswitten/Antimicrobial-Peptides/tar.gz/refs/heads/master
Extract into ./data/ preserving the top-level folder names.

Then run:
  python 01_curate.py                       # curation + audit log
  python 02b_classify_one.py {LR,SVM,RF,ET,GB,MLP}   # one call per algorithm
  python 02c_ensemble.py                    # soft-voting ensemble + production refit
  python 03b_regress_one.py {0..5}          # one call per ESKAPE target
  python 03b_regress_one.py hem             # haemolysis model
  python 04_analysis.py                     # SHAP, permutation, similarity bins, correlations
  python 05_design_ga.py                    # genetic algorithm
  python 06_leads.py                        # diversity selection, reference checks
  python 07a_fig1.py; python 07_fig1to3.py; python 08_fig4to6.py; python 09_fig7to8.py
  python 10_supplementary.py                # supplementary workbook

All random seeds are fixed (2024 for models, 7 for the genetic algorithm).
Single-core runtime is roughly 45-60 minutes end to end.
