"""Reproduce every figure, table and statistic of the paper in one go.

    python run_all.py

Order matters: the statistics scripts write the CSVs that Figure 3 and the LaTeX tables read.
"""
import subprocess
import sys
import time
from pathlib import Path

ROOT = Path(__file__).resolve().parent
STEPS = [
    ("Table 1 / S2  dataset counts", "tables_and_statistics/dataset_summary.py"),
    ("Tables 3,5,7 + DeLong/Holm", "tables_and_statistics/compute_test_metrics.py"),
    ("Tables 4, S1, S3, S4 inputs", "tables_and_statistics/compute_extra_metrics.py"),
    ("LaTeX tables", "tables_and_statistics/make_latex_tables.py"),
    ("Figure 1  workflow", "fig1_workflow/make_fig1_workflow.py"),
    ("Figures 2 & 7  ROC curves", "fig2_fig7_roc_curves/make_roc_curves.py"),
    ("Figure 3  Sn vs Sp", "fig3_sensitivity_specificity/make_fig3_sn_sp.py"),
    ("Figures 4-6, 8-10  logos + Table 6", "fig4_to_fig10_sequence_logos/make_logos.py"),
]
for name, script in STEPS:
    t = time.time(); print(f"\n=== {name}  ({script})")
    subprocess.run([sys.executable, str(ROOT / script)], check=True)
    print(f"    done in {time.time() - t:.0f}s")
print("\nAll outputs are in outputs/figures, outputs/results and outputs/latex_tables")
