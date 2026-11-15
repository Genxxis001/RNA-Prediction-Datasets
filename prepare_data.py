"""Put the raw data into data/human and data/mouse.

Usage:
  python prepare_data.py --human-zip Human_RNA_30_prediction_files.zip \
                         --human-fasta-dir <folder with the 6 human_hg38_*.fasta> \
                         --mouse-zip Mouse_RNA_6_modifications_complete.zip
"""
import argparse
import glob
import os
import re
import shutil
import zipfile

import config


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--human-zip", required=True)
    ap.add_argument("--human-fasta-dir", required=True)
    ap.add_argument("--mouse-zip", required=True)
    a = ap.parse_args()
    for z, dest in [(a.human_zip, config.DATA_HUMAN), (a.mouse_zip, config.DATA_MOUSE)]:
        dest.mkdir(parents=True, exist_ok=True)
        with zipfile.ZipFile(z) as f:
            f.extractall(dest)
    fdir = config.DATA_HUMAN / "fasta"; fdir.mkdir(exist_ok=True)
    n = 0
    for f in glob.glob(os.path.join(a.human_fasta_dir, "*human_hg38_*_balanced_*.fasta")):
        name = re.sub(r"^.*?(human_hg38_)", r"\1", os.path.basename(f))  # drop upload prefixes
        shutil.copy(f, fdir / name); n += 1
    print(f"human FASTA copied: {n}  ->  {fdir}")
    print("mouse FASTA:", len(glob.glob(str(config.DATA_MOUSE / 'fasta' / '*.fasta'))))


if __name__ == "__main__":
    main()
