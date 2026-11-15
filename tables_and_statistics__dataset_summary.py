"""Record counts per dataset and split -> Table 1 (human) and Table S2 (mouse)."""
import sys
from pathlib import Path

import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import config
from common.data_io import find_fasta, read_fasta, species_dir


def main():
    rows = []
    for sp, mods in [("human", config.TABLE_ORDER), ("mouse", config.TABLE_ORDER)]:
        for mod in mods:
            _, S, y = read_fasta(find_fasta(sp, mod))
            split = pd.read_csv(species_dir(sp) / "splits" / f"{mod}_split.csv").split.value_counts()
            rows.append(dict(species=sp, mod=mod, total=len(y), positives=int(y.sum()),
                             negatives=int((y == 0).sum()), length=sorted({len(s) for s in S}),
                             train=int(split.get("train", 0)), val=int(split.get("val", 0)),
                             test=int(split.get("test", 0))))
    R = pd.DataFrame(rows); R.to_csv(config.OUT_RES / "dataset_summary.csv", index=False)
    print(R.to_string(index=False))


if __name__ == "__main__":
    main()
