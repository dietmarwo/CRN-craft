# C8 saved results

Tables of [C8](../../C8-experiment-design.md), rerun in review on one
thread with

```bash
python tutorials/optimization/run.py C8 --preset smoke --threads 1 \
    --native-pair-samples 100 --expanded-bounds 100 135 6200 \
    --data-dir DATA --output OUT
```

(no random seed is involved in these tables: they are fixed-design
solves; about 10 minutes wall time). The figures in `../../img/C8-*.png`
are drawn from these files only, by `draw(tables, folder)` in
`publications/experiment_design.py`.

| File | Content |
|---|---|
| `burst-sweep.csv` | Bursting gene, Figure 5: `grid` (`paper-60`: the article's 60 intervals; `script-50`: the author script's 50; `fine`: 141 intervals from 3 to 12 min), `interval_min`, `native_eopt` (ours: smallest eigenvalue, 1000 cells, natural parameters koff and transcription rate), `archive_eopt` (the author's value for the script grid, from `bursting_gene/out/E_opt_fsamp_fsp1.0_100.txt`; dimensionless units, ours / 160000) |
| `native-pair.csv` | Toggle, Figure 7, designs A = (2 h, 3 h, UV 9) and B = (3 h, 4 h, UV 6), parameter samples 0 to 99: for each design `native` (ours) and `archive` smallest eigenvalue, `sink` (largest sink probability of our solve on the archive's domain) and `matrix_gap` (relative Frobenius difference of our 7 × 7 matrix from the archived one) |
| `domain-expansion.csv` | The same columns for the ten samples whose sink exceeds 1e-3, solved on the domain `x ≤ 100, y ≤ 135, (x − 4)(y − 4) < 6200`; `archive` repeats the archived value (original domain) |

The archived values are derived from the author archive
[FSPFIM_2018 at commit 2e2f76f](https://github.com/Munsky/FSPFIM_2018/tree/2e2f76fdeeb0d9290b48aa981210422d254b07c2)
(see `THIRD_PARTY_NOTICES.md`). No study run (8 threads, DE search) is
saved here: no study run with the conventions of these tables exists.
