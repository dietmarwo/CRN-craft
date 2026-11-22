# C10 saved results

| File | Origin |
| --- | --- |
| `runs.csv` | implementer's study runs on 8 threads (`run.py C10 --preset study --evaluations 2000`, then `--evaluations 20000 --broad-only --retries 1 --resume-from`), seeds 51941 + retry; recounted from their `summary.json`, not rerun (the first 240 evaluations of both 20,000-evaluation runs were replayed in review and are identical) |
| `best-so-far.csv` | the same runs' traces, best cost so far at about 80 log-spaced evaluation counts (failed evaluations excluded); rows with method `least-squares-…` are the review rerun below |
| `local.csv` | rerun in review on one thread: `publications.pathway.local_fit` at the search tolerance, spreads 0.3, 1 and 1 decade, seeds 7, 8, 9 (the study preset runs the same fits) |
| `nominal.csv` | `publications.pathway.reference` at the nominal vector, all 16 conditions and 21 times (identical to the implementer's `reconstruction.json`) |

The figures `img/C10-*.png` are drawn from these tables by
`publications.pathway.draw`.
