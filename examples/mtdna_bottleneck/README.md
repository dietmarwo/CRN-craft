# Birth, death and partition

See [chapter 21](../../docs/guide/21-mtdna-bottleneck.md) for the
observables and model. The first cycling phase has 29 divisions
7 hours apart, followed by seven divisions 16 hours apart.
Observation at a division uses the daughter counts. Table 1
rounds some early times, so labeled division times use the
corresponding 7-hour boundaries. Starred oocytes are plotted at 100 dpc.

```sh
python examples/mtdna_bottleneck/study.py --size quick
python examples/mtdna_bottleneck/study.py --size full
python examples/mtdna_bottleneck/study.py --figures docs/guide/img/mtdna_bottleneck
```

The inference uses the appendix's natural-log copy-number residual,
weight 1,000 on normalized heteroplasmy variance, and a hard distance
threshold of 40. Copy means are exact for the model. The fitted variance
uses the first-order expansion in fluctuations of the two counts.
Use Numba for the larger campaigns. The six phases include a terminal
state with zero replication and degradation, combining the final
constraints stated in the main text and appendix.

After fitting Table 1, inherit the calibrated parameters at the ages
of the independent HB groups. Exactly sample one-age founder-family
descendants conditional on each group's mean mutant fraction. The
validation table gives both the population statistic's interval and
the wider interval for a sample of that group's size. The supplied
`holdout` column identifies the prespecified litters highlighted in the
comparison. All HB samples come from the independent data file.

`fetch.py` parses the paper XML and Fig. 4 source file. `digitize.py`
converts checked display coordinates for the reference Figure 3 curves.
`abc.csv` records acceptance and effective sample sizes. The result
directory holds per-age bands and HB comparisons.
[Paper and supplement](https://doi.org/10.7554/eLife.07464)
