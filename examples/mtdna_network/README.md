# Mitochondrial network and early drift

The four-species model and the variance calculation are introduced in
[chapter 20](../../docs/guide/20-mtdna-network.md).

```sh
python examples/mtdna_network/study.py --size quick
python examples/mtdna_network/study.py --size full --threads 4
python examples/mtdna_network/study.py --figures docs/guide/img/mtdna_network
```

The condition table is extracted by `fetch.py`. Add
`--online-summaries` to also regenerate the source variance series.
The fresh slope normally uses the first 2 days, while the archived
slope uses the first 100 days. Figure 2H varies the ratio of fusion to
fission at fixed magnitude. Figure S3D varies the magnitude at fixed
ratio. Both sweeps use the published initial counts. Selected extreme
switching magnitudes use a shorter reported window and 1,024 paths.
All other full-size conditions use 2,000 paths, with checkpoints.

The three primary intervals for simulated slope divided by the
Equation 13 slope are 0.904–1.074, 0.903–1.062 and 0.935–1.084.
Some other source conditions have wider intervals. Read the
`within_band` column alongside the interval and the number of paths.
[Paper and supplement](https://doi.org/10.1534/genetics.119.302423)
