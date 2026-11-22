# Saved C6 tables

Both tables are drawn by `draw(tables, folder)` in
`publications/stimulus_fit.py` into `img/C6-distributions.png` and
`img/C6-criteria.png`.

## models.csv

One row per fitted stimulus hypothesis. `k` counts the four basal rates
plus one stimulated value per changed rate. `nll` is the total joint
negative log-likelihood in nats over all 376 cells, scored with an RNA cutoff
of 450; `search_nll` is the same vector at the search cutoff of 300. AIC is
`2k + 2 NLL`, BIC is `k ln 376 + 2 NLL`, AICc adds `2k(k+1)/(376-k-1)`;
the deltas are taken against the best delivered row. `retries` lists the
best NLL of each DE retry and `polish` the best of the CMA-ES refinement;
the rates are those of the reported vector before and after the stimulus.

- `run=delivered`: recounted from the implementer's study run (8 threads,
  not rerun here): `run.py C6 --preset study --retries 4 --evaluations 550`
  with the default 300 CMA-ES calls, DE seeds 4497 to 4500, CMA-ES seed
  44497. `nll` was rescored in review from the saved rates and agrees to
  1e-9 nats; `retries` and `polish` are the minima of the saved traces.
- `run=review-cma-1500`: rerun in review. CMA-ES with sigma 0.3 in log10
  coordinates, population 12, 1,500 calls, seed 7001, started from the
  delivered `kon-koff` vector with the active rate unchanged by the stimulus.
- `run=review-de-2000`: rerun in review. One DE retry of 2,000 calls,
  population 16, seed 7002, from the full box; NLL at cutoff 300 only.

## distributions.csv

For each time, number of active sites (0 to 2) and RNA count (0 to 300):
`predicted` is the joint probability of the delivered `kon-koff` fit at
cutoff 450 (mass above 300 copies is at most 1.8e-05 at 25 minutes), and
`observed` the number of archive cells in that bin. The cell counts are
derived from the BayFish data files of Gómez-Schiavon et al. (2017),
doi:10.1186/s13059-017-1297-9, at
`70d1b685da4b4da8aa8a64934e200a5e3593791d` (GPL-3.0); see
`THIRD_PARTY_NOTICES.md`.
