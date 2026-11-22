# C6 · Select a stimulus response from joint counts

The observed data contain 376 cells from the four archive files.
Each row supplies one mature RNA count and the number of nonzero active
sites among two intensity measurements. The first snapshot uses the
unstimulated stationary law. The remaining snapshots use the same initial
distribution propagated at the stimulated rates.
[Paper and Figure 5](https://doi.org/10.1186/s13059-017-1297-9) · [Pinned BayFish archive](https://github.com/mgschiavon/BayFish/tree/70d1b685da4b4da8aa8a64934e200a5e3593791d)

```bash
python tutorials/optimization/publications/sources.py C6
python tutorials/optimization/run.py C6 --preset study --output results/c6
```

Compare the seven models on total joint negative log-likelihood. Keep the
published RNA decay fixed and count each independently adjustable rate in
AIC and BIC. Stationary truncated solves use the reflecting boundary at 300
RNA copies. The time-dependent solve records escaping probability, and
finalists are scored again through 450 copies. A normalized balance solve
handles candidates for which a unit anchor gives invalid probabilities.
The study also evaluates nested candidates and refines them with CMA-ES.

`joint-distributions.png` plots maximum-likelihood predictions. The paper
also used Bayesian posteriors for DIC and posterior-averaged predictions.
