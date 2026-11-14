# How do I ...?

Fifty common tasks, one per line, each with links to the sections of the guide that answer it.

## Getting started

- How do I install CRN-craft and check that it works? [1, Install the Python package](guide/01-install.md#install-the-python-package); [Check the installation](guide/01-install.md#check-the-installation)
- How does CRN-craft compare with GillesPy2 or rebop? [How CRN-craft compares](comparison.md)
- How do I write a reaction network in Python? [2, Species, parameters and reactions](guide/02-first-model.md#species-parameters-and-reactions)
- How do I write a rate that is not mass action, such as Michaelis–Menten or Hill? [2, Rates as numbers and as strings](guide/02-first-model.md#rates-as-numbers-and-as-strings)
- How do I get the rate constant right for `A + A → B`? [2, Rates as numbers and as strings](guide/02-first-model.md#rates-as-numbers-and-as-strings); [8, Reproducing a published simulation](guide/08-reproducibility.md#reproducing-a-published-simulation), step 1
- How do I check a model without running it, and read the error? [2, Checking a model](guide/02-first-model.md#checking-a-model); [Diagnostics](reference/diagnostics.md)

## Running

- How do I run 1,000 replicates on 8 threads? [3, Many trajectories in one call](guide/03-running.md#many-trajectories-in-one-call); [Seeds and threads](guide/03-running.md#seeds-and-threads); [5, Threads](guide/05-speed.md#threads)
- How do I get the same result again? [3, Seeds and threads](guide/03-running.md#seeds-and-threads)
- How do I see whether a run finished, and how do I limit its events or its wall time? [3, Check the termination](guide/03-running.md#check-the-termination); [Budgets](guide/03-running.md#budgets); [Troubleshooting, Runs that end early](guide/troubleshooting.md#runs-that-end-early)
- How do I stop a long run and keep what is done? [3, Stopping a run](guide/03-running.md#stopping-a-run); [10, Stop a run with a cancel file](guide/10-command-line.md#stop-a-run-with-a-cancel-file)

## Results

- How do I get the mean and the variance of a species over time, as a data frame or a NumPy array? [4, What a result holds](guide/04-results.md#what-a-result-holds); [Select by name](guide/04-results.md#select-by-name); [Mean, variance and standard error](guide/04-results.md#mean-variance-and-standard-error)
- How do I get the median or a 95% band? [4, Quantiles and the median](guide/04-results.md#quantiles-and-the-median)
- How do I keep single trajectories and read them? [4, Single trajectories](guide/04-results.md#single-trajectories)
- How do I export results for plotting, or plot them in Python? [4, Export and plot](guide/04-results.md#export-and-plot)
- How do I save a run with its model and settings in one file? [4, Save a run in one file](guide/04-results.md#save-a-run-in-one-file)

## Speed

- How do I run one model for many parameter values? [5, Parameter sweeps](guide/05-speed.md#parameter-sweeps)
- How do I compile a model to make it faster? [5, Compiled models](guide/05-speed.md#compiled-models)
- My model has hundreds of reactions. Which selector is fastest? [7, Static Direct and its selectors](guide/07-methods.md#static-direct-and-its-selectors); [Larger models](guide/07-methods.md#larger-models)

## Events

- How do I stop a replicate when a condition holds, such as extinction or a threshold? [6, Stop conditions](guide/06-stops-trackers-logs.md#stop-conditions)
- How do I get the exact time at which a threshold is first reached, the peak of a species, or how often it appears and dies out? [6, Trackers](guide/06-stops-trackers-logs.md#trackers); [11, First passage from saved paths](guide/11-analysis-import-export.md#first-passage-from-saved-paths)
- How do I record the single events of a run? [6, Event logs](guide/06-stops-trackers-logs.md#event-logs)
- How do I run my own Python code at every event? [6, Callbacks](guide/06-stops-trackers-logs.md#callbacks)

## Methods

- Which method should I use? [7, Which method for which model](guide/07-methods.md#which-method-for-which-model)
- How do I get the exact distribution instead of samples? [7, Exact distributions with finite state projection](guide/07-methods.md#exact-distributions-with-finite-state-projection)
- When can I trust tau leaping? [7, Approximate leaping](guide/07-methods.md#approximate-leaping)

## Reproducing

- How do I reproduce a simulation from a published paper? [8, Reproducing a published simulation](guide/08-reproducibility.md#reproducing-a-published-simulation)
- The paper gives its seed. Can I get the same trajectories? No: [8, Reproducing a published simulation](guide/08-reproducibility.md#reproducing-a-published-simulation), step 2
- How do I compare my ensemble with published means and variances, and how many replicates do I need? [8, Reproducing a published simulation](guide/08-reproducibility.md#reproducing-a-published-simulation), steps 3 and 4; [4, Mean, variance and standard error](guide/04-results.md#mean-variance-and-standard-error)
- The paper used a fixed time step. How do I match it? [7, Approximate leaping](guide/07-methods.md#approximate-leaping); [12d, Matching an hourly scheme](guide/12d-mtdna-methods.md#matching-an-hourly-scheme)
- What do I report so that others can reproduce my run? [8, What to record for a paper](guide/08-reproducibility.md#what-to-record-for-a-paper)
- Why do I get other numbers with another method, selector, thread count or release? [7, Same law, different random numbers](guide/07-methods.md#same-law-different-random-numbers); [8, What does not reproduce bit for bit](guide/08-reproducibility.md#what-does-not-reproduce-bit-for-bit)

## More kinds of models, and fitting

- How do I couple a continuous quantity, such as ATP or a drug level, to the reactions? [9, A variable between events](guide/09-continuous-variables.md#a-variable-between-events)
- How do I add a time-dependent forcing, such as a daily rhythm or a decaying dose? [9, Independent variables](guide/09-continuous-variables.md#independent-variables)
- How do I fit parameters to data? [13, Fitting parameters to data](guide/13-fitting.md); [the optimization tutorials](../tutorials/optimization/README.md)

## Command line, files and other formats

- How do I run a model from the command line or on a cluster? [10, A request file](guide/10-command-line.md#a-request-file); [Run, archive and verify](guide/10-command-line.md#run-archive-and-verify)
- How do I archive a run so that it can be verified later? [10, Run, archive and verify](guide/10-command-line.md#run-archive-and-verify); [Replay](guide/10-command-line.md#replay)
- How do I import an SBML or BioNetGen model? [11, Import SBML and BioNetGen](guide/11-analysis-import-export.md#import-sbml-and-bionetgen)
- How do I convert a model file of release 0.6.0? [11, Convert a network model of release 0.6.0](guide/11-analysis-import-export.md#convert-a-network-model-of-release-060); [Migrating from earlier releases](guide/migrating.md)
- How do I find the conservation laws of a network? [11, Conservation laws and structure](guide/11-analysis-import-export.md#conservation-laws-and-structure)
- How do I check a model against another simulator? [11, Export to GillesPy2](guide/11-analysis-import-export.md#export-to-gillespy2); [8, Reproducing a published simulation](guide/08-reproducibility.md#reproducing-a-published-simulation)
- How do I load a model from a JSON file, or save one? [14, Save and load a model file](guide/14-model-files.md#save-and-load-a-model-file)
- How do I model molecules that carry a property, such as a clone, a size or an age? [14, A model document for three strains](guide/14-model-files.md#a-model-document-for-three-strains)
- How do I give every new mutant its own identity? [14, New identities and random attributes](guide/14-model-files.md#new-identities-and-random-attributes); [14a, Every lineage with the clone model](guide/14a-mtdna-model-files.md#every-lineage-with-the-clone-model)

## mtDNA models

- How do I simulate mutant takeover in single cells? [12, The model](guide/12-mtdna.md#the-model); [12a, A first ensemble](guide/12a-mtdna-takeover.md#a-first-ensemble); [Takeover fraction and time to threshold](guide/12a-mtdna-takeover.md#takeover-fraction-and-time-to-threshold)
- How do I change copy number, half-life, initial heteroplasmy or the threshold? [12a, Variants](guide/12a-mtdna-takeover.md#variants)
- How do I follow mutant classes, and single lineages? [12b, Mutant classes](guide/12b-mtdna-events.md); [14a, Every lineage with the clone model](guide/14a-mtdna-model-files.md#every-lineage-with-the-clone-model)
- How do I let the population size vary with ATP? [12c, Variable population: ATP](guide/12c-mtdna-atp.md)
- How do I use the fast path `mtdna_direct`, and why is it not faster on my model? [12d, The fast path](guide/12d-mtdna-methods.md#the-fast-path); [14a, The fast path](guide/14a-mtdna-model-files.md#the-fast-path); [Troubleshooting, Methods and models](guide/troubleshooting.md#methods-and-models)

## Other languages and older code

- How do I use CRN-craft from R? [15, Using CRN-craft from R](guide/15-from-r.md)
- My script was written for 0.6.0 (`crn_craft.rules`, `run_spec`, `crn-craft ssa`). What do I change? [Migrating from earlier releases](guide/migrating.md)
