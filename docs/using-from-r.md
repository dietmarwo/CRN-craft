# Using CRN-craft from R

CRN-craft has no native R package. From R, call the Python package through
[reticulate](https://rstudio.github.io/reticulate/): models, ensembles and
results work as in Python, and pandas results arrive as R data.frames.
The example [`examples/r/sir.R`](../examples/r/sir.R) was tested with:

- R 4.3.3
- reticulate 1.47.0
- crn-craft 0.6.0
- Python 3.13

The example imports `crn_craft.rules`, which release 0.7.0 keeps unchanged:
it writes the same documents and prints the same numbers as 0.6.0. The
top-level module `crn_craft` has the same calls under the same names, for
example `crn <- import("crn_craft")` and `crn$Model("sir")`. It differs in
what it writes: `crn$ensemble` writes the 1.0 contracts, asks for the
quantiles 0.025, 0.5 and 0.975 (so `result$statistics` also holds `median` and
`quantiles`) and leaves the Static Direct selector to `auto`, so its paths
differ from those of `crn_craft.rules`. The example has not been run again in R
on 0.7.0.

## One-time setup

```r
install.packages("reticulate")
reticulate::virtualenv_create("crn-craft", packages = "crn-craft")
```

Use a current reticulate from CRAN. Ubuntu's apt package
(`r-cran-reticulate` 1.35) left results as Python objects with NumPy 2 and
pandas 3; version 1.47 converts them to data.frames and arrays.

## Example: SIR model

```r
library(reticulate)
use_virtualenv("crn-craft", required = TRUE)
rules <- import("crn_craft.rules")

# Counts must be R integers (999L); rates and parameters are doubles.
sir <- rules$Model("sir")
invisible(sir$add_species(list(S = 999L, I = 1L, R = 0L)))
invisible(sir$add_parameter(list(beta = 2e-4, gamma = 0.01)))
invisible(sir$add_reaction("beta * S * I", list("S", "I"), list("I", "I"), name = "infection"))
invisible(sir$add_reaction("gamma * I", list("I"), list("R"), name = "recovery"))
invisible(sir$set_stop("I == 0"))  # a replicate stops when the epidemic is over

# Vectors go in as list(...): a length-1 R vector would become a scalar.
result <- rules$ensemble(
  sir, end = 200, report_times = list(0, 50, 100, 200),
  replicates = 200L, seed = 1L, threads = 2L,
  method = "static_direct", keep_paths = "all"
)$simulate()

# pandas results arrive as data.frames (rows = report times or replicates).
mean <- result$statistics$mean
records <- result$records

# Kept paths: an xarray DataArray -> R array with named dimensions,
# [replicate, time, variable].
paths <- result$paths
dims <- unlist(paths$dims)
arr <- paths$values
coord <- function(name) as.vector(paths$coords$get(name)$values)
dimnames(arr) <- setNames(lapply(dims, coord), dims)
long <- as.data.frame.table(arr, responseName = "count")  # tidy data.frame

print(mean)
print(table(records$termination))
print(arr[1:3, , "I"])
print(head(long, 3))
```

```text
          I       R       S
0     1.000   0.000 999.000
50  675.960  99.960 224.080
100 538.365 431.440  30.195
200 197.345 772.695  29.960

      completed         stopped     event_limit  sampling_limit  instance_limit
            194               6               0               0               0
wall_time_limit       cancelled          failed    output_limit
              0               0               0               0
         time
replicate 0  50 100 200
        0 1 705 583 218
        1 1 825 511 176
        2 1 802 500 185
  replicate time variable count
1         0    0        I     1
2         1    0        I     1
3         2    0        I     1
```

Six of the 200 epidemics die out within the first 24 days, five of them
before any infection. There the stop condition `I == 0` holds, and these
replicates are recorded as `stopped`.

## Event logs

`event_columns()` returns a DataArray with the dimensions
`[replicate, entry, variable]`, so `result$event_columns()$values` is an R
array in that order. Its coordinates `time`, `rule` and `rule_name` are
matrices over replicate and entry:

```r
log <- result$event_columns()  # with event_log = list(max_entries = 5L)
values <- log$values                         # [replicate, entry, variable]
times <- log$coords$get("time")$values       # [replicate, entry]
names <- log$coords$get("rule_name")$values  # [replicate, entry]
one <- result$event_frame(0L)                # one replicate's log as a data.frame
```

## Pitfalls

- **Counts must be R integers.** Write `999L`, `replicates = 200L`,
  `seed = 1L` and `threads = 2L`. A plain `200` is a double, and it is
  rejected with "must be an integer" or "'float' object cannot be
  interpreted as an integer".
- **Pass vectors as `list(...)`.** reticulate turns a length-1 R vector into
  a single number, so `report_times = c(10)` fails with "'float' object is
  not iterable". Write `report_times = list(10)` instead. The same applies
  to reactant and product lists: write `list("I")`, not `"I"`.
- **Read multi-indexed frames as arrays.** reticulate drops the two-level
  index (replicate, time) of `paths_frame()`, and the (replicate, entry)
  index of `event_frame()`. Use `result$paths` and
  `result$event_columns()`, as above; `event_frame(id)` for one replicate is
  fine.
- **Read coordinates with `coords$get(...)`.** `paths$variable` returns the
  whole data array instead, because the name clashes with an attribute
  that xarray defines itself.
- **A model is not a list.** Set a stop with `model$set_stop("...")`; the
  [user manual](mtdna-user-manual.md#12-building-models-in-python) lists
  the other builder calls.
- **Large seeds.** R integers are 32-bit, so seeds above 2,147,483,647 cannot
  be written as `L` literals.
- **Version differences.** In release 0.5.0, `paths` was ordered
  `[time, variable, replicate]`. The example reads the dimension names from
  the DataArray, so it works with both orders; only `arr[1:3, , "I"]`
  assumes 0.6.0.

## Without Python

The `crn-craft` command-line tool (`cargo install crn-craft-cli`, then
`crn-craft run MODEL REQUEST`) reads and writes JSON, which R can handle
with jsonlite. You then write the model in the
[model language](rules-language.md) and the run request in JSON yourself;
the request must name the model file and its SHA-256 hash
([user manual, section 15](mtdna-user-manual.md#15-the-command-line-and-json-requests)).
