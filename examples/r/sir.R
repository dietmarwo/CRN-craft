# CRN-craft from R through reticulate: an SIR ensemble with kept paths.
#
# One-time setup:
#   install.packages("reticulate")
#   reticulate::virtualenv_create("crn-craft", packages = "crn-craft")

library(reticulate)
use_virtualenv(Sys.getenv("CRN_CRAFT_VENV", "crn-craft"), required = TRUE)
crn <- import("crn_craft")

# Counts must be R integers (999L); rates and parameters are doubles.
sir <- crn$Model("sir")
invisible(sir$add_species(list(S = 999L, I = 1L, R = 0L)))
invisible(sir$add_parameter(list(beta = 2e-4, gamma = 0.01)))
invisible(sir$add_reaction("beta * S * I", list("S", "I"), list("I", "I"), name = "infection"))
invisible(sir$add_reaction("gamma * I", list("I"), list("R"), name = "recovery"))
invisible(sir$set_stop("I == 0"))  # a replicate stops when the epidemic is over

# Vectors go in as list(...): a length-1 R vector would become a scalar.
result <- crn$ensemble(
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
