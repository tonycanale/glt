# stepstester_delta.R
# Manual checker for OLDSTUFF/update.delta.R.
#
# Run this from the package root:
#   source("stepstester_delta.R")
#   stepstester_delta()

source("OLDSTUFF/aux.R")
source("OLDSTUFF/update.delta.R")
res <- update.delta(Delta,Eta,tau,y,hyperpar)
check_proposal_plot(Delta, res)

stepstester_delta <- function(
  Delta_in = Delta,
  Eta_star_in = Eta,
  tau_in = tau,
  y_in = y,
  hyperpar_in = hyperpar,
  seed = 42,
  nrep = 100,
  plot = TRUE
) {
  Delta_start <- Delta_in
  storage.mode(Delta_start) <- "integer"

  set.seed(seed)
  Delta_updated <- update.delta(
    Delta = Delta_start,
    Eta_star = Eta_star_in,
    tau = tau_in,
    y = y_in,
    hyperpar = hyperpar_in
  )

  storage.mode(Delta_updated) <- "integer"

  if (isTRUE(plot)) {
    check_proposal_plot(Delta_start, Delta_updated)
  }

  seeds <- seed + seq_len(nrep) - 1L

  summary_mat <- vapply(
    seeds,
    function(s) {
      set.seed(s)
      out <- update.delta(
        Delta = Delta_start,
        Eta_star = Eta_star_in,
        tau = tau_in,
        y = y_in,
        hyperpar = hyperpar_in
      )

      c(
        changed = as.integer(!identical(out, Delta_start)),
        changed_entries = sum(out != Delta_start),
        same_pivots = as.integer(identical(get_pivots(out), get_pivots(Delta_start)))
      )
    },
    numeric(3)
  )

  list(
    seed = seed,
    single_run = list(
      Delta_start = Delta_start,
      Delta_updated = Delta_updated,
      changed = !identical(Delta_updated, Delta_start),
      changed_entries = sum(Delta_updated != Delta_start)
    ),
    repeated_runs = list(
      nrep = nrep,
      change_rate = mean(summary_mat["changed", ]),
      mean_changed_entries = mean(summary_mat["changed_entries", ]),
      pivot_preservation_rate = mean(summary_mat["same_pivots", ])
    )
  )
}
