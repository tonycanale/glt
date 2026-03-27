# stepstester_pivots.R
# Wrapper and checks for update.pivots_cpp against OLDSTUFF/update.pivots.R.
#
# Run this AFTER loading the package in the console:
#   devtools::load_all()
#   source("stepstester_pivots.R")

source("OLDSTUFF/aux.R")
source("OLDSTUFF/update.pivots.R")

# Silence the plotting side effect in the old reference implementation.
check_proposal_plot <- function(Delta, Deltastar) {
  invisible(NULL)
}

update.pivots.new <- getFromNamespace("update.pivots.new", "gltfactor")

hyperpar <- list(alpha = 2, beta = 1, kappa = 1, sigma = 1, a_s = 1, b_s = 1,
                 a_nu = 2, b_nu = 2)
set.seed(1)
Tt <- 20; H <- 5; m <- 10
Eta <- matrix(rnorm(Tt * H), H, Tt)
y <- matrix(rnorm(Tt * m), Tt, m)

Bad_Delta <- rbind(
  matrix(0L, 3, 5),
  diag(1L, 5),
  matrix(1L, 2, 5)
)
Delta <- Bad_Delta
storage.mode(Delta) <- "integer"

# ── Single-step comparison ────────────────────────────────────────────────────
set.seed(42)
ref <- update.pivots(y, Delta, Eta, hyperpar)

set.seed(42)
new <- update.pivots.new(y, Delta, Eta, hyperpar)

cat("=== Single-step comparison ===\n")
cat("Delta match   :", isTRUE(all.equal(ref$Delta, new$Delta)), "\n")
cat("Eta match     :", isTRUE(all.equal(ref$Eta,   new$Eta)),   "\n")
cat("ell match     :", isTRUE(all.equal(apply(ref$Delta, 2, match, x = 1), new$ell)), "\n")
cat("ref changed   :", !isTRUE(all.equal(ref$Delta, Delta)), "\n")
cat("new changed   :", !isTRUE(all.equal(new$Delta, Delta)), "\n")

# ── Monte Carlo comparison under matched seeds ───────────────────────────────
N <- 5000
seeds <- sample.int(.Machine$integer.max, N)

results <- vapply(seeds, function(s) {
  set.seed(s)
  r <- update.pivots(y, Delta, Eta, hyperpar)

  set.seed(s)
  n <- update.pivots.new(y, Delta, Eta, hyperpar)

  c(
    ref_changed = as.integer(!isTRUE(all.equal(r$Delta, Delta))),
    new_changed = as.integer(!isTRUE(all.equal(n$Delta, Delta))),
    exact_match = as.integer(isTRUE(all.equal(r$Delta, n$Delta)))
  )
}, numeric(3))

cat("\n=== Monte Carlo rates over", N, "draws ===\n")
cat("Change rate      — ref:", mean(results["ref_changed", ]),
    " new:", mean(results["new_changed", ]), "\n")
cat("Exact Delta match:", mean(results["exact_match", ]), "\n")
