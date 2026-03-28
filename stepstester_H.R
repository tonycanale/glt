# stepstester.R
# Wrapper matching the old update.H signature, delegating to update_H_cpp.
# Use this file to compare the C++ implementation against OLDSTUFF/update.H.R.
#
# Run this AFTER loading the package in the console:
#   devtools::load_all()
#   source("stepstester.R")

# ── Reproducible setup (mirrors OLDSTUFF/tests.R) ────────────────────────────
source("OLDSTUFF/aux.R")       # provides marg.lik / marg.lik.i
source("OLDSTUFF/update.H.R")  # provides the reference update.H()

# Grab update_H_cpp from the loaded (not necessarily installed) namespace
update_H_cpp <- getFromNamespace("update_H_cpp", "gltfactor")

hyperpar <- list(alpha = 2, beta = 1, kappa = 1, sigma = 1, a_s = 1, b_s = 1,
                 a_nu = 2, b_nu = 2)
set.seed(1)
Tt <- 20; H <- 5; m <- 10
Eta    <- matrix(rnorm(Tt * H), H, Tt)
y      <- matrix(rnorm(Tt * m), Tt, m)
Lambda <- matrix(rnorm(m * H) * rbinom(m * H, 1, prob = 0.7), m, H)
Lambda[1:6, 4] <- 0
Lambda[1:8, 5] <- 0
Lambda[10, 5]  <- rnorm(1)
Delta  <- matrix(as.integer(Lambda != 0), m, H)
y      <- t(Lambda %*% Eta) + matrix(rnorm(Tt * m, sd = 0.2), Tt, m)
nu     <- 0.5

# ── Wrapper: same signature as the old update.H() ────────────────────────────
update.H.new <- function(y, Delta, Eta, hyperpar, nu, q = 0.5) {
  Delta_int <- Delta
  storage.mode(Delta_int) <- "integer"
  update_H_cpp(
    y        = y,
    Delta_in = Delta_int,
    Eta_in   = Eta,
    hyperpar = hyperpar,
    nu       = nu,
    q        = q
  )
}

# ── Single-step comparison ────────────────────────────────────────────────────
set.seed(42)
ref <- update.H(y, Delta, Eta, hyperpar, nu, q = 0.5)

set.seed(42)
new <- update.H.new(y, Delta, Eta, hyperpar, nu, q = 0.5)

cat("=== Single-step comparison ===\n")
cat("ref$increased :", ref$increased,  " | new$increased :", new$increased,  "\n")
cat("ref$accepted  :", ref$accepted,   " | new$accepted  :", new$accepted,   "\n")
cat("ref H         :", ncol(ref$Delta), " | new H         :", ncol(new$Delta), "\n")
cat("Delta match   :", isTRUE(all.equal(ref$Delta,  new$Delta)),  "\n")
cat("Eta match     :", isTRUE(all.equal(ref$Eta,    new$Eta)),    "\n")

# ── Monte Carlo acceptance-rate comparison ────────────────────────────────────
N <- 10000
results <- replicate(N, {
  q_val <- if (H > 1 && H < m) 0.5 else 1.0
  r  <- update.H(y, Delta, Eta, hyperpar, nu, q = q_val)
  Di <- Delta; storage.mode(Di) <- "integer"
  n  <- update_H_cpp(y, Di, Eta, hyperpar, nu, q_val)
  c(ref_accept   = as.integer(r$accepted),
    new_accept   = as.integer(n$accepted),
    ref_increase = as.integer(r$increased),
    new_increase = as.integer(n$increased))
}, simplify = TRUE)

cat("\n=== Monte Carlo rates over", N, "draws ===\n")
cat("Acceptance rate  — ref:", mean(results["ref_accept",  ]),
                    " new:", mean(results["new_accept",  ]), "\n")
cat("Birth-move rate  — ref:", mean(results["ref_increase",]),
                    " new:", mean(results["new_increase",]), "\n")
