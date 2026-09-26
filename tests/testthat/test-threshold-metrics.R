make_scores <- function() {
  # Known hosts score high, unlabeled score low: strong signal.
  list(
    pred = c(0.9, 0.8, 0.85, 0.2, 0.1, 0.3, 0.15, 0.05),
    is_host = c(1, 1, 1, 0, 0, 0, 0, 0)
  )
}

test_that("host_threshold_metric computes a, b, c, f_gamma correctly", {
  s <- make_scores()
  m <- host_threshold_metric(s$pred, s$is_host, gamma = 1)
  expect_s3_class(m, "phylopred_threshold_metric")
  expect_equal(nrow(m), length(unique(s$pred)))
  # At the lowest threshold everything is predicted: a = 1, b = 1.
  expect_equal(m$a[1], 1)
  expect_equal(m$b[1], 1)
  expect_equal(m$c, 1 - m$a)
  expect_equal(m$f_gamma, m$a / (m$b + m$c))
  # Optimal threshold retains all hosts with minimal breadth (t = 0.8).
  expect_equal(attr(m, "optimal_threshold"), 0.8)
  expect_equal(attr(m, "optimal_f"), 1 / (3 / 8))
  expect_output(print(m), "Optimal threshold")
})

test_that("gamma shifts the optimum toward compact predictions", {
  s <- make_scores()
  m_low <- host_threshold_metric(s$pred, s$is_host, gamma = 0.5)
  m_high <- host_threshold_metric(s$pred, s$is_host, gamma = 10)
  expect_true(attr(m_high, "optimal_threshold") >= attr(m_low, "optimal_threshold"))
  expect_error(host_threshold_metric(s$pred, s$is_host, gamma = 0), "positive")
})

test_that("threshold metric validates inputs", {
  s <- make_scores()
  expect_error(host_threshold_metric("a", s$is_host), "numeric")
  expect_error(host_threshold_metric(s$pred, s$is_host[-1]), "same length")
  expect_error(host_threshold_metric(s$pred, rep(0, 8)), "at least one known host")
  expect_error(host_threshold_metric(s$pred, rep(1, 8)), "non-host")
  expect_error(host_threshold_metric(c(NA, s$pred[-1]), s$is_host), "NA")
})

test_that("phylo_partial_roc detects signal vs random predictions", {
  s <- make_scores()
  roc <- phylo_partial_roc(s$pred, s$is_host, n_boot = 100, seed = 42)
  expect_s3_class(roc, "phylopred_partial_roc")
  expect_true(roc$auc_ratio > 1)
  expect_true(roc$p_value < 0.5)
  expect_output(print(roc), "Partial AUC ratio")

  set.seed(7)
  pred_rand <- runif(200)
  is_host_rand <- c(rep(1, 40), rep(0, 160))
  roc_rand <- phylo_partial_roc(pred_rand, is_host_rand,
                                error_rate = 0.5, n_boot = 100, seed = 42)
  expect_true(abs(roc_rand$auc_ratio - 1) < 0.35)
})

test_that("phylo_partial_roc is reproducible and validates inputs", {
  s <- make_scores()
  r1 <- phylo_partial_roc(s$pred, s$is_host, n_boot = 50, seed = 1)
  r2 <- phylo_partial_roc(s$pred, s$is_host, n_boot = 50, seed = 1)
  expect_equal(r1$boot_ratios, r2$boot_ratios)
  expect_error(phylo_partial_roc(s$pred, s$is_host, error_rate = 1), "0, 1")
  expect_error(phylo_partial_roc(s$pred, s$is_host, n_boot = 0), ">= 1")
  expect_error(phylo_partial_roc(s$pred, s$is_host, boot_prop = 0), "0, 1")
})
