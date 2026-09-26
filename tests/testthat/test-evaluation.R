make_eval_interactions <- function() {
  data.frame(
    parasite = c("pA 1", "pA 1", "pA 2", "pB 1", "pB 1", "pB 2"),
    host = c("host1", "host2", "host1", "host2", "host3", "host3"),
    stringsAsFactors = FALSE
  )
}

test_that("evaluate_phylopred_model returns sensible metrics", {
  pairs <- prepare_pair_data(make_eval_interactions(), make_simple_phydist())
  ev <- evaluate_phylopred_model(pairs)
  expect_s3_class(ev, "phylopred_evaluation")
  expect_true(ev$mcfadden_r2 >= 0 && ev$mcfadden_r2 <= 1)
  expect_true(ev$tjur_r2 >= -1 && ev$tjur_r2 <= 1)
  expect_true(ev$auc >= 0 && ev$auc <= 1)
  expect_true(is.na(ev$cv_auc))
  expect_equal(ev$n_parasites, 4)
  expect_output(print(ev), "McFadden")
})

test_that("evaluate_phylopred_model grouped CV works and is reproducible", {
  pairs <- prepare_pair_data(make_eval_interactions(), make_simple_phydist())
  ev1 <- evaluate_phylopred_model(pairs, cv = TRUE, k = 2, seed = 1)
  ev2 <- evaluate_phylopred_model(pairs, cv = TRUE, k = 2, seed = 1)
  expect_equal(ev1$cv_auc, ev2$cv_auc)
  expect_false(is.na(ev1$cv_auc))
  expect_error(evaluate_phylopred_model(pairs, cv = TRUE, k = 10), "at least 'k'")
  expect_error(evaluate_phylopred_model(pairs, k = 1), ">= 2")
})

test_that(".auc matches a known value", {
  y <- c(0, 0, 1, 1)
  score <- c(0.1, 0.4, 0.35, 0.8)
  expect_equal(phylopred:::.auc(y, score), 0.75)
})

test_that("phylopred_confidence_ellipse works from glm and bootstrap", {
  pairs <- prepare_pair_data(make_eval_interactions(), make_simple_phydist())
  fit <- stats::glm(suscept ~ phydist, data = pairs, family = stats::binomial())
  ell <- phylopred_confidence_ellipse(fit, n_points = 50)
  expect_equal(nrow(ell), 50)
  expect_equal(colnames(ell), c("(Intercept)", "phydist"))
  expect_equal(attr(ell, "center"), stats::coef(fit)[c("(Intercept)", "phydist")])

  boot <- cluster_bootstrap_phylopred(pairs, n_boot = 20, seed = 3)
  ell_b <- phylopred_confidence_ellipse(boot, n_points = 30)
  expect_equal(nrow(ell_b), 30)

  expect_error(phylopred_confidence_ellipse(fit, terms = c("a", "b")), "not found")
  expect_error(phylopred_confidence_ellipse(list()), "must be a glm")
  expect_error(phylopred_confidence_ellipse(fit, level = 2), "between 0 and 1")
})
