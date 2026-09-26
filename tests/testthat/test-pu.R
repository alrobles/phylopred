test_that("fit_phylopred_pu validates inputs", {
  expect_error(fit_phylopred_pu(list(), phy_dist, "a"), "data.frame")
  expect_error(
    fit_phylopred_pu(beetleTreeInteractions, unname(phy_dist), "a"),
    "square named matrix"
  )
  expect_error(
    fit_phylopred_pu(beetleTreeInteractions, phy_dist, c("a", "b")),
    "single parasite"
  )
  expect_error(
    fit_phylopred_pu(beetleTreeInteractions, phy_dist, "not_a_parasite"),
    "No known hosts"
  )
})

test_that("fit_phylopred_pu predicts hosts and feeds the metrics", {
  parasite <- names(sort(table(beetleTreeInteractions[[1]]),
                         decreasing = TRUE))[1]
  pu <- fit_phylopred_pu(beetleTreeInteractions, phy_dist,
                      parasite = parasite, max_iter = 10, seed = 1)
  expect_s3_class(pu, "phylopred_pu")
  expect_equal(nrow(pu$prediction), nrow(phy_dist))
  expect_true(all(pu$prediction$pred >= 0 & pu$prediction$pred <= 1))
  expect_true(sum(pu$prediction$is_host) > 0)
  expect_output(print(pu), "PU-learning")

  m <- host_threshold_metric(pu$prediction$pred, pu$prediction$is_host)
  expect_s3_class(m, "phylopred_threshold_metric")
  roc <- phylo_partial_roc(pu$prediction$pred, pu$prediction$is_host,
                           error_rate = 0.95, n_boot = 20, seed = 1)
  expect_s3_class(roc, "phylopred_partial_roc")
})
