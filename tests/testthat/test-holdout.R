test_that("validate_host_prediction validates inputs", {
  expect_error(validate_host_prediction(list(), phy_dist, "a"), "data.frame")
  expect_error(
    validate_host_prediction(beetleTreeInteractions, phy_dist, "a",
                             prop_holdout = 1),
    "0, 1"
  )
  expect_error(
    validate_host_prediction(beetleTreeInteractions, phy_dist,
                             "not_a_parasite"),
    "at least 3 known hosts"
  )
})

test_that("glm holdout validation returns AUCs and is reproducible", {
  parasite <- names(sort(table(beetleTreeInteractions[[1]]),
                         decreasing = TRUE))[1]
  v1 <- validate_host_prediction(beetleTreeInteractions, phy_dist, parasite,
                                 method = "glm", n_rep = 3, seed = 1)
  v2 <- validate_host_prediction(beetleTreeInteractions, phy_dist, parasite,
                                 method = "glm", n_rep = 3, seed = 1)
  expect_s3_class(v1, "phylopred_holdout")
  expect_length(v1$auc, 3)
  expect_true(all(v1$auc >= 0 & v1$auc <= 1))
  expect_equal(v1$auc, v2$auc)
  expect_output(print(v1), "Holdout AUC")
})

test_that("pu holdout validation works", {
  parasite <- names(sort(table(beetleTreeInteractions[[1]]),
                         decreasing = TRUE))[1]
  v <- validate_host_prediction(beetleTreeInteractions, phy_dist, parasite,
                                method = "pu", n_rep = 2, seed = 1,
                                max_iter = 5)
  expect_length(v$auc, 2)
  expect_true(all(v$auc >= 0 & v$auc <= 1))
})
