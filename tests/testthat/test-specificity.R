test_that("taxonomic_step_matrix follows the Poulin & Mouillot 1-5 scale", {
  tax <- data.frame(
    species = c("A a", "A b", "B a", "C a", "D a"),
    genus  = c("A", "A", "B", "C", "D"),
    family = c("F1", "F1", "F1", "F2", "F3"),
    order  = c("O1", "O1", "O1", "O1", "O2"),
    stringsAsFactors = FALSE
  )
  D <- taxonomic_step_matrix(tax)
  expect_true(is.matrix(D))
  expect_equal(rownames(D), tax$species)
  expect_equal(unname(diag(D)), c(0, 0, 0, 0, 0))
  expect_equal(D["A a", "A b"], 1)   # same genus
  expect_equal(D["A a", "B a"], 2)   # same family, different genus
  expect_equal(D["A a", "C a"], 3)   # same order, different family
  expect_equal(D["A a", "D a"], 4)   # different order (no class column -> capped)
  expect_true(isSymmetric(unname(D)))
})

test_that("poulin_std matches the S_TD formula on a toy distance matrix", {
  D <- make_simple_phydist() # host1-host2=1, host1-host3=2, host2-host3=1
  m <- poulin_std(D, c("host1", "host2", "host3"))
  expect_s3_class(m, "phylopred_std")
  expect_equal(m$n_hosts, 3)
  expect_equal(m$std, mean(c(1, 2, 1)))
  expect_equal(m$var_std, sum((c(1, 2, 1) - mean(c(1, 2, 1)))^2) * 2 / (3 * 2))
})

test_that("poulin_std is undefined (NA) for a single host", {
  D <- make_simple_phydist()
  expect_warning(m <- poulin_std(D, "host1"), "single-host")
  expect_equal(m$n_hosts, 1)
  expect_true(is.na(m$std))
  expect_true(is.na(m$var_std))
})

test_that("poulin_std leaves var_std NA for exactly two hosts", {
  D <- make_simple_phydist()
  m <- poulin_std(D, c("host1", "host2"))
  expect_equal(m$n_hosts, 2)
  expect_equal(m$std, 1)
  expect_true(is.na(m$var_std))
})

test_that("poulin_std drops hosts absent from phydist with a warning", {
  D <- make_simple_phydist()
  expect_warning(m <- poulin_std(D, c("host1", "host2", "ghost")),
                  "not found in 'phydist'")
  expect_equal(m$n_hosts, 2)
})

test_that("host_specificity_index reproduces poulin_std per parasite", {
  inter <- make_interactions()
  D <- make_simple_phydist()
  out <- host_specificity_index(inter, D)
  expect_s3_class(out, "phylopred_specificity")
  expect_named(out, c("parasite", "n_hosts", "std", "var_std"))

  row_pa1 <- out[out$parasite == "pA 1", ]
  expect_equal(row_pa1$n_hosts, 2)
  expect_equal(row_pa1$std, 1) # host1-host2 distance

  row_pa2 <- out[out$parasite == "pA 2", ]
  expect_equal(row_pa2$n_hosts, 1)
  expect_true(is.na(row_pa2$std))
})

test_that("host_specificity_index supports a group column", {
  inter <- make_interactions()
  D <- make_simple_phydist()
  out <- host_specificity_index(inter, D, group = "genus")
  expect_true("group" %in% colnames(out))
  expect_setequal(unique(out$group), c("pA", "pB"))
})
