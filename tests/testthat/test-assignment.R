test_that("assign_subtypes labels a clear single-signature winner correctly", {
  scores <- data.frame(
    Sig1_UCell = c(0.9, 0.1, 0.1, 0.1, 0.1),
    Sig2_UCell = c(0.1, 0.1, 0.1, 0.1, 0.1)
  )
  rownames(scores) <- paste0("cell", 1:5)

  res <- assign_subtypes(scores, quantile_prob = 0.75)

  expect_equal(res$Subtype[1], "Sig1")
  expect_true(all(res$Subtype[2:5] %in% c("Unknown", "Sig1", "Sig2")))
  expect_false(any(grepl("_UCell$", res$Subtype)))
})

test_that("assign_subtypes labels cells passing multiple signatures as Mixed", {
  scores <- data.frame(
    Sig1_UCell = c(0.9, 0.1),
    Sig2_UCell = c(0.9, 0.1)
  )
  rownames(scores) <- c("cellA", "cellB")

  res <- assign_subtypes(scores, quantile_prob = 0.5)

  expect_true("Mixed" %in% res$Subtype)
})

test_that("assign_subtypes labels cells passing nothing as Unknown", {
  scores <- data.frame(
    Sig1_UCell = rep(0.5, 10),
    Sig2_UCell = rep(0.5, 10)
  )
  scores$Sig1_UCell[1] <- 0.99  # one outlier so quantile has a real threshold
  rownames(scores) <- paste0("cell", 1:10)

  res <- assign_subtypes(scores, quantile_prob = 0.99)

  expect_true(all(res$Subtype[-1] == "Unknown"))
})

test_that("assign_subtypes rejects invalid quantile_prob", {
  scores <- data.frame(Sig1_UCell = c(0.1, 0.2))
  expect_error(assign_subtypes(scores, quantile_prob = 0), "between 0 and 1")
  expect_error(assign_subtypes(scores, quantile_prob = 1.5), "between 0 and 1")
})

test_that("assign_subtypes returns expected columns", {
  scores <- data.frame(Sig1_UCell = c(0.1, 0.9), Sig2_UCell = c(0.5, 0.5))
  res <- assign_subtypes(scores)
  expect_true(all(c("Cell", "n_signatures_passed", "Subtype") %in% colnames(res)))
  expect_equal(nrow(res), 2)
})
