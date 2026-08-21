test_that("read_gmt parses a well-formed gmt file", {
  tmp <- tempfile(fileext = ".gmt")
  writeLines(c(
    "SigA\tdesc A\tGENE1\tGENE2\tGENE3",
    "SigB\tdesc B\tGENE4"
  ), tmp)

  sigs <- read_gmt(tmp)

  expect_named(sigs, c("SigA", "SigB"))
  expect_equal(sigs$SigA, c("GENE1", "GENE2", "GENE3"))
  expect_equal(sigs$SigB, "GENE4")
})

test_that("read_gmt errors on missing file", {
  expect_error(read_gmt("does_not_exist.gmt"), "not found")
})

test_that("read_gmt errors on malformed lines", {
  tmp <- tempfile(fileext = ".gmt")
  writeLines(c("SigA\tonly_description_no_genes"), tmp)
  expect_error(read_gmt(tmp), "Malformed")
})

test_that("read_gmt errors on duplicate signature names", {
  tmp <- tempfile(fileext = ".gmt")
  writeLines(c(
    "SigA\tdesc\tGENE1",
    "SigA\tdesc\tGENE2"
  ), tmp)
  expect_error(read_gmt(tmp), "Duplicate")
})

test_that("built-in HGSOC panel loads and has expected signatures", {
  sigs <- load_builtin_signatures("HGSOC")
  expect_setequal(names(sigs), c("Differentiated", "Immunoreactive", "Mesenchymal", "Proliferative"))
  expect_true("MUC16" %in% sigs$Differentiated)
  expect_true("MKI67" %in% sigs$Proliferative)
})

test_that("load_builtin_signatures errors informatively for unknown panels", {
  expect_error(load_builtin_signatures("LUNG_CANCER"), "manual")
})
