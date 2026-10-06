test_that("no design-stage object holds the primary outcome", {
  fx <- fixture_pipeline()
  y <- fx$d$outcomes$y
  for (nm in c("lock", "design", "sim", "nc", "dossier"))
    expect_false(contains_outcome(fx[[nm]], "y", y), label = nm)
})

test_that("the simulation runs on design data with no outcome column at all", {
  fx <- fixture_pipeline()
  expect_false("y" %in% names(fx$lock$data))
  expect_s3_class(fx$sim, "cr_simulation")
})
