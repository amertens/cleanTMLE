test_that(".hash ignores formula environments and sees content", {
  f1 <- local(~ a + b)
  f2 <- local({ x <- 1; ~ a + b })
  expect_identical(.hash(list(r = f1)), .hash(list(r = f2)))
  expect_false(identical(.hash(data.frame(a = 1:3)),
                         .hash(data.frame(a = c(1L, 2L, 4L)))))
  expect_identical(.hash(function(x) x + 1), .hash(function(x) x + 1))
})

test_that(".make_folds is balanced and seeded, and K = 1 gives one fold", {
  f <- .make_folds(10, 3, seed = 4)
  expect_setequal(unique(f), 1:3)
  expect_true(max(table(f)) - min(table(f)) <= 1)
  expect_identical(f, .make_folds(10, 3, seed = 4))
  expect_identical(.make_folds(5, 1, seed = 4), rep(1L, 5))
})

test_that(".bound keeps probabilities inside (eps, 1 - eps)", {
  expect_equal(.bound(c(0, 0.5, 1), 0.01), c(0.01, 0.5, 0.99))
})

test_that(".hash sees a function's body and its captured variables", {
  expect_false(identical(.hash(function(x) x + 1), .hash(function(x) x + 2)))
  mk <- function(delta) function(a, X) stats::plogis(delta * a)
  expect_identical(.hash(mk(1)), .hash(mk(1)))
  expect_false(identical(.hash(mk(1)), .hash(mk(2))))
  # a captured function contributes its text only, so a self-reference cannot recurse
  mk2 <- function(g) function(a) g(a)
  expect_identical(.hash(mk2(function(z) z)), .hash(mk2(function(z) z)))
  expect_false(identical(.hash(mk2(function(z) z)), .hash(mk2(function(z) z + 1))))
})
