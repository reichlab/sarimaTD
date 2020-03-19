library(sarimaTD)

context("seasonal difference")

################################################################################
## Tests for seasonal differencing
################################################################################

test_that("do_difference/invert_difference: d = 1, D = 0", {
  ## data all positive
  y <- 1:104
  d <- 1
  D <- 0
  
  y_diff <- do_difference(y = y, d = d, D = D, frequency = 52)
  y_orig <- invert_difference(dy = y_diff[(d+D+1):104], y = y[d], d = d, D = D, frequency = 52)
  
  expect_equal(y_diff, ts(c(rep(NA, d+D), rep(1, 104-(d+D))), freq = 52))
  expect_equal(ts(y[(d+D+1):104], freq = 52), y_orig)
})

test_that("do_difference/invert_difference: d = 2, D = 0", {
  ## data all positive
  y <- 1:104
  d <- 2
  D <- 0
  
  y_diff <- do_difference(y = y, d = d, D = D, frequency = 52)
  y_orig <- invert_difference(dy = y_diff[(d+D+1):104], y = y[1:d], d = d, D = D, frequency = 52)
  
  expect_equal(y_diff, ts(c(rep(NA, d+D), rep(0, 104-(d+D))), freq = 52))
  expect_equal(ts(y[(d+D+1):104], freq = 52), y_orig)
})


test_that("do_difference/invert_difference: d = 0, D = 1", {
	## data all positive
	y <- 1:104

  y_diff <- do_difference(y = y, d = 0, D = 1, frequency = 52)
	y_orig <- invert_difference(dy = y_diff[53:104], y = y[1:52], d = 0, D = 1, frequency = 52)

	expect_equal(y_diff, ts(c(rep(NA, 52), rep(52, 52)), freq = 52))
	expect_equal(ts(y[53:104], freq = 52), y_orig)
})



test_that("do_difference/invert_difference: d = 0, D = 2", {
  ## data all positive
  y <- 1:156
  
  y_diff <- do_difference(y = y, d = 0, D = 2, frequency = 52)
  y_orig <- invert_difference(dy = y_diff[105:156], y = y[1:104], d = 0, D = 2, frequency = 52)
  
  expect_equal(y_diff, ts(c(rep(NA, 104), rep(0, 52)), freq = 52))
  expect_equal(ts(y[105:156], freq = 52), y_orig)
})



test_that("do_difference/invert_difference: d = 2, D = 2", {
  ## data all positive
  y <- rlnorm(n = 208)
  d <- 2
  D <- 1
  
  y_diff <- do_difference(y = y, d = d, D = D, frequency = 52)
  y_orig <- invert_difference(dy = y_diff[(d+52*D+1):length(y)], y = y[seq_len(d + 52*D)], d = d, D = D, frequency = 52)
  
  expect_equal(ts(y[(d+52*D+1):length(y)], freq = 52), y_orig)
})
