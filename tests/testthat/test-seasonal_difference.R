library(sarimautils)

context("seasonal difference")

################################################################################
## Tests for seasonal differencing
################################################################################

test_that("do_seasonal_difference/invert_seasonal_difference", {
	## data all positive
	y <- 1:104

  y_diff <- do_seasonal_difference(y = y, ts_frequency = 52)
	y_orig <- invert_seasonal_difference(dy = y_diff[53:104], y = y[1:52], ts_frequency = 52)

	expect_equal(y_diff, ts(c(rep(NA, 52), rep(52, 52)), freq = 52))
	expect_equal(ts(y[53:104], freq = 52), y_orig)
})
