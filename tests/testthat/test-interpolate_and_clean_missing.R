library(sarimautils)

context("clean missing and infinite values")

################################################################################
## Tests for interpolating and cleaning missing and infinite values
################################################################################

test_that("interpolate_and_clean_missing NAs", {
	## data all positive
	y <- 1:104
	y[3:30] <- NA
	y <- c(rep(NA, 10), y)

  y_out <- interpolate_and_clean_missing(y)

	expect_equal(y_out, 1:104)
})

test_that("interpolate_and_clean_missing Inf", {
	## data all positive
	y <- 1:104
	y[3:30] <- Inf
	y <- c(rep(Inf, 10), y)

  y_out <- interpolate_and_clean_missing(y)

	expect_equal(y_out, 1:104)
})

test_that("interpolate_and_clean_missing -Inf", {
	## data all positive
	y <- 1:104
	y[3:30] <- -Inf
	y <- c(rep(-Inf, 10), y)

  y_out <- interpolate_and_clean_missing(y)

	expect_equal(y_out, 1:104)
})

test_that("interpolate_and_clean_missing -Inf NA Inf", {
	## data all positive
	y <- 1:104
	y[3:10] <- -Inf
	y[11:20] <- NA
	y[21:10] <- Inf
	y <- c(rep(-Inf, 10), rep(NA, 10), rep(Inf, 10), y)

  y_out <- interpolate_and_clean_missing(y)

	expect_equal(y_out, 1:104)
})
