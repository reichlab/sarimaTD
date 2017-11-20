library(sarimautils)

context("initial transform")

################################################################################
## Tests for initial transformations
################################################################################

test_that("do_initial_transform/invert_initial_transform: transformation box-cox", {
	## data all positive
	y <- 1:10

	bc_params <- car::powerTransform(y, family = "bcnPower")
  y_trans <- do_initial_transform(
		y = y,
		transformation = "box-cox",
		bc_params = bc_params)

  y_trans2 <- car::bcnPower(
    U = y,
	  lambda = bc_params$lambda,
	  gamma = bc_params$gamma)

	y_out <- invert_initial_transform(
		y = y_trans,
		transformation = "box-cox",
		bc_params = bc_params)

	expect_equal(y_trans, y_trans2)
	expect_equal(y, y_out)



	## data all non-negative
	y <- 0:10

	bc_params <- car::powerTransform(y, family = "bcnPower")
  y_trans <- do_initial_transform(
		y = y,
		transformation = "box-cox",
		bc_params = bc_params)

  y_trans2 <- car::bcnPower(
    U = y,
	  lambda = bc_params$lambda,
	  gamma = bc_params$gamma)

	y_out <- invert_initial_transform(
		y = y_trans,
		transformation = "box-cox",
		bc_params = bc_params)

	expect_equal(y_trans, y_trans2)
	expect_equal(y, y_out)

	## data negative and positive
	y <- (-10):10

	bc_params <- car::powerTransform(y, family = "bcnPower")
  y_trans <- do_initial_transform(
		y = y,
		transformation = "box-cox",
		bc_params = bc_params)

  y_trans2 <- car::bcnPower(
    U = y,
	  lambda = bc_params$lambda,
	  gamma = bc_params$gamma)

	y_out <- invert_initial_transform(
		y = y_trans,
		transformation = "box-cox",
		bc_params = bc_params)

	expect_equal(y_trans, y_trans2)
	expect_equal(y, y_out)
})

test_that("do_initial_transform/invert_initial_transform: transformation log", {
	## data all positive
	y <- 1:10

  y_trans <- do_initial_transform(
		y = y,
		transformation = "log")

  y_trans2 <- log(y)

	y_out <- invert_initial_transform(
		y = y_trans,
		transformation = "log")

	expect_equal(y_trans, y_trans2)
	expect_equal(y, y_out)
})

test_that("do_initial_transform/invert_initial_transform: transformation none", {
	## data all positive
	y <- 1:10

  y_trans <- do_initial_transform(
		y = y,
		transformation = "none")

  y_trans2 <- y

	y_out <- invert_initial_transform(
		y = y_trans,
		transformation = "none")

	expect_equal(y_trans, y_trans2)
	expect_equal(y, y_out)
})
