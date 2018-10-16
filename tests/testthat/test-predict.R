library(sarimaTD)

context("prediction")

################################################################################
## Tests for prediction
################################################################################

test_that("predict.sarimaTD, variable used for specification of transform during estimation", {
  # this test is related to Issue https://github.com/reichlab/sarimaTD/issues/15,
  # where errors were thrown on prediction if a variable was used as an argument to fit_sarima
  
  set.seed(1)
	y <- rlnorm(104)
	ts_frequency <- 52L
	transformation <- "log"
	seasonal_difference <- TRUE
	
	sarima_fit <- fit_sarima(
	  y = y,
	  ts_frequency = ts_frequency,
	  transformation = transformation,
	  seasonal_difference = seasonal_difference)
	
	simulate_result <- tryCatch(
  	sims <- simulate(
  	  object = sarima_fit,
  	  nsim = 1,
  	  seed = NULL,
  	  newdata = y,
  	  h = 1
  	),
    error = function(e) {"Threw an error!"}
	)
  
  expect_false(identical(simulate_result, "Threw an error!"))
})
