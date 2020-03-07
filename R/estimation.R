## functions for SARIMA estimation

#' Estimate SARIMA model
#'
#' @param y a univariate time series or numeric vector.
#' @param ts_frequency frequency of time series.  Must be provided if y is not
#'   of class "ts".  See the help for stats::ts for more.
#' @param transformation character specifying transformation type:
#'   "box-cox", "log", "forecast-box-cox", or "none".  See details for more.
#' @param bc_gamma numeric offset used in Box-Cox transformation; gamma is
#'   added to all observations before transforming.  Default value of 0.5
#'   allows us to use the Box-Cox transform (which requires positive inputs)
#'   in case of observations of 0, and also ensures that the de-transformed
#'   values will always be at least -0.5, so that they round up to non-negative
#'   values.
#' @param seasonal_difference boolean; take a seasonal difference before passing
#'   to auto.arima?
#' @param d order of first differencing argument to auto.arima.
#' @param D order of seasonal differencing argument to auto.arima.
#' @param ... arguments passed on to forecast::auto.arima
#'
#' @return a SARIMA model fit
#'
#' @details This function is a wrapper around forecast::auto.arima, providing
#'   some useful defaults for preliminary transformations of the data.
#'   Formal and informal experimentation has shown these preliminary
#'   transformations to be helpful with a few infectious disease time series
#'   data sets.  Note that if any transformation was specified or the
#'   seasonal_difference argument was TRUE in the call to this function, only
#'   prediction/forecast utilities provided by the sarimaTD package can be
#'   used!  We have found that using the default arguments for transformation,
#'   seasonal_difference, d, and D, yields good performance.
#'
#' @export
fit_sarima <- function(
  y,
  ts_frequency,
  transformation = "box-cox",
  bc_gamma = 0.5,
  seasonal_difference = TRUE,
  d = NA,
  D = NA,
  ...) {
  ## Validate arguments
  if(!(is.numeric(y) || is.ts(y))) {
    stop("The argument y must be a numeric vector or object of class 'ts'.")
  }

  if(!is.ts(y) && missing(ts_frequency)) {
    stop("If y is not an object of class 'ts', the ts_frequency argument must be supplied.")
  }

  if(is.ts(y)) {
    ts_frequency <- frequency(y)
  }

  ## Initial transformation, if necessary
  if(identical(transformation, "box-cox")) {
    est_bc_params <- car::powerTransform(y + bc_gamma, family = "bcPower")
    est_bc_params <- list(
      lambda = est_bc_params$lambda,
      gamma = bc_gamma)
  }
  transformed_y <- do_initial_transform(
    y = y,
    transformation = transformation,
    bc_params = est_bc_params)

  ## Initial seasonal differencing, if necessary
  if(seasonal_difference) {
    differenced_y <- do_seasonal_difference(
      y = transformed_y,
      ts_frequency = ts_frequency)
  } else {
    differenced_y <- ts(transformed_y, frequency = ts_frequency)
  }

  ## Get SARIMA fit
  if(identical(transformation, "forecast-box-cox")) {
    ## box-cox transformation done by auto.arima
    lambda <- forecast::BoxCox.lambda(differenced_y)

    sarima_fit <- forecast::auto.arima(differenced_y,
      d = d,
      D = D,
      stationary = TRUE,
      lambda = lambda,
      optim.method = "L-BFGS-B",
      ...)
  } else {
    ## box-cox or log transformation (already done) or no transformation
    sarima_fit <- forecast::auto.arima(differenced_y,
      d = d,
      D = D,
      stationary = TRUE,
      optim.method = "L-BFGS-B",
      ...)
  }

  sarima_fit$sarimaTD_call <- match.call()
  for(param_name in c("y", "ts_frequency", "transformation", "seasonal_difference", "d", "D")) {
    sarima_fit[[paste0("sarimaTD_used_", param_name)]] <- get(param_name)
  }
  if(identical(transformation, "box-cox")) {
    sarima_fit$sarimaTD_est_bc_params <- est_bc_params
  }
  
  class(sarima_fit) <- c("sarimaTD", class(sarima_fit))

  return(sarima_fit)
}
