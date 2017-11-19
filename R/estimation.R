## functions for SARIMA estimation

#' Estimate SARIMA model
#'
#' @param y a univariate time series or numeric vector.
#' @param ts_frequency frequency of time series.  Must be provided if y is not
#'   of class "ts".  See the help for stats::ts for more.
#' @param transformation character specifying transformation type:
#'   "box-cox", "log", "forecast-box-cox", or "none".  See details for more.
#' @param seasonal_difference boolean; take a seasonal difference before passing
#'   to auto.arima?
#' @param auto.arima_d order of first differencing argument to auto.arima.
#' @param auto.arima_D order of seasonal differencing argument to auto.arima.
#'
#' @return a SARIMA model fit
#'
#' @details This function is a wrapper around forecast::auto.arima, providing
#'   some useful defaults for preliminary transformations of the data.
#'   Formal and informal experimentation has shown these preliminary
#'   transformations to be helpful with a few infectious disease time series
#'   data sets.  Note that if any transformation was specified or the
#'   seasonal_difference argument was TRUE in the call to this function, only
#'   prediction/forecast utilities provided by the sarima_utils package can be
#'   used!  We have found that using the default arguments for transformation,
#'   seasonal_difference, d, and D, yields good performance.
#'
#' @export
fit_sarima <- function(
  y,
  ts_frequency,
  transformation = "box-cox",
  seasonal_difference = TRUE,
  d = NA,
  D = NA) {
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
    est_bc_params <- car::powerTransform(y, family = "bcnPower")
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
      lambda = lambda)
  } else {
    ## box-cox or log transformation (already done) or no transformation
    sarima_fit <- forecast::auto.arima(differenced_y,
      d = d,
      D = D,
      stationary = TRUE)
  }

  sarima_fit$sarima_utils_call <- match.call()
  if(identical(transformation, "box-cox")) {
    sarima_fit$sarima_utils_est_bc_params <- est_bc_params
  }

  return(sarima_fit)
}
