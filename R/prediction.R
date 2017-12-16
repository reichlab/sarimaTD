## functions for SARIMA prediction

#' Simulate predictive trajectories from an ARIMA model
#'
#' This is directly taken from the forecast.Arima function from the forecast
#' package, but it returns the simulated trajectories which were not returned in
#' the original function definition.  This function calls
#' forecast:::simulate.Arima, which is NOT exported from the forecast package!!
#'
#' @param object an Arima fit object (with class "Arima")
#' @param h number of time steps forwards to simulate
#' @param npaths number of sample trajectories to simulate
#'
#' @return an npaths by h matrix with simulated values
#'
#' @export
sample_predictive_trajectories_arima <- function(
  object,
  h = ifelse(object$arma[5] > 1, 2 * object$arma[5], 10),
  bootstrap = FALSE,
  npaths = 1000,
  ...) {
    sim <- matrix(NA, nrow = npaths, ncol = h)

    for (i in 1:npaths) {
        try(sim[i, ] <- forecast:::simulate.Arima(object, nsim = h), silent = TRUE)
    }

    return(sim)
}

#' A wrapper around sample_predictive_trajectories_arima
#'
#' This function does a few things worth noting.  It linearly interpolates
#' missing values on the interior of the time series, which is necessary for
#' some model fits with moving average components.  It also handles
#' transformations and possible seasonal differencing that may have been done
#' in the fit_sarima function.
#'
#' @param sarima_fit a sarima fit, probably as returned by fit_sarima
#' @param h number of time steps forwards to simulate
#' @param npaths number of sample trajectories to simulate
#'
#' @return an npaths by h matrix with simulated values
#'
#' @export
sample_predictive_trajectories_arima_wrapper <- function(
  sarima_fit,
  newdata,
  h = 1,
  npaths = 1
) {
  if(is.null(sarima_fit$sarima_utils_call)) {
    transformation <- "none"
    seasonal_difference <- FALSE
  } else {
    transformation <- sarima_fit$sarima_utils_call$transformation
    seasonal_difference <-
      sarima_fit$sarima_utils_call$seasonal_difference
    ts_frequency <- sarima_fit$arma[5]
  }

  ## Update SARIMA fit object with transformed and seasonally differenced data
  ## Initial transformation, if necessary
  if(identical(transformation, "box-cox")) {
    est_bc_params <- sarima_fit$sarima_utils_est_bc_params
  } else {
    est_bc_params <- NULL
  }
  transformed_y <- do_initial_transform(
    y = newdata,
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

  ## Drop leading missing values, fill in internal missing values via linear
  ## interpolation.  This is necessary to ensure non-missing predictions if the
  ## sarima model has a moving average component.
  ## deal with internal missing or infinite values in y that can
  ## result in simulated trajectories of all NAs if the model has a moving
  ## average component.  Here we do this by linear interpolation.
  ##
  ## Another (better?) solution would be to write a version of stats::filter
  ## that does the "right thing" with NAs if filter coefficients are 0 and then
  ## use that function in forecast:::myarima.sim
  interpolated_y <- interpolate_and_clean_missing(differenced_y)

  ## Update fit with transformed, differenced, and interpolated data
  updated_sarima_fit <- forecast::Arima(
      interpolated_y,
      model = sarima_fit)

  ## Sample trajectories on the transformed and differenced scale
  raw_trajectory_samples <- sample_predictive_trajectories_arima(
      updated_sarima_fit,
      h = h,
      npaths = npaths)

  ## Sampled trajectories are of seasonally differenced transformed time series
  ## Get to trajectories for originally observed time series ("orig") by
  ## adding seasonal lag of incidence and inverting the transformation
  orig_trajectory_samples <- raw_trajectory_samples
  if(seasonal_difference) {
    for(i in seq_len(npaths)) {
      orig_trajectory_samples[i, ] <-
        invert_seasonal_difference(
          dy = raw_trajectory_samples[i, ],
          y = transformed_y,
          ts_frequency = ts_frequency)
      orig_trajectory_samples[i, ] <-
        invert_initial_transform(
          y = orig_trajectory_samples[i, ],
          transformation = transformation,
          bc_params = est_bc_params)
    }
  } else {
    for(i in seq_len(npaths)) {
      orig_trajectory_samples[i, ] <-
        invert_initial_transform(
          y = raw_trajectory_samples[i, ],
          transformation = transformation,
          bc_params = est_bc_params)
    }
  }

  return(orig_trajectory_samples)
}
