## utility functions for SARIMA fits

#' Do an initial transformation of time series data.
#'
#' @param y a univariate time series or numeric vector.
#' @param transformation character specifying transformation type:
#'   "box-cox", "log", "forecast-box-cox", or "none".
#' @param bc_params (required if transformation is "box-cox").  Parameters for
#'   the Box-Cox transformation.  This should be the result of a call to
#'   car::powerTransform(y, family = "bcnPower")
#'
#' @return a transformed object of the same class as y
#'
#' @export
do_initial_transform <- function(y, transformation, bc_params) {
  ## Update SARIMA fit object with transformed and seasonally differenced data
  if(identical(transformation, "log")) {
    transformed_y <- log(y)
  } else if(identical(transformation, "box-cox")) {
    if(missing(bc_params)) {
      stop("bc_params must be provided if transformation is 'box-cox'")
    }
    transformed_y <- car::bcnPower(
      U = y,
      lambda = bc_params$lambda,
      gamma = bc_params$gamma)
  } else if(identical(transformation, "none") ||
      identical(transformation, "forecast-box-cox")) {
    transformed_y <- y
  } else {
    stop("Invalid transformation: must be one of 'box-cox', 'log', 'forecast-box-cox', or 'none'.")
  }

  return(transformed_y)
}

#' Invert an initial transformation of time series data.
#'
#' @param y a univariate time series or numeric vector.
#' @param transformation character specifying transformation type:
#'   "box-cox", "log", "forecast-box-cox", or "none".
#' @param bc_params (required if transformation is "box-cox").  Parameters for
#'   the Box-Cox transformation.  This should be the result of a call to
#'   car::powerTransform(y, family = "bcnPower")
#'
#' @return a transformed object of the same class as y
#'
#' @export
invert_initial_transform <- function(y, transformation, bc_params) {
  if(identical(transformation, "log")) {
    detransformed_y <- exp(y)
  } else if(identical(transformation, "box-cox")) {
    if(missing(bc_params)) {
      stop("bc_params must be provided if transformation is 'box-cox'")
    }
    detransformed_y <- invert_bcn_transform(b = y,
      lambda = bc_params$lambda,
      gamma = bc_params$gamma)
  } else if(identical(transformation, "none") ||
      identical(transformation, "forecast-box-cox")) {
    detransformed_y <- y
  } else {
    stop("Invalid transformation: must be one of 'box-cox', 'log', 'forecast-box-cox', or 'none'.")
  }

  return(detransformed_y)
}

invert_bcn_transform <- function(b, lambda, gamma) {
  ## Two steps: 1) undo box-cox 2) get y from z
  warning("Function invert_bcn_transform is untested!!")

  ## 1) undo box-cox: straightforward
  z <- (lambda * b + 1)^(1 / lambda)

  ## 2) get y from z: I couldn't find an exact formula, doing it numerically...
  y <- sapply(z, function(z_i) {
      temp <- optim(
        par = 0, # not a good default, could do better?  plot 0.5 * (y_i + (y_i^2 + gamma^2)^(0.5))
        fn = function(y_i) {
          (z_i - 0.5 * (y_i + (y_i^2 + gamma^2)^(0.5)))^2
        },
        method = "L-BFGS-B"
      )
      return(temp$par)
    })

  return(y)
}

#' Do first-order seasonal differencing (go from original time series values to
#' seasonally differenced time series values).
#'
#' @param y a univariate time series or numeric vector.
#' @param ts_frequency frequency of time series.  Must be provided if y is not
#'   of class "ts".  See the help for stats::ts for more.
#'
#' @return a seasonally differenced time series object (of class 'ts'),
#'   padded with leading NAs.
#'
#' @export
do_seasonal_difference <- function(y, ts_frequency) {
  differenced_y <- ts(c(rep(NA, ts_frequency),
    y[seq(from = ts_frequency + 1, to = length(y))] -
      y[seq(from = 1, to = length(y) - ts_frequency)]),
    frequency = ts_frequency)
  return(differenced_y)
}

#' Invert first-order seasonal differencing (go from seasonally differenced time
#' series values to original time series values).
#'
#' @param dy a first-order seasonally differenced univariate time series with
#'   values like y_{t} - y_{t - ts_frequency}
#' @param y a univariate time series or numeric vector with values like
#'   y_{t - ts_frequency}.
#' @param ts_frequency frequency of time series.  Must be provided if y is not
#'   of class "ts".  See the help for stats::ts for more.
#'
#' @details y may have longer length than dy.  It is assumed that dy "starts"
#'   one time index after y "ends": that is, if y is of length T then
#'   dy[1] = y[T + 1] - y[T + 1 - ts_frequency]
#'
#' @return a time series object (of class 'ts')
#'
#' @export
invert_seasonal_difference <- function(dy, y, ts_frequency) {
  return(dy + y[length(y) + seq_along(dy) - ts_frequency])
}

## deal with internal missing or infinite values in y that can
## result in simulated trajectories of all NAs if the model has a moving
## average component.  Here we do this by linear interpolation.
##
## Another (better?) solution would be to write a version of stats::filter
## that does the "right thing" with NAs if filter coefficients are 0 and then
## use that function in forecast:::myarima.sim
interpolate_and_clean_missing <- function(y) {
  if(any(is.na(y) | is.infinite(y)) && !is.na(tail(y, 1))) {
    ## drop leading NAs
    if(is.na(y[1] | is.infinite(y[1]))) {
      num_leading_nas <- rle(is.na(y) | is.infinite(y))$lengths[1]
      y <- y[- seq_len(num_leading_nas)]
    }

    ## interpolate internal NAs
    while(any(is.na(y) | is.infinite(y))) {
      na_rle <- rle(is.na(y) | is.infinite(y))
      na_run_start_ind <- na_rle$lengths[1] + 1
      na_run_end_ind <- na_run_start_ind + na_rle$lengths[2] - 1
      y[na_run_start_ind:na_run_end_ind] <-
        approx(
          x = c(na_run_start_ind - 1, na_run_end_ind + 1),
          y = y[c(na_run_start_ind - 1, na_run_end_ind + 1)],
          xout = na_run_start_ind:na_run_end_ind,
          method = "linear"
          )$y
    }
  }

  return(y)
}
