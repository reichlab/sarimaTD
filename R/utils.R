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
#    transformed_y <- car::bcnPower(
#      U = y,
#      lambda = bc_params$lambda,
#      gamma = bc_params$gamma)
    transformed_y <- car::bcPower(
      U = y + bc_params$gamma,
      lambda = bc_params$lambda)
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
#    detransformed_y <- invert_bcn_transform(b = y,
#      lambda = bc_params$lambda,
#      gamma = bc_params$gamma)
    detransformed_y <- invert_bc_transform(b = y,
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

#' Invert a Box-Cox (with negatives allowed) transformation.  See car::bcnPower
#' for the original transformation.
#'
#' @param b a univariate numeric vector.
#' @param lambda exponent for Box-Cox transformation
#' @param gamma offset for Box-Cox transformation
#'
#' @details Undoing the Box-Cox step is straightforward, but I could not find a
#'   closed-form way to undo the second step (but I may just be missing it).
#'   Currently, just using optim to minimize square difference from z and
#'   transformed y, as a function of y.  There must be something faster?
#'
#' @return a transformed object of the same class as y
#'
#' @export
invert_bcn_transform <- function(b, lambda, gamma) {
  ## Two steps: 1) undo box-cox 2) get y from z

  ## 1) undo box-cox: straightforward
  if(abs(lambda) <= 1e-10) {
    z <- exp(b)
  } else {
    z <- (lambda * b + 1)^(1 / lambda)
  }

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

#' Invert a Box-Cox transformation.  See car::bcPower for the original
#' transformation.
#'
#' @param b a univariate numeric vector.
#' @param lambda exponent for Box-Cox transformation
#' @param gamma offset for Box-Cox transformation
#'
#' @details Undoing the Box-Cox step is straightforward, then subtract gamma.
#'
#' @return a transformed object of the same class as y
#'
#' @export
invert_bc_transform <- function(b, lambda, gamma) {
  ## Two steps: 1) undo box-cox 2) get y from z

  ## 1) undo box-cox: straightforward
  if(abs(lambda) <= 1e-10) {
    z <- exp(b)
  } else {
    z <- (lambda * b + 1)^(1 / lambda)
  }

  ## 2) Subtract gamma
  return(z - gamma)
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
  return(ts(dy + y[length(y) + seq_along(dy) - ts_frequency],
    freq = ts_frequency))
}

#' Remove leading values that are infinite or missing, and replace all internal
#' sequences of infinite or missing values by linearly interpolating between
#' the surrounding observed (finite) values.  Used in prediction methods
#'
#' @param y a univariate time series or numeric vector
#' 
#' @return a vector of the same class as y
#'
#' @export
interpolate_and_clean_missing <- function(y) {
  if(any(is.na(y) | is.infinite(y)) && !is.na(tail(y, 1))) {
    ## drop leading NAs
    if(is.na(y[1]) | is.infinite(y[1])) {
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
