simulate.Arima <- function (object, nsim = length(object$x), seed = NULL, xreg = NULL, 
          future = TRUE, bootstrap = FALSE, innov = NULL, lambda = object$lambda, 
          ...) 
{
  if (object$arma[7] < 0) {
    stop("Value for seasonal difference is < 0. Must be >= 0")
  }
  else if ((sum(object$arma[c(3, 4, 7)]) > 0) && (object$arma[5] < 
                                                  2)) {
    stop("Invalid value for seasonal period")
  }
  if (!is.null(xreg)) {
    xreg <- as.matrix(xreg)
    nsim <- nrow(xreg)
  }
  if (is.null(innov)) {
    if (!exists(".Random.seed", envir = .GlobalEnv)) {
      runif(1)
    }
    if (is.null(seed)) {
      RNGstate <- .Random.seed
    }
    else {
      R.seed <- .Random.seed
      set.seed(seed)
      RNGstate <- structure(seed, kind = as.list(RNGkind()))
      on.exit(assign(".Random.seed", R.seed, envir = .GlobalEnv))
    }
  }
  else {
    nsim <- length(innov)
  }
  flag.s.arma <- (sum(object$arma[c(3, 4)]) > 0)
  if (sum(object$arma[c(3, 4, 7)]) > 0) {
    if (sum(object$model$phi) == 0) {
      ar <- NULL
    }
    else {
      ar <- as.double(object$model$phi)
    }
    if (sum(object$model$theta) == 0) {
      ma <- NULL
    }
    else {
      ma <- as.double(object$model$theta)
    }
    order <- c(length(ar), object$arma[6], length(ma))
    if (future) {
      model <- list(order = order, ar = ar, ma = ma, sd = sqrt(object$sigma2), 
                    residuals = residuals(object), seasonal.difference = object$arma[7], 
                    seasonal.period = object$arma[5], flag.seasonal.arma = flag.s.arma, 
                    seasonal.order = object$arma[c(3, 7, 4)])
    }
    else {
      model <- list(order = order, ar = ar, ma = ma, sd = sqrt(object$sigma2), 
                    residuals = residuals(object))
    }
    flag.seasonal.diff <- (object$arma[7] > 0)
  }
  else {
    order <- object$arma[c(1, 6, 2)]
    if (order[1] > 0) {
      ar <- object$model$phi[1:order[1]]
    }
    else {
      ar <- NULL
    }
    if (order[3] > 0) {
      ma <- object$model$theta[1:order[3]]
    }
    else {
      ma <- NULL
    }
    if (object$arma[2] != length(ma)) {
      stop("MA length wrong")
    }
    else if (object$arma[1] != length(ar)) {
      stop("AR length wrong")
    }
    if (future) {
      model <- list(order = object$arma[c(1, 6, 2)], ar = ar, 
                    ma = ma, sd = sqrt(object$sigma2), residuals = residuals(object), 
                    seasonal.difference = 0, flag.seasonal.arma = flag.s.arma, 
                    seasonal.order = c(0, 0, 0), seasonal.period = 1)
    }
    else {
      model <- list(order = object$arma[c(1, 6, 2)], ar = ar, 
                    ma = ma, sd = sqrt(object$sigma2), residuals = residuals(object))
    }
    flag.seasonal.diff <- FALSE
  }
  x <- object$x <- forecast::getResponse(object)
  if (is.null(tsp(x))) {
    x <- ts(x, frequency = 1, start = 1)
  }
  if (!is.null(lambda)) {
    x <- BoxCox(x, lambda)
    lambda <- attr(x, "lambda")
  }
  n <- length(x)
  if (bootstrap) {
    res <- na.omit(c(model$residuals) - mean(model$residuals, 
                                             na.rm = TRUE))
    e <- sample(res, nsim, replace = TRUE)
  }
  else if (is.null(innov)) {
    e <- rnorm(nsim, 0, model$sd)
  }
  else if (length(innov) == nsim) {
    e <- innov
  }
  else {
    stop("Length of innov must be equal to nsim")
  }
  use.drift <- is.element("drift", names(object$coef))
  usexreg <- (!is.null(xreg) | use.drift)
  xm <- oldxm <- 0
  if (!is.null(xreg)) {
    xreg <- as.matrix(xreg)
    if (nrow(xreg) < nsim) {
      stop("Not enough rows in xreg")
    }
    else {
      xreg <- xreg[1:nsim, ]
    }
  }
  if (use.drift) {
    if (NCOL(xreg) == 1) {
      xreg <- NULL
    }
    else {
      xreg <- xreg[, !is.element(colnames(xreg), "drift"), 
                   drop = FALSE]
    }
    dft <- as.matrix(1:nsim)
    if (future) {
      dft <- dft + n
    }
    xreg <- cbind(drift = dft, xreg)
  }
  narma <- sum(object$arma[1L:4L])
  if (length(object$coef) > narma) {
    if (names(object$coef)[narma + 1L] == "intercept") {
      xreg <- cbind(intercept = rep(1, nsim), xreg)
      object$xreg <- cbind(intercept = rep(1, n), object$xreg)
    }
    if (!is.null(xreg)) {
      xm <- if (narma == 0) {
        drop(as.matrix(xreg) %*% object$coef)
      }
      else {
        drop(as.matrix(xreg) %*% object$coef[-(1L:narma)])
      }
      oldxm <- if (narma == 0) {
        drop(as.matrix(object$xreg) %*% object$coef)
      }
      else {
        drop(as.matrix(object$xreg) %*% object$coef[-(1L:narma)])
      }
    }
  }
  if (future) {
    sim <- myarima.sim(model, nsim, x - oldxm, e = e) + xm
  }
  else {
    if (flag.seasonal.diff) {
      zeros <- object$arma[5] * object$arma[7]
      sim <- arima.sim(model, nsim, innov = e)
      sim <- diffinv(sim, lag = object$arma[5], differences = object$arma[7])[-(1:zeros)]
      sim <- ts(tail(sim, nsim) + xm)
    }
    else {
      sim <- ts(tail(arima.sim(model, nsim, innov = e), 
                     nsim) + xm)
    }
    if (nsim == n) 
      tsp(sim) <- tsp(x)
    if (model$order[2] > 0 || flag.seasonal.diff) {
      sim <- sim - sim[1] + x[1]
    }
  }
  if (!is.null(lambda)) {
    sim <- InvBoxCox(sim, lambda)
  }
  return(sim)
}


myarima.sim <- function (model, n, x, e, ...) {
  start.innov <- residuals(model)
  innov <- e
  data <- x
  first.nonmiss <- which(!is.na(x))[1]
  if (first.nonmiss > 1) {
    tsp.x <- tsp(x)
    start.x <- tsp.x[1] + (first.nonmiss - 1)/tsp.x[3]
    x <- window(x, start = start.x)
    start.innov <- window(start.innov, start = start.x)
  }
  model$x <- x
  n.start <- length(x)
  x <- ts(c(start.innov, innov), start = 1 - n.start, frequency = model$seasonal.period)
  flag.noadjust <- FALSE
  if (is.null(tsp(data))) {
    data <- ts(data, frequency = 1, start = 1)
  }
  if (!is.list(model)) {
    stop("'model' must be list")
  }
  if (n <= 0L) {
    stop("'n' must be strictly positive")
  }
  p <- length(model$ar)
  q <- length(model$ma)
  d <- 0
  D <- model$seasonal.difference
  m <- model$seasonal.period
  if (!is.null(ord <- model$order)) {
    if (length(ord) != 3L) {
      stop("'model$order' must be of length 3")
    }
    if (p != ord[1L]) {
      stop("inconsistent specification of 'ar' order")
    }
    if (q != ord[3L]) {
      stop("inconsistent specification of 'ma' order")
    }
    d <- ord[2L]
    if (d != round(d) || d < 0) {
      stop("number of differences must be a positive integer")
    }
  }
  # COMMENTING OUT CHECK FOR STATIONARITY.
  # In some versions of R polyroot seems to return incorrect results that result
  # in a failing check for stationary models.
  #if (p) {
  #  minroots <- min(Mod(polyroot(c(1, -model$ar))))
  #  if (minroots <= 1) {
  #    stop("'ar' part of model is not stationary")
  #  }
  #}
  if (length(model$ma)) {
    x <- stats::filter(x, c(1, model$ma), method = "convolution", 
                       sides = 1L)
    x[seq_along(model$ma)] <- 0
  }
  len.ar <- length(model$ar)
  if (length(model$ar) && (len.ar <= length(data))) {
    if ((D != 0) && (d != 0)) {
      diff.data <- diff(data, lag = 1, differences = d)
      diff.data <- diff(diff.data, lag = m, differences = D)
    }
    else if ((D != 0) && (d == 0)) {
      diff.data <- diff(data, lag = model$seasonal.period, 
                        differences = D)
    }
    else if ((D == 0) && (d != 0)) {
      diff.data <- diff(data, lag = 1, differences = d)
    }
    else {
      diff.data <- data
    }
    x.new.innovations <- x[(length(start.innov) + 1):length(x)]
    x.with.data <- c(diff.data, x.new.innovations)
    for (i in (length(diff.data) + 1):length(x.with.data)) {
      lagged.x.values <- x.with.data[(i - len.ar):(i - 
                                                     1)]
      ar.coefficients <- model$ar[length(model$ar):1]
      sum.multiplied.x <- sum((lagged.x.values * ar.coefficients)[abs(ar.coefficients) > 
                                                                    .Machine$double.eps])
      x.with.data[i] <- x.with.data[i] + sum.multiplied.x
    }
    x.end <- x.with.data[(length(diff.data) + 1):length(x.with.data)]
    x <- ts(x.end, start = 1, frequency = model$seasonal.period)
    flag.noadjust <- TRUE
  }
  else if (length(model$ar)) {
    x <- stats::filter(x, model$ar, method = "recursive")
  }
  if ((d == 0) && (D == 0) && (flag.noadjust == FALSE)) {
    if (n.start >= 20) {
      xdiff <- (model$x - x[1:n.start])[n.start - (19:0)]
    }
    else {
      xdiff <- model$x - x[1:n.start]
    }
    if (all(sign(xdiff) == 1) || all(sign(xdiff) == -1)) {
      xdiff <- xdiff[length(xdiff)]
    }
    else {
      xdiff <- mean(xdiff)
    }
    x <- x + xdiff
  }
  if ((n.start > 0) && (flag.noadjust == FALSE)) {
    x <- x[-(1:n.start)]
  }
  if ((D > 0) && (d == 0)) {
    i <- length(data) - D * m + 1
    seasonal.xi <- data[i:length(data)]
    length.s.xi <- length(seasonal.xi)
    x <- diffinv(x, lag = m, differences = D, xi = seasonal.xi)[-(1:length.s.xi)]
  }
  else if ((d > 0) && (D == 0)) {
    x <- diffinv(x, differences = d, xi = data[length(data) - 
                                                 (d:1) + 1])[-(1:d)]
  }
  else if ((d > 0) && (D > 0)) {
    delta.four <- diff(data, lag = m, differences = D)
    regular.xi <- delta.four[(length(delta.four) - D):length(delta.four)]
    x <- diffinv(x, differences = d, xi = regular.xi[length(regular.xi) - 
                                                       (d:1) + 1])[-(1:d)]
    i <- length(data) - D * m + 1
    seasonal.xi <- data[i:length(data)]
    length.s.xi <- length(seasonal.xi)
    x <- diffinv(x, lag = m, differences = D, xi = seasonal.xi)
    x <- x[-(1:length.s.xi)]
  }
  x <- ts(x[1:n], frequency = frequency(data), start = tsp(data)[2] + 
            1/tsp(data)[3])
  return(x)
}
