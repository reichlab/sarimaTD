#' sarimaTD interface via ForecastFramework
#'
#' Allows the user to use functionality provided by ForecastFramework with
#' sarimaTD models.
#'
#' @section Usage:
#' \preformatted{sarimaTD_model = sarimaTD_FF$new()
#'
#' sarimaTD_model$fit(data,
#'     period = 52,
#'     transformation = "box-cox",
#'     seasonal_difference = TRUE)
#'
#' sarimaTD_model$forecast(newdata,
#'     steps = 1,
#'     nsim = 1000)
#'
#' sarimaTD_model$data()
#'
#' sarimaTD_model$frequency()
#'
#' sarimaTD_model$nsim(nsim)
#'
#' sarimaTD_model$models()
#' }
#'
#' @section Arguments:
#' \code{data} A data set to fit sarimaTD models to, in a format compatible with
#'     ForecastFramework::IncidenceMatrix
#'
#' \code{frequency} frequency of time series.  Must be provided if y is not
#'     of class "ts".  See the help for stats::ts for more.
#'
#' \code{transformation} character specifying transformation type:
#'     "box-cox", "log", "forecast-box-cox", or "none".  See details of
#'     sarimaTD::fit_sarima for more.
#'
#' \code{seasonal_difference} boolean; take a seasonal difference before passing
#'     to auto.arima?
#'
#' \code{newdata} new data from which to predict forward, in a format compatible
#'     with ForecastFramework::IncidenceMatrix
#'
#' \code{steps} If \code{TRUE}, force the Python process to terminate
#'   using a system call.
#'
#' \code{nsim} Number of simulations to generate during forecasting
#'
#' @section Methods:
#' \code{$new()} Initialize a sarimaTD model object.
#'
#' \code{$fit()} Fit a set of models to the provided data, one model per
#'   location.
#'
#' \code{$forecast()} Generate forecasts for each location.
#'
#' \code{$data()} Get the data set used for the model fit.
#'
#' \code{$nsim()} Get or set the number of simulations to generate during
#'     forecasting.
#'
#' \code{$frequency()} Get the time series frequency of the data used for the
#'     model fit.
#'
#' \code{$transformation()} Get the transformation performed before fitting the
#'     model.
#'
#' \code{$seasonal_difference()} Get the indicator of whether or not first-order
#'     seasonal differencing is performed before fitting the model.
#'
#' \code{$models()} Get a list of sarimaTD model fit objects, one for each
#'     location in the training data.
#'
#' @name sarimaTD_FF
#' @examples
#' vignette(package="sarimaTD", topic="sarimaTD_FF")
NULL

#' @export
sarimaTD_FF <- R6::R6Class(
    inherit = ForecastFramework::ForecastModel,
    private = list(
        .data = NULL, ## data, as an IncidenceMatrix
        .frequency = integer(0), ## seasonal period for data
        .transformation = character(0), ## transformation to use
        .seasonal_difference = logical(0), ## seasonal differencing
        .nsim = integer(0), ## how many simulated trajectories for forecast
        .models = list() ## models that are fit separately for each location
    ),
    public = list(
        ## data will be MatrixData
        fit = function(data) {
            if("fit" %in% private$.debug){browser()}
            ## store data and check to make sure it's the right class
            private$.data <- IncidenceMatrix$new(data)
            
            ## for each location/row
            for (row_idx in 1:private$.data$nrow) {
                ## create a y vector with incidence at time t
                y <- private$.data$subset(rows = row_idx, mutate = FALSE)
                
                ## convert y vector to time series data type
                y_ts <- ts(as.vector(y$mat), frequency = private$.frequency)
                
                ## fit sarimaTD with 'fit_sarima()' from sarimaTD package
                ## if requested, fit_sarima() performs transformation and
                ## seasonal differencing
                private$.models[[row_idx]] <- fit_sarima(y = y_ts,
                    transformation = private$.transformation,
                    seasonal_difference = private$.seasonal_difference)
            }
        },
        forecast = function(
                newdata = private$.data,
                steps = 1) {
            ## include for debugging
            if("forecast" %in% private$.debug){browser()}
            
            ## number of models (locations) to forecast
            nmodels <- length(private$.models)
            
            ## define an array to store the simulated forecasts
            sim_forecasts <- array(dim = c(nmodels, steps, private$.nsim))
            dimnames(sim_forecasts) <- list(newdata$rnames, 1:steps, NULL)
            
            ## iterate through each province and forecast with simulate.sarimaTD
            for(model_idx in 1:length(private$.models)) {
                tmp_arima <-  simulate(object = private$.models[[model_idx]],
                    nsim = private$.nsim,
                    seed = 1,
                    newdata = as.vector(newdata$mat[model_idx,]),
                    h = steps
                )
                ## transpose simulate() output to be consistent with
                ## ForecastFramework
                tmp_arima <- t(tmp_arima)
                sim_forecasts[model_idx, , ] <- tmp_arima
            }
            private$output <- SimulatedIncidenceMatrix$new(sim_forecasts)
            return(IncidenceForecast$new(private$output,
                forecastTimes = rep(TRUE, steps)))
        },
        initialize = function(
                nsim = 1000,
                frequency = 52,
                transformation = "box-cox",
                seasonal_difference = TRUE) {
            ## this code is run during sarimaTD_FF$new()
            private$.nsim <- nsim
            private$.frequency <- frequency
            private$.transformation <- transformation
            private$.seasonal_difference <- seasonal_difference
        },
        predict = function(newdata) {
            stop("predict method has not been written.")
        }
    ),
    active = list(
        ## This list determines how users can access pieces of the model object
        data = function(value) {
            if(!missing(value))
                stop("Writing directly to the data is not allowed.")
            return(private$.data)
        },
        nsim = function(nsim) {
            if(!missing(nsim)) {
                if(!is.integer(nsim)) {
                    stop("nsim must be an integer.")
                }
                
                private$.nsim <- nsim
                return(invisible(private$.nsim))
            } else {
                return(private$.nsim)
            }
        },
        frequency = function(value) {
            if(!missing(value))
                stop("The time series frequency can only be set via a call to the new function.")
            return(private$.frequency)
        },
        transformation = function(value) {
            if(!missing(value))
                stop("The transformation can only be set via a call to the new function.")
            return(private$.transformation)
        },
        seasonal_difference = function(value) {
            if(!missing(value))
                stop("The time series differencing can only be set via a call to the new function.")
            return(private$.seasonal_difference)
        },
        models = function(value) {
            if(!missing(value))
                stop("Writing directly to the models is not allowed.")
            return(private$.models)
        }
    )
)
