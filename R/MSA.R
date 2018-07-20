#' constructor for MSA
#'
#' This is the constructor.
#' @name initializeMSA
#' @export
setMethod("initialize",
          signature = signature(.Object = "MSA"),
          function(.Object, ..., id, name, description, pro, variable, part, appraiser, data, usl = NA_real_, lsl = NA_real_, tolerance = NA_real_, sigma = 6, alphaLim = 0.05, errorTerm = "interaction", digits = 4) {

            part <- as.factor(part)
            appraiser <- as.factor(appraiser)

            callNextMethod(.Object, ..., id = id, name = name, description = description, pro = pro, variable = variable, part = part, appraiser = appraiser, data = data, usl = usl, lsl = lsl, tolerance = tolerance, sigma = sigma, alphaLim = alphaLim, errorTerm = errorTerm, digits = digits)
          })
