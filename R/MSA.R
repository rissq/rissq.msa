#' constructor for MSA
#'
#' This is the constructor.
#' @name initializeMSA
#' @export
setMethod("initialize",
          signature = signature(.Object = "MSA"),
          function(.Object, ..., id, name, description, pro, variable, part, appraiser, data, usl, lsl, tolerance, sigma, alphaLim, digits) {

            if(missing(usl)) {
              usl = NA_real_
            }

            if(missing(lsl)) {
              lsl = NA_real_
            }

            if(missing(tolerance)) {
              tolerance = usl - lsl
            }

            if(missing(sigma)) {
              sigma = 6
            }

            if(missing(alphaLim)) {
              alphaLim = 0.05
            }

            if(missing(digits)) {
              digits = 4
            }

            callNextMethod(.Object, ..., id = id, name = name, description = description, pro = pro, variable = variable, part = part, appraiser = appraiser, data = data, usl = usl, lsl = lsl, tolerance = tolerance, sigma = sigma, alphaLim = alphaLim, digits = digits)
          })


setMethod("show", "MSA", function(object) {

  callNextMethod(object)
})
