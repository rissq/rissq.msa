#' constructor for MSA
#'
#' This is the constructor.
#' @name initializeMSA
#' @export
setMethod("initialize",
          signature = signature(.Object = "MSA"),
          function(.Object, ..., id, name, description, pro, characteristic, data, tolerance, sigma, alphaLim) {

            if(missing(tolerance)) {
              tolerance = characteristic@U - characteristic@L
            }

            if(missing(sigma)) {
              sigma = 6
            }

            if(missing(alphaLim)) {
              alphaLim = 0.05
            }

            callNextMethod(.Object, ..., id = id, name = name, description = description, pro = pro, characteristic = characteristic, data = data, tolerance = tolerance, sigma = sigma, alphaLim = alphaLim)
          })


setMethod("show", "MSA", function(object) {

  callNextMethod(object)
})
