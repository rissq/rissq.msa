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

#' Plot rar resume
#' @name plotRar
#' @export
setMethod("plotRar",
          signature = signature(object = "MSA"),
          function(object){

            lay <- rbind(c(1,2),
                         c(3,4),
                         c(5,NA))

            g <- grid.arrange(
              ggplotComponentOfVariationChart(object),
              ggplotVariableByPartChart(object),
              ggplotRangeChart(object, gridLayout = TRUE),
              ggplotVariableByAppraiserChart(object),
              ggplotMeanChart(object, gridLayout = TRUE),
              layout_matrix = lay,
              top = object@name,
              bottom = object@description
            )

            g2 <- cowplot::ggdraw(g) +
              theme(plot.background = element_rect(fill="white", color = NA))

            plot(g2)
          })
