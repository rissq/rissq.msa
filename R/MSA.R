#' constructor for MSA
#'
#' This is the constructor.
#' @name initializeMSA
#' @export
setMethod("initialize",
          signature = signature(.Object = "MSA"),
          function(.Object, ..., id, name, description, pro, characteristic, data, tolerance, sigma, alphaLim) {

            if(missing(tolerance) || is.na(tolerance)) {
              tolerance = characteristic@U - characteristic@L
            }

            if(missing(sigma) || is.na(sigma)) {
              sigma = 6
            }

            if(missing(alphaLim) || is.na(alphaLim)) {
              alphaLim = 0.05
            }

            callNextMethod(.Object, ..., id = id, name = name, description = description, pro = pro, characteristic = characteristic, data = data, tolerance = tolerance, sigma = sigma, alphaLim = alphaLim)
          })

setMethod("show", "MSA", function(object) {

  callNextMethod(object)
})

#' Summarize a MSA object
#' @name summary
#' @export
setMethod("summary", "MSA", function(object, ...) {

  cat("Analysis, ", object@name)
  cat("\n", "ID, ", object@id)
  cat("\n", "Description, ", object@description)
  cat("\n")
  cat("\n", "Levels part, ", object@lvlPart)
  cat("\n", "Levels appraiser, ", object@lvlAppr)
  cat("\n", "N, ", object@n)
  cat("\n")
  cat("\n", "Complete model:\n")
  print(object@anova[[1]])
  cat("\n", "Alpha for removing interaction, ", object@alphaLim)
  cat("\n")
  cat("\n", "Gage R&R:\n")
  print(object@varianceComponents)
  cat("\n", "Number of distinct categories, ", object@numberCategories)
  cat("\n")

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

            return(plot(g2))
          })
