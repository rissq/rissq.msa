#' Generic for Anova MSA Study.
#' @name anovaMSA
#' @export
setGeneric (
  name = "anovaMSA",
  def = function(object){ standardGeneric("anovaMSA") }
)

#' Generic for Gage rar MSA Study.
#' @name rar
#' @export
setGeneric (
  name = "rar",
  def = function(object){ standardGeneric("rar") }
)

#' Generic for plotting Component of Variation Chart
#' @name plotComponentOfVariation
#' @export
setGeneric (
  name = "plotComponentOfVariation",
  def = function(object){ standardGeneric("plotComponentOfVariation") }
)

#' Generic for plotting Variable by Part Chart
#' @name plotVariableByPart
#' @export
setGeneric (
  name = "plotVariableByPart",
  def = function(object){ standardGeneric("plotVariableByPart") }
)

#' Generic for plotting Variable by Appraiser Chart
#' @name plotVariableByAppraiser
#' @export
setGeneric (
  name = "plotVariableByAppraiser",
  def = function(object){ standardGeneric("plotVariableByAppraiser") }
)

#' Generic for plotting Interaction Chart
#' @name plotInteraction
#' @export
setGeneric (
  name = "plotInteraction",
  def = function(object){ standardGeneric("plotInteraction") }
)

#' Generic for plotting Mean Chart
#' @name plotMean
#' @export
setGeneric (
  name = "plotMean",
  def = function(object){ standardGeneric("plotMean") }
)

#' Generic for plotting Range Chart
#' @name plotRange
#' @export
setGeneric (
  name = "plotRange",
  def = function(object){ standardGeneric("plotRange") }
)
