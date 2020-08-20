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
#' @name plotComponentOfVariationChart
#' @export
setGeneric (
  name = "plotComponentOfVariationChart",
  def = function(object){ standardGeneric("plotComponentOfVariationChart") }
)

#' Generic for plotting Variable by Part Chart
#' @name plotVariableByPartChart
#' @export
setGeneric (
  name = "plotVariableByPartChart",
  def = function(object){ standardGeneric("plotVariableByPartChart") }
)

#' Generic for plotting Variable by Appraiser Chart
#' @name plotVariableByAppraiserChart
#' @export
setGeneric (
  name = "plotVariableByAppraiserChart",
  def = function(object){ standardGeneric("plotVariableByAppraiserChart") }
)

#' Generic for plotting Interaction Chart
#' @name plotInteractionChart
#' @export
setGeneric (
  name = "plotInteractionChart",
  def = function(object){ standardGeneric("plotInteractionChart") }
)

#' Generic for plotting Mean Chart
#' @name plotMeanChart
#' @export
setGeneric (
  name = "plotMeanChart",
  def = function(object){ standardGeneric("plotMeanChart") }
)

#' Generic for plotting Range Chart
#' @name plotRangeChart
#' @export
setGeneric (
  name = "plotRangeChart",
  def = function(object){ standardGeneric("plotRangeChart") }
)
