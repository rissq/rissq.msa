#' @export .MSA
#' @exportClass MSA
.MSA <- setClass("MSA",
                 slots = c(
                   # variable = "character",
                   # part = "factor",
                   # appraiser = "factor",
                   # usl = "numeric",
                   # lsl = "numeric",
                   tolerance = "numeric",
                   sigma = "numeric",
                   alphaLim = "numeric",
                   digits = "numeric",
                   lvlPart = "numeric",
                   lvlAppr = "numeric",
                   n = "numeric",
                   anova = "list",
                   varianceComponents = "matrix",
                   numberCategories = "numeric"
                 ),
                 contains = "Analysis"
)

#' @export .NestedMSA
#' @exportClass NestedMSA
.NestedMSA <- setClass("NestedMSA",
                 contains = "MSA"
)

#' @export .CrossedMSA
#' @exportClass CrossedMSA
.CrossedMSA <- setClass("CrossedMSA",
                 contains = "MSA"
)
