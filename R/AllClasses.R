#' @export .MSA
#' @exportClass MSA
.MSA <- setClass("MSA",
                 slots = c(
                   variable = "character",
                   part = "factor",
                   appraiser = "factor",
                   data = "data.frame",
                   usl = "numeric",
                   lsl = "numeric",
                   tolerance = "numeric",
                   sigma = "numeric",
                   alphaLim = "numeric",
                   errorTerm = "character",
                   digits = "numeric",
                   lvlPart = "numeric",
                   lvlAppr = "numeric",
                   n = "numeric",
                   anova = "list"
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
