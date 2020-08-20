#' @export .MSA
#' @exportClass MSA
.MSA <- setClass("MSA",
                 slots = c(
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

#' @export .BaseMSA
#' @exportClass BaseMSA
.BaseMSA <- setClass("BaseMSA",
                        contains = "MSA"
)

#' @export .NestedMSA
#' @exportClass NestedMSA
.NestedMSA <- setClass("NestedMSA",
                        contains = "MSA"
)

#' @export .CrossedMSA
#' @exportClass CrossedMSA
.CrossedMSA <- setClass("CrossedMSA",
                          slots = c(
                            anovaReduced = "list"
                          ),
                          contains = "MSA"
)
