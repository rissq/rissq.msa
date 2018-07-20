#' constructor for NestedMSA
#'
#' This is the constructor.
#' @name initializeNestedMSA
#' @export
setMethod("initialize",
          signature = signature(.Object = "NestedMSA"),
          function(.Object, ..., id, name, description, pro = NA, variable, part, appraiser, data, usl = NA_real_, lsl = NA_real_, tolerance = NA_real_, sigma = 6, alphaLim = 0.05, errorTerm = "interaction", digits = 4) {

            lvlPart = nlevels(data[[part]])
            n = nrow(data)/lvlPart

            callNextMethod(.Object, ..., id = id, name = name, description = description, pro = pro, variable = variable, part = part, appraiser = appraiser, data = data, usl = usl, lsl = lsl, sigma = sigma, alphaLim = alphaLim, errorTerm = errorTerm, digits = digits, lvlPart = lvlPart, n = n)
          })

#' Anova study for Nested MSA
#' @name anovaStudy
#' @export
setMethod("anovaStudy",
          signature = signature(object = "NestedMSA"),
          function(object){
            ## Complete model (with interaction)
            modelf <- as.formula(paste(object@variable, "~", object@appraiser, "/", object@part))

            model <- aov(modelf, data = object@pro@data)
            modelm <- summary(model)

            if (object@errorTerm == "interaction"){
              modelm[[1]][1, 4] <- modelm[[1]][1, 3]/modelm[[1]][2, 3]
              modelm[[1]][1, 5] <- pf(modelm[[1]][1, 4],
                                      modelm[[1]][1, 1],
                                      modelm[[1]][3, 1], lower.tail = FALSE)
            }

            rownames(modelm[[1]])[3] <- "Repeatability"

            modelm[[1]] <- rbind(modelm[[1]],
                                 c(colSums(modelm[[1]][, 1:2]), rep(NA, 3)))

            object@anova <- list(modelm)
            return(object)
          })
