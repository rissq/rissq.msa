#' constructor for NestedMSA
#'
#' This is the constructor.
#' @name initializeNestedMSA
#' @export
setMethod("initialize",
          signature = signature(.Object = "NestedMSA"),
          function(.Object, ..., id, name, description, pro, characteristic, data, tolerance, sigma, alphaLim) {

            headers <- names(data@data)

            part <- as.factor(headers[1])
            appraiser <- as.factor(headers[2])

            data@data[[part]] <- factor(data@data[[part]])
            data@data[[appraiser]] <- factor(data@data[[appraiser]])

            #Number of parts levels and replicates
            lvlPart = nlevels(data@data[[part]])
            lvlAppr = nlevels(data@data[[appraiser]])
            n = nrow(data@data)/lvlPart

            .Object <- callNextMethod(.Object, ..., id = id, name = name, description = description, pro = pro, characteristic = characteristic, data = data, tolerance = tolerance, sigma = sigma, alphaLim = alphaLim, lvlPart = lvlPart, lvlAppr = lvlAppr, n = n)

            #If data is rady calculations are made on the initialization method call
            .Object <- anovaMSA(.Object)
            .Object <- rar(.Object)
          })

#' Anova study for Nested MSA
#' @name anovaMSA
#' @export
setMethod("anovaMSA",
          signature = signature(object = "NestedMSA"),
          function(object){
            headers <- names(object@data@data)

            part <- headers[1]
            appraiser <- headers[2]
            variable <- headers[3]

            ## Complete model (with interaction)
            modelf <- as.formula(paste(variable, "~", appraiser, "/", part))

            model <- aov(modelf, data = object@data@data)
            modelm <- summary(model)

            rownames(modelm[[1]])[3] <- "REPEATIBILITY"

            #Total row for Df and SumSq
            modelm[[1]] <- rbind(modelm[[1]], c(colSums(modelm[[1]][, 1:2]), rep(NA, 3)))

            rownames(modelm[[1]])[4] <- "TOTAL"

            object@anova <- list(modelm[[1]])

            show(object)
            return(object)
          })


#' Gage rar for Nested MSA
#' @name rar
#' @export
setMethod("rar",
          signature = signature(object = "NestedMSA"),
          function(object){

            ## if anova is not stored we calcute it
            if(!length(object@anova)) {
              object <- anovaMSA(object)
            }

            object@varianceComponents = matrix(ncol = 6, nrow = 5)

            rownames(object@varianceComponents) <- c("Total Gage R&R", "  Repeatability", "  Reproducibility", "Part-To-Part", "Total Variation")

            colnames(object@varianceComponents) <- c("VarComp", "%Contrib", "StdDev", "StudyVar", "%StudyVar", "%Tolerance")

            #Variance Components Table
            #Repeatibility
            object@varianceComponents[2, 1] <- object@anova[[1]][3, 3]
            #Total reproducibility
            object@varianceComponents[3, 1] <- max(c((object@anova[[1]][1, 1] * (object@anova[[1]][1, 3] - object@anova[[1]][2, 3]) / (object@lvlPart * object@n)), 0))
            #Part to part
            object@varianceComponents[4, 1] <- max(c((object@anova[[1]][2, 3] - object@anova[[1]][3, 3]) / object@n, 0))
            #Totat Gage rar
            object@varianceComponents[1, 1] <- object@varianceComponents[2, 1] + object@varianceComponents[3, 1]
            #Total variation
            object@varianceComponents[5, 1] <- object@varianceComponents[1, 1] + object@varianceComponents[4, 1]
            #%Contrib
            object@varianceComponents[, 2] <- round(100 * (object@varianceComponents[, 1]/object@varianceComponents[5, 1]), 2)
            #Standard Deviation
            object@varianceComponents[, 3] <- sqrt(object@varianceComponents[, 1])
            #Study Variation edited from 5.15 to variable
            object@varianceComponents[, 4] <- object@varianceComponents[, 3] * object@sigma
            #
            object@varianceComponents[, 5] <- round(100 * (object@varianceComponents[, 3]/object@varianceComponents[5, 3]), 2)
            #
            object@varianceComponents[, 6] <- round(100 * (object@varianceComponents[, 4]/(object@tolerance)), 2)

            #Number of distinct categories
            object@numberCategories <- max(c(1, floor((object@varianceComponents[4, 4]/object@varianceComponents[1, 4])*1.41)))

            show(object)
            return(object)
          })

#' MSA resume Chart
#' @name plot
#' @export
setMethod("plot",
          signature = signature(object = "NestedMSA"),
          function(object){

            show(object)
            return(object)
          })

#' Components of Variation Chart
#' @name plotComponentOfVariation
#' @export
setMethod("plotComponentOfVariation",
          signature = signature(object = "NestedMSA"),
          function(object){

            show(object)
            return(object)
          })

#' Variable by Part
#' @name plotVariableByPart
#' @export
setMethod("plotVariableByPart",
          signature = signature(object = "NestedMSA"),
          function(object){

            show(object)
            return(object)
          })

#' Variable by Appraiser
#' @name plotVariableByAppraiser
#' @export
setMethod("plotVariableByAppraiser",
          signature = signature(object = "NestedMSA"),
          function(object){

            show(object)
            return(object)
          })

#' Control Chart
#' @name plotComponentOfVariation
#' @export
setMethod("plotControl",
          signature = signature(object = "NestedMSA"),
          function(object){

            show(object)
            return(object)
          })

#' Mean Chart
#' @name plotMean
#' @export
setMethod("plotMean",
          signature = signature(object = "NestedMSA"),
          function(object){

            show(object)
            return(object)
          })

#' Range Chart
#' @name plotRange
#' @export
setMethod("plotRange",
          signature = signature(object = "NestedMSA"),
          function(object){

            show(object)
            return(object)
          })
