#' constructor for NestedMSA
#'
#' This is the constructor.
#' @name initializeNestedMSA
#' @export
setMethod("initialize",
          signature = signature(.Object = "NestedMSA"),
          function(.Object, ..., id, name, description, pro, variable, part, appraiser, data, usl, lsl, tolerance, sigma, alphaLim, digits) {

            #Checks missing names for columns
            if(missing(part) || missing(appraiser) || missing(variable)) {
              stop("Part, appraiser and measured variable must be specified.")
            }

            #Set all names to upper case to avoid problems when introducing column names for variables manually
            data <- setNames(data, toupper(names(data)))
            part <- toupper(part)
            appraiser <- toupper(appraiser)
            variable <- toupper(variable)

            #Checks if names are correctly introduced
            if (is.data.frame(data)){
              #Check for result variable
              if (!(variable %in% names(data))) {
                stop(variable, " is not a valid column name for ", deparse(substitute(data)))
              }

              #Checks for part, and convert column values to factors
              if (part %in% names(data)) {
                data[[part]] <- factor(data[[part]])
              } else{
                stop(part, " is not a valid column name for", deparse(substitute(data)))
              }

              #Checks for appraiser, and convert column values to factors
              if (appraiser %in% names(data)) {
                data[[appraiser]] <- factor(data[[appraiser]])
              } else{
                stop(appraiser, "is not a valid column name for", deparse(substitute(data)))
              }
            } else {
              stop("A data.frame object is needed as data argument")
            }

            part <- as.factor(part)
            appraiser <- as.factor(appraiser)

            #Number of parts levels and replicates
            lvlPart = nlevels(data[[part]])
            lvlAppr = nlevels(data[[appraiser]])
            n = nrow(data)/lvlPart

            callNextMethod(.Object, ..., id = id, name = name, description = description, pro = pro, variable = variable, part = part, appraiser = appraiser, data = data, usl = usl, lsl = lsl, sigma = sigma, alphaLim = alphaLim, digits = digits, lvlPart = lvlPart, lvlAppr = lvlAppr, n = n)
          })

#' Anova study for Nested MSA
#' @name anovaMSA
#' @export
setMethod("anovaMSA",
          signature = signature(object = "NestedMSA"),
          function(object){
            ## Complete model (with interaction)
            modelf <- as.formula(paste(object@variable, "~", object@appraiser, "/", object@part))

            model <- aov(modelf, data = object@pro@data)
            modelm <- summary(model)

            rownames(modelm[[1]])[3] <- "REPEATIBILITY"

            #Total row for Df and SumSq
            modelm[[1]] <- rbind(modelm[[1]], c(colSums(modelm[[1]][, 1:2]), rep(NA, 3)))

            rownames(modelm[[1]])[4] <- "TOTAL"

            object@anova <- list(modelm[[1]])

            return(object)
          })


#' Gage rar for Nested MSA
#' @name rar
#' @export
setMethod("rar",
          signature = signature(object = "NestedMSA"),
          function(object){
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

            return(object)
          })
