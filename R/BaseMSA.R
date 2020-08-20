#' constructor for BaseMSA
#'
#' This is the constructor.
#' @name initializeBaseMSA
#' @export
setMethod("initialize",
          signature = signature(.Object = "BaseMSA"),
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

            #If data is ready calculations are made on the initialization method call
            if (missing(pro)) {
              message("[BaseMSA: validation] .Pro to be analysed is not presented. Anova and R&R methods must be executed manually after getting it." )
            } else if (missing(characteristic)) {
              message("[BaseMSA: validation] .Characteristic to be analysed is not presented. Anova and R&R methods must be executed manually after getting it." )
            } else if (missing(data)) {
              message("[BaseMSA: validation] .ProData to be analysed is not presented. Anova and R&R methods must be executed manually after getting it." )
            } else {
              .Object <- anovaMSA(.Object)
              .Object <- rar(.Object)
            }

          })

#' Anova study for Base MSA
#' @name anovaMSA
#' @export
setMethod("anovaMSA",
          signature = signature(object = "BaseMSA"),
          function(object){
            headers <- names(object@data@data)

            part <- headers[1]
            appraiser <- headers[2]
            variable <- headers[3]

            ## Complete model (with interaction)
            modelf <- as.formula(paste(variable, "~", part))

            model <- aov(modelf, data = object@data@data)
            modelm <- summary(model)

            rownames(modelm[[1]])[2] <- "REPEATIBILITY"

            #Total row for Df and SumSq
            modelm[[1]] <- rbind(modelm[[1]], c(colSums(modelm[[1]][, 1:2]), rep(NA, 3)))

            rownames(modelm[[1]])[3] <- "TOTAL"

            object@anova <- list(modelm[[1]])

            return(object)
          })


#' Gage rar for Base MSA
#' @name rar
#' @export
setMethod("rar",
          signature = signature(object = "BaseMSA"),
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
            object@varianceComponents[3, 2] <- NA
            #Part to part
            object@varianceComponents[4, 1] <- max(c((object@anova[[1]][2, 3] - object@anova[[1]][3, 3]) / object@n, 0))
            #Totat Gage rar
            object@varianceComponents[1, 1] <- object@varianceComponents[2, 1]
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

            object@varianceComponents <- object@varianceComponents[-c(3:5), ]

            return(object)
          })

#' Components of Variation Chart
#' @name plotComponentOfVariationChart
#' @export
setMethod("plotComponentOfVariationChart",
          signature = signature(object = "BaseMSA"),
          function(object){

            ## Set rows and cols to take from components of variation table to be printed
            rows <- c(1,2,3)
            rlabels <- c("G.R&R", "Repeat", "Part2Part")

            if ((!is.na(object@characteristic@U) && !is.na(object@characteristic@U)) || !is.na(object@tolerance)) {
              cols <- c(2, 5, 6)
              clabels <- c("%Contribution", "%Study Var", "%Tolerance")
            } else{
              cols <- c(2, 5)
              clabels <- c("%Contribution", "%Study Var")
            }

            chartData <- object@varianceComponents[rows,cols]
            rownames(chartData) <- rlabels
            colnames(chartData) <- clabels

            plot <- barplot(height = t(chartData), beside = TRUE,
                            main = "Components of Variation",
                            legend.text = clabels, ylim = c(0,100), axes = TRUE)

            grid()

            return(plot)
          })

#' Variable by Part
#' @name plotVariableByPartChart
#' @export
setMethod("plotVariableByPartChart",
          signature = signature(object = "BaseMSA"),
          function(object){
            headers <- names(object@data@data)

            part <- headers[1]
            variable <- headers[3]

            ## Formula for the chart
            f <- as.formula(paste(variable, "~",  part))

            plot <- stripchart(f, data = object@data@data, vertical = TRUE,
                               method = "jitter", main = paste(variable, "by", part),
                               xlab = part)

            grid()
          })

#' Variable by Appraiser
#' @name plotVariableByAppraiserChart
#' @export
setMethod("plotVariableByAppraiserChart",
          signature = signature(object = "BaseMSA"),
          function(object){
            headers <- names(object@data@data)

            appraiser <- headers[2]
            variable <- headers[3]

            ## Formula for the chart
            f <- as.formula(paste(variable, "~",  appraiser))

            plot <- stripchart(f, data = object@data@data, vertical = TRUE,
                               method = "jitter", main = paste(variable, "by", appraiser))

            grid()
          })

#' Mean Chart
#' @name plotMeanChart
#' @export
setMethod("plotMeanChart",
          signature = signature(object = "BaseMSA"),
          function(object){
            headers <- names(object@data@data)

            part <- headers[1]
            appraiser <- headers[2]
            variable <- headers[3]

            ## Mean and range formula
            f <- as.formula(paste(variable, "~", appraiser, "+", part))

            rangeFunction <- function(x) {
              max(x) - min(x)
            }

            xmean <- aggregate(f, data = object@data@data, mean)
            xrange <- aggregate(f, data = object@data@data, rangeFunction)

            averageRange <- mean(xrange[[variable]])

            meanbar <- mean(object@data@data[[variable]], na.rm = TRUE)

            ucl <- meanbar + (3/(ss.cc.getd2(object@n)*sqrt(object@n)))*averageRange
            lcl <- meanbar - (3/(ss.cc.getd2(object@n)*sqrt(object@n)))*averageRange

            graphLimits <- c(min(range(xmean[[variable]])[1], lcl),
                             max(range(xmean[[variable]])[2], ucl)) +
              c(-1, 1) * 0.1* diff(range(xmean[[variable]]))

            ## Formula for chart
            #chartf <- as.formula(paste(variable, "~", part))

            ## Plotting
            distinctAppraisers <- unique(xmean[[appraiser]])
            minX <- min(xmean[[variable]])
            maxX <- max(xmean[[variable]])

            ## Save te previous layout to restore it after printing the plot
            par_temp = par()
            par(mfrow = c(1, length(distinctAppraisers)), mar=c(5,4,7,2)+0.1, xpd=FALSE)

            for (i in seq_along(distinctAppraisers)) {
              filterIndexs <- xmean[,1] == distinctAppraisers[[i]]

              data = xmean[filterIndexs,]

              ## To avoid boxplot to be printed instead of xyplot
              data[[part]] <- as.numeric(levels(data[[part]]))[data[[part]]]

              plot(x = data[[part]], y = data[[variable]], ylim = graphLimits, type = "b", pch = 1, ylab = paste(variable), xlab = paste(part), col = "blue")

              title(distinctAppraisers[[i]], line = 1)

              grid()

              abline(h = meanbar, col = 'grey', lty = 2)      # mean
              abline(h = ucl, col = "red")
              abline(h = lcl, col = "red")

              text(y = meanbar, x = 1.35, expression(bold(bar(X))), cex=1, pos=3, col="grey")

              text(y = ucl, x = 1.35, "UCL", cex=1, pos=3, col="red")
              text(y = lcl, x = 1.35, "LCL", cex=1, pos=1, col="red")
            }

            mtext(paste("Mean Chart by", appraiser), side = 3, line = -4, outer = TRUE, font = 2)

            ## Restore previous layout
            par(par_temp)
          })

#' Range Chart
#' @name plotRangeChart
#' @export
setMethod("plotRangeChart",
          signature = signature(object = "BaseMSA"),
          function(object){
            headers <- names(object@data@data)

            part <- headers[1]
            appraiser <- headers[2]
            variable <- headers[3]

            ## Mean and range formula
            f <- as.formula(paste(variable, "~", appraiser, "+", part))

            rangeFunction <- function(x) {
              max(x) - min(x)
            }

            xrange <- aggregate(f, data = object@data@data, rangeFunction)

            averageRange <- mean(xrange[[variable]])

            d3 <- ss.cc.getd3(object@n)
            d2 <- ss.cc.getd2(object@n)

            ## Range limits
            url <- averageRange*(1 + 3 * (d3/d2))
            lrl <- max(averageRange * (1 - 3 * (d3/d2)), 0)

            ## Graph limits
            graphLimits <- c(min(range(xrange[[variable]])[1], lrl),
                             max(range(xrange[[variable]])[2], url)) +
              c(-1, 1) * 0.1 * diff(range(xrange[[variable]]))

            ## Ploting
            distinctAppraisers <- unique(xrange[[appraiser]])
            minX <- min(xrange[[variable]])
            maxX <- max(xrange[[variable]])

            ## Save te previous layout to restore it after printing the plot
            par_temp = par()
            par(mfrow = c(1, length(distinctAppraisers)), mar=c(5,4,7,2)+0.1, xpd=FALSE)

            for (i in seq_along(distinctAppraisers)) {
              filterIndexs <- xrange[,1] == distinctAppraisers[[i]]

              data = xrange[filterIndexs,]

              ## To avoid boxplot to be printed instead of xyplot
              data[[part]] <- as.numeric(levels(data[[part]]))[data[[part]]]

              plot(x = data[[part]], y = data[[variable]], ylim = graphLimits, type = "b", pch = 1, ylab = paste(variable), xlab = paste(part), col = "blue")

              title(distinctAppraisers[[i]], line = 1)

              grid()

              abline(h = averageRange, col = 'grey', lty = 2)    # average Range
              abline(h = url, col = "red")    # upper Range limit
              abline(h = lrl, col = "red")    # lower Range limit

              text(y = averageRange, x = 1.35, expression(bold(bar(R))), cex=1, pos=3, col="grey")

              text(y = url, x = 1.35, "URL", cex=1, pos=3, col="red")
              text(y = lrl, x = 1.35, "LRL", cex=1, pos=1, col="red")
            }

            mtext(paste("Range Chart by", appraiser), side = 3, line = -4, outer = TRUE, font = 2)

            ## Restore previous layout
            par(par_temp)
          })
