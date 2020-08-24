#' constructor for CrossedMSA
#'
#' This is the constructor.
#' @name initializeCrossedMSA
#' @export
setMethod("initialize",
          signature = signature(.Object = "CrossedMSA"),
          function(.Object, ..., id, name, description, pro, characteristic, data, tolerance, sigma, alphaLim) {

            headers <- names(data@data)

            part <- as.factor(headers[1])
            appraiser <- as.factor(headers[2])

            data@data[[part]] <- factor(data@data[[part]])
            data@data[[appraiser]] <- factor(data@data[[appraiser]])

            ## Number of parts levels and replicates
            lvlPart = nlevels(data@data[[part]])
            lvlAppr = nlevels(data@data[[appraiser]])
            n = nrow(data@data)/(lvlPart * lvlAppr)

            ## Checks for unbalanced designs
            if (abs(n - round(n))) {
              stop("[CrossedMSA: validation] Not balanced design.")
            }

            ## At least 2 measures by combination must be given
            if (n < 2) {
              stop("[CrossedMSA: validation] Not enough number of replications ")
            }


            .Object <- callNextMethod(.Object, ..., id = id, name = name, description = description, pro = pro, characteristic = characteristic, data = data, tolerance = tolerance, sigma = sigma, alphaLim = alphaLim, lvlPart = lvlPart, lvlAppr = lvlAppr, n = n)

            #If data is ready calculations are made on the initialization method call
            if (missing(pro)) {
              message("[CrossedMSA: validation] .Pro to be analysed is not presented. Anova and R&R methods must be executed manually after getting it." )
            } else if (missing(characteristic)) {
              message("[CrossedMSA: validation] .Characteristic to be analysed is not presented. Anova and R&R methods must be executed manually after getting it." )
            } else if (missing(data)) {
              message("[CrossedMSA: validation] .ProData to be analysed is not presented. Anova and R&R methods must be executed manually after getting it." )
            } else {
              .Object <- anovaMSA(.Object)
              .Object <- rar(.Object)
            }
          })

#' Anova study for Crossed MSA
#' @name anovaMSA
#' @export
setMethod("anovaMSA",
          signature = signature(object = "CrossedMSA"),
          function(object){
            headers <- names(object@data@data)

            part <- headers[1]
            appraiser <- headers[2]
            variable <- headers[3]

            ## Complete model (with interaction)
            modelf <- as.formula(paste(variable, "~", appraiser, "*", part))

            model <- aov(modelf, data = object@data@data)
            modelm <- summary(model)

            #
            # HERE WILL BE ERROR TERM IMPLEMENTATION
            #

            rownames(modelm[[1]])[4] <- "REPEATIBILITY"

            ## Total row for Df and SumSq
            modelm[[1]] <- rbind(modelm[[1]], c(colSums(modelm[[1]][, 1:2]), rep(NA, 3)))

            rownames(modelm[[1]])[5] <- "TOTAL"

            object@anova <- list(modelm[[1]])

            ## F test for interaction
            p <- modelm[[1]][3,5]

            ## REDUCED MODEL
            ## If P is bigger than the alpha limit interaction is removed
            if (p > object@alphaLim) {
              reducedModelf <- as.formula(paste(variable, "~", part, "+", appraiser))

              reducedModel <- aov(reducedModelf, data = object@data@data)
              reducedModelm <- summary(reducedModel)

              rownames(reducedModelm[[1]])[3] <- "REPEATIBILITY"

              ## Total row for Df and SumSq
              reducedModelm[[1]] <- rbind(reducedModelm[[1]], c(colSums(reducedModelm[[1]][, 1:2]), rep(NA, 3)))

              rownames(reducedModelm[[1]])[4] <- "TOTAL"

              object@anovaReduced <- list(reducedModelm[[1]])
            }

            return(object)
          })


#' Gage rar for Crossed MSA
#' @name rar
#' @export
setMethod("rar",
          signature = signature(object = "CrossedMSA"),
          function(object){

            headers <- names(object@data@data)

            part <- headers[1]
            appraiser <- headers[2]
            variable <- headers[3]

            ## if anova is not stored we calcute it
            if(!length(object@anova)) {
              object <- anovaMSA(object)
            }

            object@varianceComponents = matrix(ncol = 6, nrow = 7)

            rownames(object@varianceComponents) <- c("Total Gage R&R", "  Repeatability", "  Reproducibility", paste("Reproducibility by ", appraiser), paste("Reproducibility by ", part, ":", appraiser), "Part-To-Part", "Total Variation")

            colnames(object@varianceComponents) <- c("VarComp", "%Contrib", "StdDev", "StudyVar", "%StudyVar", "%Tolerance")

            ## Variance Components Table
            #Repeatibility
            object@varianceComponents[2, 1] <- object@anova[[1]][4, 3]
            #Appraiser reproducibility
            object@varianceComponents[4, 1] <- max(c((object@anova[[1]][2, 3] - object@anova[[1]][3, 3])/(object@lvlPart * object@n), 0))
            #Part:Appraiser reproducibility
            object@varianceComponents[5, 1] <- max(c((object@anova[[1]][3, 3] - object@anova[[1]][4, 3])/object@n, 0))
            #Total reproducibility
            object@varianceComponents[3, 1] <- object@varianceComponents[4, 1] + object@varianceComponents[5, 1]
            #Part to part
            object@varianceComponents[6, 1] <- max(c((object@anova[[1]][1, 3] - object@anova[[1]][3, 3]) / (object@lvlAppr / object@n), 0))
            #Totat Gage rar
            object@varianceComponents[1, 1] <- object@varianceComponents[2, 1] + object@varianceComponents[3, 1]
            #Total variation
            object@varianceComponents[7, 1] <- object@varianceComponents[1, 1] + object@varianceComponents[6, 1]

            #%Contrib
            object@varianceComponents[, 2] <- round(100 * (object@varianceComponents[, 1]/object@varianceComponents[7, 1]), 2)
            #Standard Deviation
            object@varianceComponents[, 3] <- sqrt(object@varianceComponents[, 1])
            #Study Variation edited from 5.15 to variable
            object@varianceComponents[, 4] <- object@varianceComponents[, 3] * object@sigma
            #%Study Variation
            object@varianceComponents[, 5] <- round(100 * (object@varianceComponents[, 3]/object@varianceComponents[7, 3]), 2)
            #StudyVar/Tolerance
            object@varianceComponents[, 6] <- round(100 * (object@varianceComponents[, 4]/(object@tolerance)), 2)

            #Number of distinct categories
            object@numberCategories <- max(c(1, floor((object@varianceComponents[6, 4]/object@varianceComponents[1, 4])*1.41)))

            ## If reduced model, remove interaction
            p <- object@anova[[1]][3, 5]

            if (p > object@alphaLim){
              object@varianceComponents <- object@varianceComponents[-c(5), ]
            }

            return(object)
          })

#' Plot rar resume
#' @name plotRar
#' @export
setMethod("plotRar",
          signature = signature(object = "CrossedMSA"),
          function(object){

            lay <- rbind(c(1,2),
                         c(3,4),
                         c(5,6))

            g <- grid.arrange(
              ggplotComponentOfVariationChart(object),
              ggplotVariableByPartChart(object),
              ggplotRangeChart(object, gridLayout = TRUE),
              ggplotVariableByAppraiserChart(object),
              ggplotMeanChart(object, gridLayout = TRUE),
              ggplotInteractionChart(object),
              layout_matrix = lay,
              top = object@name,
              bottom = object@description
            )

            g2 <- cowplot::ggdraw(g) +
              theme(plot.background = element_rect(fill="white", color = NA))

            plot(g2)
          })

#' Components of Variation Chart
#' @name plotComponentOfVariationChart
#' @export
setMethod("plotComponentOfVariationChart",
          signature = signature(object = "CrossedMSA"),
          function(object){

            ## Set rows and cols to take from components of variation table to be printed
            p <- object@anova[[1]][3, 5]

            if (p > object@alphaLim){
              rows <- c(1,2,3,5)
            } else {
              rows <- c(1,2,3,6)
            }

            rlabels <- c("G.R&R", "Repeat", "Reprod", "Part2Part")

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

#' Components of Variation Chart with ggplot2
#' @name ggplotComponentOfVariationChart
#' @export
setMethod("ggplotComponentOfVariationChart",
          signature = signature(object = "CrossedMSA"),
          function(object){

            ## Set rows and cols to take from components of variation table to be printed
            p <- object@anova[[1]][3, 5]

            if (p > object@alphaLim){
              rows <- c(1,2,3,5)
            } else {
              rows <- c(1,2,3,6)
            }

            rlabels <- c("G.R&R", "Repeat", "Reprod", "Part2Part")

            if ((!is.na(object@characteristic@U) && !is.na(object@characteristic@U)) || !is.na(object@tolerance)) {
              cols <- c(2, 5, 6)
              clabels <- c("Contribution", "Stud.Var", "Tolerance")
            } else{
              cols <- c(2, 5)
              clabels <- c("Contribution", "Study.Var")
            }

            chartData <- object@varianceComponents[rows,cols]
            rownames(chartData) <- rlabels
            colnames(chartData) <- clabels

            chartData <- cbind(rownames(chartData), data.frame(chartData, row.names=NULL))
            colnames(chartData) <- c("Variance.Component",clabels)

            contributionData <- data.frame(chartData[,c(1,2)], "Contribution")
            studyVarData <- data.frame(chartData[,c(1,3)], "StudyVar")

            colnames(contributionData) <- c("Component", "Value", "Rate")
            colnames(studyVarData) <- c("Component", "Value", "Rate")

            if (length(clabels) == 3) {
              toleranceVarData <- data.frame(chartData[,c(1,3)], "Tolerance")
              colnames(toleranceVarData) <- c("Component", "Value", "Rate")

              data <- rbind.data.frame(contributionData, studyVarData, toleranceVarData)

            } else {
              data <- rbind.data.frame(contributionData, studyVarData)
            }

            gg <- ggplot(data = data, aes_string("Component", "Value", fill = "Rate")) +
              geom_col(position = "dodge") +
              labs(title="Components of Variation by rate")

            gg <- gg +
              geom_hline(yintercept = 10, lty = 2, col = "grey") +
              geom_hline(yintercept = 30, lty = 2, col = "grey")

            print(gg)
          })

#' Variable by Part
#' @name plotVariableByPartChart
#' @export
setMethod("plotVariableByPartChart",
          signature = signature(object = "CrossedMSA"),
          function(object){
            headers <- names(object@data@data)

            part <- headers[1]
            appraiser <- headers[2]
            variable <- headers[3]

            ## Formula for the chart
            f <- as.formula(paste(variable, "~",  part))

            stripchart(f, data = object@data@data, vertical = TRUE,
                               method = "overplot", main = paste(variable, "by", part), pch = 1)

            grid()
          })

#' Variable by Part
#' @name ggplotVariableByPartChart
#' @export
setMethod("ggplotVariableByPartChart",
          signature = signature(object = "CrossedMSA"),
          function(object){
            headers <- names(object@data@data)

            part <- headers[1]
            appraiser <- headers[2]
            variable <- headers[3]

            ## Formula for the chart
            f <- as.formula(paste(variable, "~",  part))

            meanByPart <- aggregate(f, data = object@data@data, mean)

            show(meanByPart)

            show(object@data@data)

            gg <- ggplot(data = object@data@data, aes_string(x=part, y=variable, group = 1, color = part)) +
              geom_point() +
              labs(title=paste(variable, "by", part), x=part, y=variable)

            gg <- gg + geom_line(data = meanByPart, aes_string(x=part, y=variable, group = 1), colour = "grey")

            print(gg)
          })

#' Variable by Appraiser
#' @name plotVariableByAppraiserChart
#' @export
setMethod("plotVariableByAppraiserChart",
          signature = signature(object = "CrossedMSA"),
          function(object){
            headers <- names(object@data@data)

            appraiser <- headers[2]
            variable <- headers[3]

            ## Formula for the chart
            f <- as.formula(paste(variable, "~",  appraiser))

            stripchart(f, data = object@data@data, vertical = TRUE,
                               method = "overplot", main = paste(variable, "by", appraiser), pch = 1)

            grid()
          })

#' Variable by Appraiser with ggplot2
#' @name ggplotVariableByAppraiserChart
#' @export
setMethod("ggplotVariableByAppraiserChart",
          signature = signature(object = "CrossedMSA"),
          function(object){
            headers <- names(object@data@data)

            part <- headers[1]
            appraiser <- headers[2]
            variable <- headers[3]

            ## Formula for the chart
            f <- as.formula(paste(variable, "~",  appraiser))

            meanByPart <- aggregate(f, data = object@data@data, mean)

            show(meanByPart)

            show(object@data@data)

            gg <- ggplot(data = object@data@data, aes_string(x=appraiser, y=variable, group = 1, color = appraiser)) +
              geom_point() +
              labs(title=paste(variable, "by", appraiser), x=appraiser, y=variable)

            gg <- gg + geom_line(data = meanByPart, aes_string(x=appraiser, y=variable, group = 1), colour = "grey")

            print(gg)
          })

#' Interaction
#' @name plotInteractionChart
#' @export
setMethod("plotInteractionChart",
          signature = signature(object = "CrossedMSA"),
          function(object){
            headers <- names(object@data@data)

            part <- headers[1]
            appraiser <- headers[2]
            variable <- headers[3]

            ## Formula for the chart
            f <- as.formula(paste(variable, "~",  appraiser, "+", part))

            agregatedData <- aggregate(f, data = object@data@data, mean)

            minY <- min(agregatedData[[variable]]) - 0.1 * diff(range(agregatedData[[variable]]))
            maxY <- max(agregatedData[[variable]]) + 0.1 * diff(range(agregatedData[[variable]]))

            distinctAppraisers <- unique(agregatedData[[appraiser]])

            for (i in seq_along(distinctAppraisers)) {
              filterIndexs <- agregatedData[,1] == distinctAppraisers[[i]]

              data = agregatedData[filterIndexs,]

              ## To avoid boxplot to be printed instead of xyplot
              dataX <- as.numeric(levels(data[[part]]))[data[[part]]]

              if (i == 1) {
                plot(x = dataX, y = data[[variable]], ylim = c(minY, maxY), xaxt = "n", type = "b", pch = i, lty = i, ylab = paste(variable), xlab = paste(part))
                axis(1, seq_along(data[[part]]), data[[part]])

                title(paste0(part, ":", appraiser, " Interaction"), line = -2, outer = TRUE, font = 2)
              } else {
                points(x = dataX, y = data[[variable]], pch = i)
                lines(x = dataX, y = data[[variable]], lty=i)
              }
            }

            legend(x = 1, y = maxY, legend=distinctAppraisers, pch=seq_along(distinctAppraisers),lty=seq_along(distinctAppraisers), ncol=1)

            grid()
          })

#' Interaction
#' @name ggplotInteractionChart
#' @export
setMethod("ggplotInteractionChart",
          signature = signature(object = "CrossedMSA"),
          function(object){
            headers <- names(object@data@data)

            part <- headers[1]
            appraiser <- headers[2]
            variable <- headers[3]

            ## Formula for the chart
            f <- as.formula(paste(variable, "~",  appraiser, "+", part))

            agregatedData <- aggregate(f, data = object@data@data, mean)

            minY <- min(agregatedData[[variable]]) - 0.1 * diff(range(agregatedData[[variable]]))
            maxY <- max(agregatedData[[variable]]) + 0.1 * diff(range(agregatedData[[variable]]))

            gg <- ggplot(data = agregatedData, aes_string(x=part, y=variable, group = appraiser, color = appraiser)) +
              geom_point() +
              geom_line() +
              labs(title=paste("Interaction Chart by", appraiser), x=part, y=variable)

            gg <- gg + ylim(c(minY, maxY))

            print(gg)
          })

#' Mean Chart
#' @name plotMeanChart
#' @export
setMethod("plotMeanChart",
          signature = signature(object = "CrossedMSA"),
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

#' Mean Chart with ggplot2
#' @name plotGridMeanChart
#' @export
setMethod("ggplotMeanChart",
          signature = signature(object = "CrossedMSA"),
          function(object, gridLayout = FALSE, ...){
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


            ## Plotting
            distinctAppraisers <- unique(xmean[[appraiser]])

            print(xmean)

            gg <- ggplot(data = xmean, aes_string(x=part, y=variable, group = 1, color=appraiser)) +
              geom_point() +
              geom_line() +
              facet_wrap(as.formula(paste("~", appraiser)), ncol=length(distinctAppraisers), drop=TRUE) +
              labs(title=paste("Mean Chart by", appraiser), x=part, y=variable)

            if(gridLayout == TRUE) {
              labelsSize = 2
            } else {
              labelsSize = 5
            }

            gg <- gg +
              geom_hline(yintercept = meanbar, col = 'grey', lty = 2) +
              geom_text(aes(2, meanbar, label = "MEAN"), colour = "grey", vjust = "top", nudge_y = -0.1, size = labelsSize, alpha = 0.8)

            gg <- gg +
              geom_hline(yintercept = ucl, col = 'red') +
              geom_text(aes(2, ucl, label = "UCL"), colour = "red", vjust = "bottom", nudge_y = 0.1, size = labelsSize, alpha = 0.8)

            gg <- gg +
              geom_hline(yintercept = lcl, col = 'red') +
              geom_text(aes(2, lcl, label = "LCL"), colour = "red", vjust = "bottom", nudge_y = 0.1, size = labelsSize, alpha = 0.8)

            gg <- gg + theme(legend.position = "none") + ylim(graphLimits)

            print(gg)

          })

#' Range Chart
#' @name plotRangeChart
#' @export
setMethod("plotRangeChart",
          signature = signature(object = "CrossedMSA"),
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

#' Range Chart with ggplot2
#' @name ggplotRangeChart
#' @export
setMethod("ggplotRangeChart",
          signature = signature(object = "CrossedMSA"),
          function(object, gridLayout = FALSE, ...){
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

            gg <- ggplot(data = xrange, aes_string(x=part, y=variable, group = 1, color=appraiser)) +
              geom_point() +
              geom_line() +
              facet_wrap(as.formula(paste("~", appraiser)), ncol=length(distinctAppraisers), drop=TRUE) +
              labs(title=paste("Range Chart by", appraiser), x=part, y=variable)

            if(gridLayout == TRUE) {
              labelsSize = 2
            } else {
              labelsSize = 5
            }

            gg <- gg +
              geom_hline(yintercept = averageRange, col = 'grey', lty = 2) +
              geom_text(aes(2, averageRange, label = "MEAN"), colour = "grey", vjust = "top", nudge_y = -0.1, size = labelsSize, alpha = 0.8)

            gg <- gg +
              geom_hline(yintercept = url, col = 'red') +
              geom_text(aes(2, url, label = "URL"), colour = "red", vjust = "bottom", nudge_y = 0.1, size = labelsSize, alpha = 0.8)

            gg <- gg +
              geom_hline(yintercept = lrl, col = 'red') +
              geom_text(aes(2, lrl, label = "LRL"), colour = "red", vjust = "bottom", nudge_y = 0.1, size = labelsSize, alpha = 0.8)

            gg <- gg + theme(legend.position = "none") + ylim(graphLimits)

            print(gg)

          })
