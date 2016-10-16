splineAnalyze2<-function (Y, map, smoothness = 100, s2 = NA, mean = NA, plotRaw = FALSE, 
    plotWindows = FALSE, method = 3) 
{


    # finding x intercept for identification of inflection point
    roots <- function(data) {
        data1 <- c(NA, data[1:{length(data) - 1}])
        data2 <- data
        posneg <- which(data1 > 0 & data2 < 0) - 0.5
        negpos <- which(data1 < 0 & data2 > 0) - 0.5
        zero <- which(data == 0)
        roots <- sort(c(posneg, negpos, zero))
        return(roots)
    }
    
    
    #data input and cleaning
    rawData <- data.frame(Pos = map, Y = Y)
    data <- rawData[which(is.na(rawData$Y) == FALSE), ]
    
    #
    pspline <- smooth.Pspline(data[, 1], data[, 2], norder = 2, method = method)
    predict <- predict(pspline, seq(0, max(pspline$x), by = smoothness))
    psplinederiv <- predict(pspline, seq(0, max(pspline$x), by = smoothness), nderiv = 2)
    psplineInflection <- roots(psplinederiv) * smoothness
    
    print(paste("Total number of windows = ", length(psplineInflection) + 1))
    print(" ---- Computing window statistics ----")
    
    if (is.na(s2)){ 
        s2 <- var(Y)
        }
    if (is.na(mean)){
        mean <- mean(Y, na.rm = TRUE)
        }
        
    cat(1, "of", length(psplineInflection) + 1, "\r")
    
    #generating of windows
    Distinct <- data.frame(WindowStart = rep(NA, length(psplineInflection) + 1), WindowStop = NA, SNPcount = NA, MeanY = NA, Wstat = NA)
    Distinct$WindowStart[1] <- min(pspline$x)
    Distinct$WindowStop[1] <- psplineInflection[1]
    Distinct$SNPcount[1] <- length(which(data[, 1] <= psplineInflection[1]))
    Distinct$MeanY[1] <- mean(data[which(data[, 1] <= psplineInflection[1]), 2], na.rm = TRUE)
    Distinct$Wstat[1] <- {mean(data[which(data[, 1] <= psplineInflection[1]), 2], na.rm = TRUE) - mean}/sqrt(s2/length(which(data[, 1] <= psplineInflection[1])))
    
   
    #,.combine="rbind",.multicombine=T
     Distinct=foreach(i = 2:length(psplineInflection)) %dopar% {
        Distinct$WindowStart[i] <- psplineInflection[i - 1]
        Distinct$WindowStop[i] <- psplineInflection[i]
        Distinct$SNPcount[i] <- length(which(data[, 1] >= psplineInflection[i - 1] & data[, 1] <= psplineInflection[i]))
        Distinct$MeanY[i] <- mean(data[which(data[, 1] >= psplineInflection[i - 1] & data[, 1] <= psplineInflection[i]), 2], na.rm = TRUE)
        Distinct$Wstat[i] <- {mean(data[which(data[, 1] >= psplineInflection[i - 1] & data[, 1] <= psplineInflection[i]), 2], na.rm = TRUE) - mean}/sqrt(s2/length(which(data[, 1] >= psplineInflection[i - 1] & data[, 1] <= psplineInflection[i])))
        return(Distinct)
    }

    return(Distinct)
    
    # i <- i + 1
    # cat(i, "of", length(psplineInflection) + 1, "\r")
    # Distinct$WindowStart[i] <- psplineInflection[i - 1]
    # Distinct$WindowStop[i] <- max(data$Pos)
    # Distinct$SNPcount[i] <- length(which(data[, 1] >= psplineInflection[i - 
    #     1]))
    # Distinct$MeanY[i] <- mean(data[which(data[, 1] >= psplineInflection[i - 
    #     1]), 2], na.rm = TRUE)
    # Distinct$Wstat[i] <- {
    #     mean(data[which(data[, 1] >= psplineInflection[i - 1]), 
    #         2], na.rm = TRUE) - mean
    # }/sqrt(s2/length(which(data[, 1] >= psplineInflection[i - 
    #     1])))
    # print(" ---- done ---- ")
    # if (plotRaw == TRUE & plotWindows == TRUE) 
    #     par(mfrow = c(2, 1))
    # if (plotRaw == TRUE) {
    #     plot(data, xlab = "Position (bp)", ylab = "Raw values")
    #     lines(seq(0, max(pspline$x), by = smoothness), predict, 
    #         col = "red")
    # }
    # if (plotWindows == TRUE) {
    #     plot((Distinct$WindowStop - Distinct$WindowStart)/2 + 
    #         Distinct$WindowStart, Distinct$Wstat, xlab = "Position (bp)", 
    #         ylab = "Spline Wstat", pch = 19)
    # }
    # return(list(rawSpline = pspline, breaks = psplineInflection, 
    #     windowData = Distinct))
  }
