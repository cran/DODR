#' Creation of Matrices used in linear regression
#'
#' @param times : timepoints of measurement sequence
#' @param period: exspected period of oscillation
#'
#' @return list containing two matrices (X and Cx)
#' @noRd
matrixPreparation <- function(times, period){
    timesRad <- times / period * 2 * pi
    xMat <- cbind(1, cos(timesRad), sin(timesRad))
    cxMat <- solve(t(xMat) %*% xMat)
    rList <- list(X = xMat, CX = cxMat)
    return(rList)
}

#' Creation of Matrices used in linear regression
#'
#' @param times : timepoints of measurement sequence
#' @param period: exspected period of oscillation
#'
#' @return list containing two matrices (X and Cx)
#' @noRd
matrixPreparationReduced <- function(times, period){
    timesRad <- times / period * 2 * pi
    xMat <- cbind(cos(timesRad), sin(timesRad))
    cxMat <- solve(t(xMat) %*% xMat)
    rList <- list(X = xMat, CX = cxMat)
    return(rList)
}



#' Linear fitting without NA
#'
#' Small unflexible script to apply standard linear fitting
#' if no NA's are in the val matrix
#' @param mat matrices used for calculations
#' @param val matrix of values, one sample per col one time
#' point per row
#'
#' @return result of linear fit
#' b: fitted parameters
#' ssr : sum of squared residuals
#' df: degrees of freedom
#' @noRd
linFit <- function (mat, val) {
    bVec <- mat$CX %*% t(mat$X) %*% val
    ssr <- apply((val - mat$X %*% bVec)^2, 2, sum)
    df <- nrow(mat$X)-ncol(mat$X)
    reList <- list(b = bVec, ssr = ssr, df = df)
}


#' Linear fitting with NAs
#'
#' Small unflexible script to apply standard linear fitting
#' NA's are in the val matrix
#' @param mat matrices used for calculations
#' @param val matrix of values, one sample per col one time
#' point per row
#'
#' @return result of linear fit
#' b: fitted parameters
#' ssr : sum of squared residuals
#' df: degrees of freedom
#' @noRd
#' @import parallel
linFitSep <- function (mat, val) {
    resmat <- lapply(seq_len(ncol(val)), function(ind){
        notNA <- !is.na(val[,ind])
        std <- all(notNA)
        if(std){
            bVec <- mat$CX %*% t(mat$X) %*% val[,ind]
            ssr <- apply((val[,ind] - mat$X %*% bVec)^2, 2, sum)
            df <- nrow(mat$X)-ncol(mat$X)
        }
        else{
            X <- mat$X[notNA,]
            CX <- solve(t(X) %*% X)
            vals <- matrix(val[notNA, ind], ncol=1)
            bVec <- CX %*% t(X) %*% vals
            ssr <- apply((vals - X %*% bVec)^2, 2, sum)
            df <- nrow(X)-ncol(X)
        }
        singleret <- list(bVec, ssr, df)
    })
    do.call(cbind, lapply(resmat, '[[', 1))

    reList <- list(b = do.call(cbind, lapply(resmat, '[[', 1)),
                ssr = do.call(c, lapply(resmat, '[[', 2)),
                df = do.call(c, lapply(resmat, '[[', 3))
                )
    return(reList)
}

#' Linear fitting
#'
#' Tests whether NA's are in value matrix and starts acording script
#' @param mat matrices used for calculations
#' @param val matrix of values, one sample per col one time
#' point per row
#'
#' @return result of linear fit
#' b: fitted parameters
#' ssr : sum of squared residuals
#' df: degrees of freedom
#' @noRd
linFitt <- function (mat, val) {
    if(all(!is.na(val))){
        return(linFit(mat, val))
    }
    else{
        return(linFitSep(mat, val))
    }
}
