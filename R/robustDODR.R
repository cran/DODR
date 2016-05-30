#' Asymmetric variant of Chow method applied on harmonic regression
#'
#' @inheritParams dodr
#'
#' @details This test uses a robust Fitting drop test to test for differences
#' between two time series.
#' Therefore the time series are fitted to sine curves with a
#' fixed period length and free phase and amplitude. In one case phase and
#' amplitude have two be the same for both series, in the other case phase
#' and amplitude could differ for the two series.
#' @return data frame with following columns:
#' \itemize{
#' \item{p.value: }{P-value for difference between the two time series}
#' \item{F: }{F score from the underlying \code{\link[Rfit]{drop.test}} test}
#' \item{diff: }{Measure for the difference between the two fits}
#' }
#' @export
#' @import Rfit parallel
#' @importFrom stats coefficients
robustDODR <- function(val1, val2,
                    times1, times2,
                    period, norm = TRUE){


    #prepare val1 and val2 if they contain only one time series
    if(length(dim(val1)) == 1){
        val1 <- matrix(val1, ncol=1)
    }
    if(length(dim(val2)) == 1){
        val2 <- matrix(val2, ncol=1)
    }

    #check for input integrity

    errMsg <- checkIntegrity(val1, val2, times1, times2)
    if(! is.null(errMsg)) stop(errMsg)

    #normalisation
    if(norm){
        val1 <- apply(val1, 2, function(series){
            normval <- series/mean(series, na.rm=TRUE)
            return(normval)
        })

        val2 <- apply(val2, 2, function(series){
            normval <- series/mean(series, na.rm=TRUE)
            return(normval)
        })
    }

    #prepare matrices for Rfit calculations
    timevec <- c(times1, times2)/period*2*pi

    #vectors for combined fit
    cosB <- cos(timevec)
    sinB <- sin(timevec)

    #vetors for individual Fits
    cosD <- cosB*c(rep(0,length(times1)),rep(1,length(times2)))
    sinD <- sinB*c(rep(0,length(times1)),rep(1,length(times2)))


    cosE <- cosB*c(rep(1,length(times1)),rep(0,length(times2)))
    sinE <- sinB*c(rep(1,length(times1)),rep(0,length(times2)))
    cosF <- cosB*c(rep(0,length(times1)),rep(1,length(times2)))
    sinF <- sinB*c(rep(0,length(times1)),rep(1,length(times2)))
    #combine both value matrices
    valC <- rbind(val1, val2)

    tp1 <- length(times1)
    tp2 <- length(times2)

    muchna <- apply(val1, 2, function(x) sum(is.na(x))-tp1 >= -3) |
        apply(val2, 2, function(x)  sum(is.na(x))-tp1 >= -3)

    #apply fitting
    results <- do.call(rbind,
        mclapply(seq_len(ncol(val1)), function(index){
        if(muchna[index]){
            retlist <- c(
                NA,
                NA,
                NA)
            return(retlist)
        }

        f1.r <- suppressWarnings( rfit(valC[,index] ~ cosB + sinB ))
        f1.f <- suppressWarnings(rfit(valC[,index] ~ cosB + sinB + cosD + sinD
                        ,yhat0=f1.r$fitted))


        dt1 <- try(drop.test(f1.f, f1.r),TRUE)
        if(class(dt1)=="try-error"){
            dt1 <- list(p.value=1, F=0)
            diff <- 0
        }
        else{
            diff <- sqrt(sum(coefficients(f1.f)[4:5]^2))
        }
        retlist <- c(
            dt1$p.value,
            dt1$F,
            diff)



        return(retlist)
    }, mc.preschedule=TRUE, mc.cleanup=TRUE))

    resultdf <- data.frame(results)
    colnames(resultdf) <- c('p.value', 'F', 'diff')
    return(resultdf)

}
