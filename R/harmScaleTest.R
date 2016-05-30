#' Test for difference in scales upon harmonic regressions
#'
#' The Test uses the F-Test for variances to decide wether the two given
#' time series have a comparable noise scale or on of the series has a higher
#' noise level.
#'
#' @inheritParams dodr
#'
#' @return dataframe with following columns:
#' \itemize{
#' \item{p.value: }{P-value for difference between the two time series}
#' \item{score: }{score from the underlying F test}
#' }
#' @export
#'
#' @importFrom stats lm pf var
harmScaleTest <- function(val1, val2,
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

    t1 <- times1 /period*2*pi
    t2 <- times2 /period*2*pi
    n <- length(times2)
    m <- length(times1)
    p <- 3

    xmat1 <- cbind(1, cos(t1), sin(t1))
    xmat2 <- cbind(1, cos(t2), sin(t2))

    tp1 <- length(times1)
    tp2 <- length(times2)

    muchna <- apply(val1, 2, function(x) sum(is.na(x))-tp1 >= -3) |
        apply(val2, 2, function(x)  sum(is.na(x))-tp1 >= -3)

    results <- do.call(rbind,
        mclapply(seq_len(ncol(val1)), function(ind){

            if(muchna[ind]){
                retlist <- c(
                    NA,
                    NA)
                return(retlist)
            }

            ehatls1 <- lm(val1[,ind]~xmat1)$resid
            ehatls2 <- lm(val2[,ind]~xmat2)$resid
            ftst <- max(c(var(ehatls2),var(ehatls1)))/min(c(var(ehatls2),var(ehatls1)))
            n <- length(ehatls2)
            m <- length(ehatls1)
            pval <- 2*(1-pf(ftst,n-p,m-p))
            res <- c(pval,
                    ftst
                    )
            return(res)
        }, mc.preschedule=TRUE, mc.cleanup=TRUE))

    resultdf <- data.frame(results)
    colnames(resultdf) <- c('p.value', 'F')
    return(resultdf)
}
