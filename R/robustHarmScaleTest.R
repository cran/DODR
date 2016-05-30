#' Robust test for difference in scales upon harmonic regressions
#'
#' The Test uses the Fligner-Killeen Test for differences in scales. The test
#' remains valid if the noise distribution contains outliers or is
#' not-Gaussian
#'
#' @inheritParams dodr
#'
#' @return data frame with following columns:
#' \itemize{
#' \item{p.value: }{P-value for difference between the two time series}
#' \item{score: }{score from the underlying
#' \code{\link[npsm]{fk.test}} test}
#' \item{factor: }{Measure for the difference between the two fits}
#' }
#' @seealso \code{\link[npsm]{fk.test}}
#' @export
#' @import Rfit
#' @import npsm
#' @import parallel
robustHarmScaleTest <- function(val1, val2,
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
                    NA,
                    NA)
                return(retlist)
            }

            ehatrb1 <- suppressWarnings(rfit(val1[,ind]~xmat1)$resid)
            ehatrb2 <- suppressWarnings(rfit(val2[,ind]~xmat2)$resid)
            fk <- try(fk.test(ehatrb1,ehatrb2), TRUE)

            if(class(fk)=="try-error"){
                fk <- list(p.value=NA, score=NA, estimate=NA)
            }

            res <- c(fk$p.value,
                    fk$statistic,
                    fk$estimate)

            return(res)
    }, mc.preschedule=TRUE, mc.cleanup=TRUE))

    resultdf <- data.frame(results)
    colnames(resultdf) <- c('p.value', 'score', 'factor')
    return(resultdf)
}
