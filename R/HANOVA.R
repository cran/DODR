#' ANOVA applied on harmonic regression
#'
#' Detection of differential between rhythms in two time series using llsq
#' fits and ANOVA
#'
#' @inheritParams dodr
#'
#' @details This test uses general ANOVA to test for differences between two
#' time series. Therefore the time series are fitted to sine curves with a
#' fixed period length and free phase and amplitude. In one case phase and
#' amplitude have two be the same for both series, in the other case phase
#' and amplitude could differ for the two series.
#' @return data frame with columns:
#' \itemize{
#' \item{p.value: }{P-value for difference between the two time series}
#' \item{F: }{F score from the underlying ANOVA test}
#' \item{diff: }{Measure for the difference between the two fits}
#' }
#' @export
#' @importFrom methods as
#' @importFrom stats anova coefficients lm na.omit
HANOVA <- function (val1, val2,
                            times1, times2,
                            period, norm = TRUE,
                            verbose=options('verbose')[[1]]
                        ){

    if(verbose) message("Preparing HANOVA... ")
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


    #cosE <- cosB*c(rep(1,length(times1)),rep(0,length(times2)))
    #sinE <- sinB*c(rep(1,length(times1)),rep(0,length(times2)))
    #cosF <- cosB*c(rep(0,length(times1)),rep(1,length(times2)))
    #sinF <- sinB*c(rep(0,length(times1)),rep(1,length(times2)))
    #combine both value matrices
    valC <- rbind(val1, val2)

    tp1 <- length(times1)
    tp2 <- length(times2)

    muchna <- apply(val1, 2, function(x) sum(is.na(x))-tp1 >= -3) |
        apply(val2, 2, function(x)  sum(is.na(x))-tp1 >= -3)

    xred <- cbind(cosB, sinB)
    xfull <- cbind(cosB, sinB, cosD, sinD)

    relist <- function(obj){
        len <- ncol(obj$coefficients)
        lm.list <- lapply(seq_len(len), function(index){
            ret <- list(coefficients = obj$coefficients[,index],
                        residuals = obj$residuals[,index],
                        effects = obj$effects[,index],
                        df.residual = obj$df.residual,
                        model = obj$model,
                        call = obj$call)
            class(ret) = 'lm'
            ret <- as(ret, 'lm')
            return(ret)
        })

    }

    if(verbose) message("Fitting HANOVA... ")

    #apply fitting

    has.na <- TRUE#!all(!is.na(valC))

#     if(has.na){
#         f.rl <- mclapply(seq_len(ncol(valC)), function(ind){
#             lm(valC[,ind] ~ xred, na.action = na.omit)
#         })
#         f.fl <- mclapply(seq_len(ncol(valC)), function(ind){
#             lm(valC[,ind] ~ xfull, na.action = na.omit)
#         })
#     } else{
#
#         f.r <- lm(valC[,ser] ~ xred, na.action = na.omit)
#         f.f <- lm(valC ~ xfull, na.action = na.omit)
#
#         f.rl <- relist(f.r)
#         f.fl <- relist(f.f)
#
#     }

    if(verbose) message("ANOVA calculation HANOVA... ")

    results <- do.call(rbind,
       mclapply(seq_len(ncol(val1)), function(index){
        vals <- valC[, index]
        if(muchna[index]){
            retlist <- c(
                NA,
                NA,
                NA)
            return(retlist)
        }
        f.r <- lm(vals ~ xred, na.action = na.omit)
        f.f <- lm(vals ~ xfull, na.action = na.omit)

        dt1 <- anova(f.r, f.f)

        if("try-error" %in% class(dt1)){
            dt1 <- list('Pr(>F)'=1, F=0)
            diff <- 0
        } else{
            diff <- sqrt(sum(coefficients(f.f)[4:5]^2))
        }
        retlist <- c(
            dt1$'Pr(>F)'[2],
            dt1$F[2],
            diff)
        return(retlist)
    }, mc.preschedule=TRUE, mc.cleanup=TRUE))

    resultdf <- data.frame(results)
    colnames(resultdf) <- c('p.value', 'F', 'diff')
    return(resultdf)

}
