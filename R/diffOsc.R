#' Detection of differences in rhythmic behavior between
#' two time series sets
#'
#' @param val1 matrix: data for values in first set. One column per
#' sample, one row per time point
#' @param val2 matrix: data for values in second set. One column per
#' sample, one row per time point
#' @param times1 vector: times of first set.
#' @param times2 vector: times of second set.
#' @param method vector<string>: method(s) to detect differences.
#' Groups of related methods have additional identifiers.
#' \cr Elementary methods:
#' \code{\link{HANOVA}},
#' \code{\link{harmScaleTest}},
#' \code{\link{harmNoisePred1}},
#' \code{\link{harmNoisePred2}},
#' \code{\link{robustDODR}},
#' \code{\link{robustHarmScaleTest}})
#' \cr Groups:
#' \itemize{
#' \item{"all"}{All methods}
#' \item{"robust"}{Combination of robust detection methods
#' (robustDODR, robustHarmScaleTest)}
#' \item{"lsq"}{Combination of least square based detetion methods
#' (HANOVA, HarmScaleTest, HarmNoisePred1, HarmNoisePred2)}
#' \item{"ANOVA"}{Combination of ANOVA like methods (HANOVA, robustDODR)}
#' \item{"scaleTest"}{Combination of scaleTest methods
#' (harmScaleTest, robustHarmScaleTest)}
#' \item{"harmNoisePred"}{combination of both scaleTest methods to generate a
#' two sided test}
#' }
#' @param verbose boolean: verbosity.
#' @param period numeric: period of the oscillations. Same unit as the
#' time points in times1 and times2
#' @param norm boolean: whether to normalize the time series (division by mean),
#' prior to the analysis.
#'
#' @details This method applies a set of different methods on a pair of two
#' experiments with one measurement matrix each. Samples to compare have to
#' have the same column in both matrices. Different methods could be selected.
#'
#' @return A list containing
#' \itemize{
#' \item{'p.value.table'}{
#' A table containing the p-values for all the tests specified by
#' \code{method}. Each row contains the results for one sample. A column
#' \code{meta.p.val} is added containing the lowest p-value, corrected for
#' multiple testing using a beta-distribution based aproach.}
#' \item{details}{
#' A list containing the detailed results from the selected methods}
#' }
#' @export
#' @examples
#' library(DODR)
#'
#' #defining the parameters for two sets of oscillations
#' n=50
#' testTimes1 <- 0:15*3
#' testTimes2 <- testTimes1
#' tp <- length(testTimes1)
#' per1 <- 24
#' amp1 <- 0.3
#' ph1 <- 5
#' sd1 <- 0.1
#'
#' per2 <- per1
#' amp2 <- amp1
#' ph2 <- ph1+4
#' sd2 <- sd1
#'
#' #creating artificial oscillation sets
#' v1 <- 1 + amp1 * cos((testTimes1 - ph1)/per1*2*pi)
#' noise1 <- rnorm(length(testTimes1)*n, 0, sd1)
#' val1 <- matrix(v1 + noise1, ncol=n)
#'
#' v2 <- 1 + amp2 * cos((testTimes2 - ph2)/per2*2*pi)
#' noise2 <- rnorm(length(testTimes2)*n, 0, sd2)
#' val2 <- matrix(v2 + noise2, ncol=n)
#'
#' # run DODR
#' dodr <- dodr(val1, val2, testTimes1, testTimes2, 24, method = 'all')
#' dodr$p.value.table[1:3,]
#'
#' #create another set with alterations in noise scale
#' ph2 <- ph1
#' sd2 <- sd1 * 3
#'
#' v2 <- 1 + amp2 * cos((testTimes2 - ph2)/per2*2*pi)
#' noise2 <- rnorm(length(testTimes2)*n, 0, sd2)
#' val2 <- matrix(v2 + noise2, ncol=n)
#'
#' dodr <- dodr(val1, val2, testTimes1, testTimes2, 24, method = 'all')
#' dodr$p.value.table[1:3,]
#'
#' @importFrom stats pbeta

dodr <- function(val1, val2,
                                    times1, times2 = times1,
                                    norm = TRUE,
                                    period = 24, method = 'robust',
                                    verbose=options('verbose')[[1]]){

    if(verbose) message("Preparing ... ")
    allowedMethods <- c('HANOVA', 'harmScaleTest',
                        'harmNoisePred1', 'harmNoisePred2',
                        'robustDODR', 'robustHarmScaleTest')

    metaMethods <- c('lsq', 'robust', 'all', 'ANOVA',
                    'scaleTest', 'harmNoisePred')
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

    #translate meta method specifications
    if('all' %in% method){
        method <- c(method, allowedMethods)
    }

    if('robust' %in% method){
        method <- c(method, 'robustDODR', 'robustHarmScaleTest')
    }

    if('lsq' %in% method){
        method <- c(method, 'HANOVA', 'harmScaleTest',
                        'harmNoisePred1', 'harmNoisePred2')
    }

    if('ANOVA' %in% method){
        method <- c(method, 'HANOVA', 'robustDODR')
    }

    if('scaleTest' %in% method){
        method <- c(method, 'harmScaleTest', 'robustHarmScaleTest')
    }

    if('HarmNoisePred' %in% method){
        method <- c(method, 'harmNoisePred1', 'harmNoisePred2')
    }

    if(length(intersect(method, allowedMethods)) <= 0){
        stop("no valid method selected")
    }

    #run the different methods and combine the output

    results <- data.frame(index <- seq_len(ncol(val1)))
    details <- list()

    if(verbose) message("done\n")

    if('HANOVA' %in% method){

        if(verbose) message("HANOVA calculation ... ")

        chow <- HANOVA(val1, val2, times1, times2,
                       period, norm = TRUE, verbose=verbose)
        results <- data.frame(results, HANOVA=chow$p.value)
        details$HANOVA <- chow

        if(verbose) message("done\n")

    }

    if('harmNoisePred1' %in% method){

        if(verbose) message("harmNoisePred1 calculation ... ")

        chowNoise <- harmNoisePred(val1, val2, times1, times2,
                        period, norm = TRUE)
        results <- data.frame(results, HarmNoisePred1=chowNoise$p.value)
        details$HarmNoisePred1 <- chowNoise

        if(verbose) message("done\n")

    }

    if('harmNoisePred2' %in% method){

        if(verbose) message("harmNoisePred2 calculation ... ")

        chowNoise <- harmNoisePred(val2, val1, times2, times1,
                        period, norm = TRUE)
        results <- data.frame(results, HarmNoisePred2=chowNoise$p.value)
        details$HarmNoisePred2 <- chowNoise

        if(verbose) message("done\n")

    }

    if('harmScaleTest' %in% method){

        if(verbose) message("harmScaleTest calculation ... ")

        oara <- harmScaleTest(val1, val2, times1, times2,
                        period, norm = TRUE)
        results <- data.frame(results, HarmScaleTest=oara$p.value)
        details$HarmScaleTest <- oara

        if(verbose) message("done\n")

    }

    if('robustDODR' %in% method){

        if(verbose) message("robustDODR calculation ... ")

        chowSP <- robustDODR(val1, val2, times1, times2,
                        period, norm = TRUE)
        results <- data.frame(results, robustDODR=chowSP$p.value)
        details$robustDODR <- chowSP

        if(verbose) message("done\n")

    }

    if('robustHarmScaleTest' %in% method){

        if(verbose) message("robustHarmScaleTest calculation ... ")

        oara <- robustHarmScaleTest(val1, val2, times1, times2,
                        period, norm = TRUE)
        results <- data.frame(results, robustHarmScaleTest=oara$p.value)
        details$robustHarmScaleTest <- oara

        if(verbose) message("done\n")

    }

    if(verbose) message("Meta P-value calculation ... ")

    resmat <- as.matrix(results[,-1])
    names <- c(colnames(results[,-1]), 'meta.p.val')

    # calculate the corrected minimal P-Value for the multiple testing
    metapval <- sapply(seq_len(nrow(resmat)), function(ind){
        rowraw <- resmat[ind,]
        row <- as.numeric(rowraw)
        vals <- row[!is.na(row)]
        if(length(vals)<1){
            return(1)
        }
        return(pbeta(min(vals), 1, length(vals)))
    })

    results <- data.frame(resmat, metapval)
    colnames(results) <- names

    if(verbose) message("done\n")

    return.list <- list(p.value.table = results, details = details)

    return(return.list)

}
