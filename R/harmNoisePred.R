#' Asymmetric variant of Chow method applied on harmonic regression
#'
#' @description A harmonic regression is done on the first set of time
#' series (val1) pair. The differences of the second series to this fit are
#' calculated. By comparison of these distances with the noise estimation of
#' the first series, a decision is made whether the second series could be
#' explained as additional samples of the first series.
#'
#' @inheritParams dodr
#'
#' @return list containing the fits for both time series and the combination
#' and the pValue for differential Oscillation
#' @export
#' @aliases harmNoisePred1 harmNoisePred2
#' @import Matrix
harmNoisePred <- function (val1, val2,
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
            return(series/mean(series, na.rm=TRUE))
        })
        val2 <- apply(val2, 2, function(series){
            return(series/mean(series, na.rm=TRUE))
        })
    }





    has.na <- !all(!is.na(val1))|!all(!is.na(val2))

    if(has.na){
        resset <- do.call(rbind, lapply(seq_len(ncol(val1)), function(set){

            num1 <- !is.na(val1[,set])
            num2 <- !is.na(val2[,set])

            if(sum(num1)<4 | sum(num2)<4){
                return(c(NA, NA))
            }

            mat1 <- try(matrixPreparation(times1[num1], period), TRUE)
            mat2 <-try( matrixPreparation(times2[num2], period), TRUE)
            if(class(mat1)=="try-error" | class(mat2)=="try-error"){
                return(c(NA, NA))
            }

            #apply linear fitting on first time series
            fit1 <- linFit(mat1, val1[num1,set])

            #calculate difference matrix
            d <- val2[num2,set] - mat2$X %*% fit1$b
            I <- as.matrix(Diagonal(nrow(mat2$X)))
            uni <- matrix(1,ncol=1,nrow=nrow(mat2$X))

            score <- matrix(1,nrow=1, ncol=ncol(d)) %*%
                ((t(d) %*% solve(I+mat2$X %*% mat1$CX %*% t(mat2$X)) %*% d) *
                     as.matrix(Diagonal(ncol(d)))) /
                (fit1$ssr/(fit1$df)*nrow(mat2$X))

            pval <- 1-pf(score[1,],nrow(d),fit1$df)

            return(c(pval, score[,1]))
        }))

        resultdf <- data.frame(resset)
        colnames(resultdf) <- c('p.value', 'F')
        return(resultdf)
    } else{

        #Prepare matrices for linear fitting

        mat1 <- matrixPreparation(times1, period)
        mat2 <- matrixPreparation(times2, period)

        indexSeries <- seq_len(ncol(val1))

        samplesets <- split(indexSeries,
                            cut(indexSeries, seq(0,max(indexSeries)+100,100),)
        )

        resset <- do.call(rbind, lapply(samplesets, function(set){
            #apply linear fitting on first time series
            fit1 <- linFit(mat1, val1[,set])

            #calculate difference matrix
            d <- val2[,set] - mat2$X %*% fit1$b
            I <- as.matrix(Diagonal(nrow(mat2$X)))
            uni <- matrix(1,ncol=1,nrow=nrow(mat2$X))

            score <- matrix(1,nrow=1, ncol=ncol(d)) %*%
                ((t(d) %*% solve(I+mat2$X %*% mat1$CX %*% t(mat2$X)) %*% d) *
                     as.matrix(Diagonal(ncol(d)))) /
                (fit1$ssr/(fit1$df)*nrow(mat2$X))

            pval <- 1-pf(score[1,],ncol(d),fit1$df)

            return(cbind(pval, score[1,]))
        }))

        resultdf <- data.frame(resset)
        colnames(resultdf) <- c('p.value', 'F')
        return(resultdf)
    }




}
