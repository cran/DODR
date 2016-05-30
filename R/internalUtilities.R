#' Check if
#'
#' @param val1 matrix containing measurements columwise from set1
#' @param val2 matrix containing measurements columwise from set2
#' @param times1 vector of time points for set1
#' @param times2 vector of time points for set2
#'
#' @return Error Text if val-matrices are not propper
#' @noRd
checkIntegrity <- function(val1, val2, times1, times2){
    dimV1 <- dim(val1)
    dimV2 <- dim(val2)

    if(length(dimV1) != 2)return("val1 is no 2 dimensional matrix")
    if(length(dimV2) != 2)return("val2 is no 2 dimensional matrix")

    if(dimV1[1] != length(times1)){
        return("nuber of rows in val1 doen't match number of times in times1")
    }
    if(dimV2[1] != length(times2)){
        return("nuber of rows in val2 doen't match number of times in times2")
    }

    if(dimV1[2] != dimV2[2]){
    return("number of samples in val1 don't match number of saples in val2")
    }
    return(NULL)
}

