#' Annotate IPSS-M result
#'
#' Consolidates the IPSS-M results to create single IPSSM score and IPSSM category columns. If the maximal range of IPSS-M risk score (worst-best) is below the allowed limit, then the IPSS-M results under the mean scenario are reported. Of course if no missing data, the IPSS-M values are fully certains.
#'
#' @param patientResult a patient result \code{data.frame}, as the output of \code{IPSSMmain} function.
#' @param range.max threshold for the allowed maximal range in IPSS-M risk (worst - best) to confidently report the mean scenario as the main result. Default 1.

#'
#' @return An annotated patient \code{data.frame}, same number of rows/patients as in \code{patientResult}, with additonal columns labelled \code{IPSSMscore}, \code{IPSSMcat}, and \code{Comment}.

#' @export
#'
#' @examples
#' dd <- read.csv(system.file("extdata", "IPSSMexample.csv", package = "ipssm"), header = TRUE)
#' dd.process <- IPSSMprocess(patientInput = dd)
#' dd.res <- IPSSMmain(patientProcess = dd.process)
#' dd.annot <- IPSSMannotate(patientResult = dd.res)
#' print(dd.annot[, c("ID", "IPSSMscore", "IPSSMcat", "Comment")])
#'
IPSSMannotate <- function(patientResult, range.max = 1) {
  cat("Annotating IPSS-M results\n")
  
  # IPSSM score and cat unified results
  patientResult$IPSSMscore <- patientResult$IPSSMscore_mean
  patientResult$IPSSMcat <- as.vector(patientResult$IPSSMcat_mean)
  patientResult$Comment <- NA
  
  # Start by filling in for cases where the distance between best and worst IPSS-M risk score
  # is lower than the allowed maximan range
  i <- which((patientResult$IPSSMscore_worst - patientResult$IPSSMscore_best <= range.max)&
               as.vector(patientResult$IPSSMcat_worst) == as.vector(patientResult$IPSSMcat_best))
  if (length(i) > 0) {
    patientResult$Comment[i] <- "range (worst - best) below limit"
  }

  # Continue by filling in for cases where the distance between best and worst IPSS-M risk score
  # is lower than the allowed maximan range but categories are different
  ii <- which((patientResult$IPSSMscore_worst - patientResult$IPSSMscore_best <= range.max)&
               as.vector(patientResult$IPSSMcat_worst) != as.vector(patientResult$IPSSMcat_best))
  if (length(ii) > 0) {
    patientResult$Comment[ii] <- "IPSS-M results less confident due to missing mutation data"
  }
 
  # Continue with cases where the range is above limit
  # But the categories remain unchanged
  j <- which((patientResult$IPSSMscore_worst - patientResult$IPSSMscore_best > range.max) &
               as.vector(patientResult$IPSSMcat_worst) == as.vector(patientResult$IPSSMcat_best))
  if (length(j) > 0) {
    patientResult$Comment[j] <- "range (worst - best) above limit, but stable category"
  }
  
  # Continue with cases where the range is above limit
  # Return warning when the categories change
  jj <- which((patientResult$IPSSMscore_worst - patientResult$IPSSMscore_best > range.max) &
                as.vector(patientResult$IPSSMcat_worst) != as.vector(patientResult$IPSSMcat_best))
  if (length(jj) > 0) {
    patientResult$Comment[jj] <- "IPSS-M results less confident due to missing mutation data"
  }
  
  
  # Change the label when no uncertainty
  k <- which(patientResult$IPSSMscore_worst - patientResult$IPSSMscore_best == 0)
  if (length(k) > 0) {
    patientResult$Comment[k] <- "no uncertainty (no missing data)"
  }
  
  return(patientResult)
}