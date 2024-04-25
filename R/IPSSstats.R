#' IPSS Stats
#' 
#' Function to add IPSSM and IPSS R stats based on risk_category
#' 
#' @return A patient statistics \code{data.frame}.
#' 
#' @export
#' 
IPSSstats <- function(patientResult){
load('ipssm_stats.Rda')
load('ipssr_stats.Rda')

if ('IPSSMcat_mean' %in% names(patientResult)) {
  patientResult <- merge(patientResult,ipssm_stats,by.x='IPSSMcat_mean',by.y='risk_category')
  print('IPSS-M model statistics added based on risk category')
} else {
  print('IPSS-M model statistics not included')
}

if ('IPSSRcat' %in% names(patientResult)) {
  patientResult <- merge(patientResult,ipssr_stats,by.x='IPSSRcat',by.y='risk_category')
  print('IPSS-R model statistics added based on risk category')
} else {
  print('IPSS-R model statistics not included')
}

return(patientResult)

}