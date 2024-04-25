#' IPSS Survival
#' 
#' Function to add IPSSM and IPSS R survival stats based on calculated risk categories
#' 
#' @return A patient survival rates \code{data.frame}.
#' 
#' @export
#' 

IPSSsurvival <- function(patientResult){
  load('ipssm_surv.Rda')
  load('ipssr_surv.Rda')
  
  if ('IPSSMcat_mean' %in% names(patientResult)) {
    patientResult <- merge(patientResult,ipssm_surv,by.x='IPSSMcat_mean',by.y='risk_category')
    print('IPSS-M model survival summary statistics added based on risk category')
  } else {
    print('IPSS-M model survival summary statistics not included')
  }
  
  if ('IPSSRcat' %in% names(patientResult)) {
    patientResult <- merge(patientResult,ipssr_surv,by.x='IPSSRcat',by.y='risk_category')
    print('IPSS-R model survival summary statistics added based on risk category')
  } else {
    print('IPSS-R model survival summary statistics not included')
  }
  
  return(patientResult)
  
}