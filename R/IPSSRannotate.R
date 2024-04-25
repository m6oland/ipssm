#' Annotate IPSS-R result
#'
#' Consolidates the IPSS-R results to create single IPSSM score and IPSSR category columns. If the maximal range of IPSS-R risk score (worst-best) is below the allowed limit, then the IPSS-R results under the mean scenario are reported. Of course if no missing data, the IPSS-R values are fully certain.
#'
#' @param patientResult a patient result \code{data.frame}, as the output of \code{IPSSMmain} function.
#' @param range.max threshold for the allowed maximal range in IPSS-R risk (worst - best) to confidently report the mean scenario as the main result. Default 1.

#'
#' @return An annotated patient \code{data.frame}, same number of rows/patients as in \code{patientResult}, with additonal columns labelled \code{IPSSRscore}, \code{IPSSRcat}, and \code{Comment}.

#' @export
#'

IPSSRannotate <- function(patientResult,calculation='ALL') {
  #Calculation Parameter:
  #ALL Calculate IPSS-M IPSS-R and IPSS-RA scores 
  #IPSS-M Calculate IPSS-M only ANC and AGE are not needed
  #IPSS-RA Calculate IPSS-RA only AGE and ANC are needed Genetic Factors are Not
  #IPSS-R Calculate IPSS-R only no AGE adjustment AGE not needed
  
  cat("Annotating IPSS-R results\n")
  
  if (!'CYTOVEC' %in% colnames(patientResult)) {
    # Cytogenetics as a numerical vector
    patientResult$CYTOVEC <- car::recode(patientResult$CYTO_IPSSR,
       "'Very Good'=0; 'Good'=1; 'Intermediate'=2; 'Poor'=3; 'Very Poor'=4",
       as.factor = FALSE
    )
  }
  
  # IPSSM score and cat unified results
  patientResult$IPSSRscore <- NA
  patientResult$IPSSRcat <- NA

  # IPSSR Score for Clinical Data
  patientResult$HB_RS <- as.numeric(as.character(cut(patientResult$HB,breaks=c(-Inf,8,10,Inf),labels=c('1.5','1','0'),right = FALSE)))
  patientResult$ANC_RS  <- as.numeric(as.character(cut(patientResult$ANC,breaks=c(-Inf,0.8,Inf),labels=c('0.5','0'),right = FALSE)))
  patientResult$BM_BLAST_RS  <- as.numeric(as.character(cut(patientResult$BM_BLAST,breaks=c(-Inf,2,4.99999,10,Inf),labels=c(0,1,2,3),right = TRUE)))
  patientResult$PLT_RS <- as.numeric(as.character(cut(patientResult$PLT,breaks=c(-Inf,50,100,Inf),labels=c(1,0.5,0),right = FALSE)))
  
  # SUM IPSSR Individual Scores to get Total IPSSR Score
  patientResult$IPSSRscore <-  patientResult$BM_BLAST_RS + patientResult$CYTOVEC + patientResult$HB_RS + patientResult$PLT_RS + patientResult$ANC_RS 
  
  #Categorize Results based on IPSSR
  patientResult$IPSSRcat <- cut(patientResult$IPSSRscore,breaks=c(-Inf,1.5,3,4.5,6,Inf),labels=c('Very Low','Low','Intermediate','High','Very High'),right = TRUE)

  
  if(calculation=='ALL' | calculation=='IPSS-RA'){
    # Adjust Score for AGE
    patientResult$IPSSRAadjustment <- (patientResult$AGE - 70)*(abs(0.05-(patientResult$IPSSRscore*0.005)))
    
    # Calculate IPSSRA Score by adjusting IPSSR Score
    patientResult$IPSSRAscore <- patientResult$IPSSRscore + patientResult$IPSSRAadjustment
    
    #Categorize Results based IPSSRA scores
    patientResult$IPSSRAcat <- cut(patientResult$IPSSRAscore,breaks=c(-Inf,1.5,3,4.5,6,Inf),labels=c('Very Low','Low','Intermediate','High','Very High'),right = TRUE)
    
  }
    
  return(patientResult)
}