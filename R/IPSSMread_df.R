#' Read User Input Data accepts input from data frame
#'
#' @param patientInput a patient \code{data.frame}, one patient per row, and variables as columns.
#' @param genesRes a \code{vector} containing the names of the residual genes.
#' @param maxvafloh maximum variant allele frequency (vaf) threshold to determine theres Loss of Heterozigosity (LOH).
#' @param Nref the average reference value for min(Nres,2) where Nres is the number of mutated residual genes.
#'
#' @return A processed patient \code{data.frame}, same number of rows/patients as in \code{patientInput}, and with the processed variables as additional columns.
#'
#' @export
#'
IPSSMread_df <- function(df,calculation='ALL') {
  #Calculation Parameter:
  #ALL Calculate IPSS-M IPSS-R and IPSS-RA scores 
  #IPSS-M Calculate IPSS-M only ANC and AGE are not needed
  #IPSS-RA Calculate IPSS-RA only AGE and ANC are needed Genetic Factors are Not
  #IPSS-R Calculate IPSS-R only no AGE adjustment AGE not needed
  
  d <- df
  
  #load field name mappings
  UPPER <- c("ID","AGE","ANC","HB","PLT","BM_BLAST","DEL5Q","DEL7_7Q","COMPLEX","CYTO_IPSSR","DEL17_17P",
  "TP53MUT","TP53MAXVAF","TP53LOH","MLL_PTD","FLT3","ASXL1","BCOR","BCORL1","CBL","CEBPA","DNMT3A","ETV6",
  "EZH2","IDH1","IDH2","KRAS","NF1","NPM1","NRAS","RUNX1","SETBP1","SF3B1","SRSF2","STAG2","U2AF1","ETNK1",
  "GATA2","GNB1","PHF6","PPM1D","PRPF8","PTPN11","WT1")
  
  CASE <- c("ID","AGE","ANC","HB","PLT","BM_BLAST","del5q","del7_7q","complex","CYTO_IPSSR","del17_17p",
  "TP53mut","TP53maxvaf","TP53loh","MLL_PTD","FLT3","ASXL1","BCOR","BCORL1","CBL","CEBPA","DNMT3A","ETV6",
  "EZH2","IDH1","IDH2","KRAS","NF1","NPM1","NRAS","RUNX1","SETBP1","SF3B1","SRSF2","STAG2","U2AF1","ETNK1",
  "GATA2","GNB1","PHF6","PPM1D","PRPF8","PTPN11","WT1")
  
  flds <- data.frame(UPPER,CASE)
  
  #convert field to upper lookup proper case for input and rename column
  for (x in colnames(d)) {
    colnames(d)[colnames(d) == x] = flds[flds$UPPER == toupper(x),'CASE']
  }
  
  #Add CALCULATION to records
  d$CALCULATION = calculation
  
  #Add UUID to records
  #d$UUID = suppressWarnings(apply(d,1,UUIDgenerate))
  
  #print("Did you make sure that the unit for Hemoglobin (HB) is g per dL. We expect values between 4 and 20 g per dL.")
  #print("Did you make sure that the unit for Hemoglobin (HB) is g per dL.")
  #print("Did you make sure that the unit for Platelets (PLT) is Giga per L.")
  #print("Did you make sure that the unit for Bone Marrow Blast (BM_BLAST) is percentage.")
  
  # ~~~~~~~~~ Variables Validation: | Presence | Format (numerical binary etc) | Range
  check.failed.numerical <- function(x) {
    return(any(!((is.numeric(x) | is.na(x)))))
  }
  check.failed.binary <- function(x) {
    return(any(!((x %in% c(0, 1) | is.na(x)))))
  }
  check.failed.char.given <- function(x, mychar) {
    return(any(!((x %in% mychar | is.na(x)))))
  }
  
  # ~~~ CLINICAL
  
  if (calculation == 'IPSS-M') {
    clinical.var <- c("HB", "PLT", "BM_BLAST") #ANC and AGE not required
  } else if (calculation == 'IPSS-R') {
    clinical.var <- c("HB", "PLT", "BM_BLAST", "ANC") #ANC required AGE not required
  } else if (calculation == 'IPSS-RA' | calculation == "ALL") {
    clinical.var <- c("HB", "PLT", "BM_BLAST", "ANC", "AGE") #ANC and AGE required
  }
  
  if (any(!clinical.var %in% colnames(d))) {
    stop(paste(paste(clinical.var, collapse = " "), "should be columns of your input data"))
  }
  for (cc in clinical.var) {
    if (check.failed.numerical(d[, cc])) {
      stop(paste(cc, "should have numerical values"))
    }
  }
  for (cc in clinical.var) {
    if (any(d[, cc] < 0, na.rm = T)) {
      stop(paste(cc, "contains negative values. Please re-check your input data!"))
    }
  }
  if (any(d[, "HB"] < 4, na.rm = T)) {
    warning("You have HB values smaller than 4 g/dL. This is suspicious.")
  }
  if (any(d[, "HB"] > 20, na.rm = T)) {
    warning("You have HB values larger than 20 g/dL. This is suspicious.")
  }
  
  if (any(d[, "PLT"] > 2000, na.rm = T)) {
    warning("You have PLT values larger than 2000 1e9/L. This is suspicious.")
  }
  
  for (x in 1:nrow(d)) {
    if (d[x,"BM_BLAST"] > 0 & d[x,"BM_BLAST"] < 1){
      print("You have BM_BLAST in decimal form, input was changed to percentage")
      d[x,"BM_BLAST"] = d[x,"BM_BLAST"] * 100
    }
  }
  
  if (any(d[, "BM_BLAST"] > 30, na.rm = T)) {
    warning("You have BM_BLAST values larger than 30 %. This is suspicious.")
  }
  
  #AGE only needed for IPSS-RA or ALL calculations
  if (calculation == 'IPSS-RA' | calculation == 'ALL'){
    if (any(d[, "AGE"] < 18, na.rm = T)) {
      warning("You have AGE values smaller than 18 Years. This is suspicious.")
    }
    if (any(d[, "AGE"] > 120, na.rm = T)) {
      warning("You have AGE values larger than 120 Years. This is suspicious.")
    }
  }
  
  if (calculation != 'IPSS-M'){
    if (any(d[, "ANC"] > 15, na.rm = T)) {
      warning("You have ANC values larger than 15 1e9/L. This is suspicious.")
    }
  }
  
  # ~~~ CYTO IPSSR Categorical
  cyto.var <- c("CYTO_IPSSR" # cytogenetics categorical: Very Good, Good, Intermediate, Poor, Very Poor
  )
  if (any(!cyto.var %in% colnames(d))) {
    stop(paste(paste(cyto.var, collapse = " "), "should be columns of your input data"))
  }
  
  if (check.failed.char.given(as.vector(d[, "CYTO_IPSSR"]), c("Very Good", "Good", "Intermediate", "Poor", "Very Poor"))) {
    stop(paste("CYTO_IPSSR", "should have values from Very Good, Good, Intermediate, Poor, Very Poor"))
  }
  
  # Validate only where IPSS-M or ALL are Calculated
  if (calculation == 'IPSS-M' | calculation == 'ALL'){
    # ~~~ CYTO
    
    cyto.var <- c(
      "del5q", "del7_7q", "del17_17p", "complex" # cytogenetics: 0/1
    )
    if (any(!cyto.var %in% colnames(d))) {
      stop(paste(paste(cyto.var, collapse = " "), "should be columns of your input data"))
    }
    for (cc in cyto.var[1:4]) {
      if (check.failed.binary(d[, cc])) {
        stop(paste(cc, "should have binary 0/1 values"))
      }
    }
    
    # TODO: ADD VALIDATION TEST For del5q/CYTO etc
    
    # ~~~ TP53
    tp53.var <- c(
      "TP53mut", # TP53 number of mutations categorical: 0, 1, 2 or more
      "TP53maxvaf", # continuous %: default 0
      "TP53loh" # TP53 LOH: 0/1/NA
    )
    if (any(!tp53.var %in% colnames(d))) {
      stop(paste(paste(tp53.var, collapse = " "), "should be columns of your input data"))
    }
    d$TP53mut <- as.character(d$TP53mut)
    if (check.failed.char.given(as.vector(d[, "TP53mut"]), c("0", "1", "2 or more"))) {
      stop(paste("TP53mut", "should have values from 0, 1, 2 or more"))
    }
    if (check.failed.numerical(d[, "TP53maxvaf"])) {
      stop(paste("TP53maxvaf", "should have numerical values"))
    }
    if (any(d[, "TP53maxvaf"] < 0, na.rm = T)) {
      stop(paste("TP53maxvaf", "contains negative values. Please re-check your input data!"))
    }
    
    for (x in 1:nrow(d)) {
      if (d[x,"TP53maxvaf"] >= 1 & d[x,"TP53maxvaf"] <= 100){
        print("You have TP53maxvaf in percentage form, input was changed to decimal")
        d[x,"TP53maxvaf"] = d[x,"TP53maxvaf"] / 100
      }
    }
    
    if (any(d[, "TP53maxvaf"] > 1, na.rm = T)) {
      stop(paste("TP53maxvaf", "contains values>1. We expect values between 0 and 1."))
    }
    
    if (check.failed.binary(d[, "TP53loh"])) {
      stop(paste("TP53loh", "should have binary 0/1 values"))
    }
    
    # ~~~ GENES
    gene.var <- c(
      "MLL_PTD", # all main effect genes below: 0/1/NA
      "FLT3",
      "ASXL1",
      "CBL",
      "DNMT3A",
      "ETV6",
      "EZH2",
      "IDH2",
      "KRAS",
      "NPM1",
      "NRAS",
      "RUNX1",
      "SF3B1",
      "SRSF2",
      "U2AF1",
      "BCOR", # all residual genes below: 0/1/NA
      "BCORL1",
      "CEBPA",
      "ETNK1",
      "GATA2",
      "GNB1",
      "IDH1",
      "NF1",
      "PHF6",
      "PPM1D",
      "PRPF8",
      "PTPN11",
      "SETBP1",
      "STAG2",
      "WT1"
    )
    if (any(!gene.var %in% colnames(d))) {
      stop(paste(paste(gene.var, collapse = " "), "should be columns of your input data"))
    }
    for (cc in gene.var) {
      if (check.failed.binary(d[, cc])) {
        stop(paste(cc, "should have binary 0/1 values"))
      }
    }
  }
  
  print("Data successfully imported.")
  
  return(d)
}