#' Logistic regression on HBD probability or FLOD score
#' 
#' @param x an atlas object
#' @param expl_var the explanatory variable 'FLOD' or 'pHBD'
#' @param covar_df a dataframe or a matrix containing covariates
#' @param covar covariates of interest such as 'age', 'sex' , ...
#' if missing, all covariates of the dataframe are considered
#' @param n.cores number of cores for parallelization calculation (default = 1)
#' @param test which test to run. Default is bilateral. Use \code{"right"} to test for deleterious effects.
#' @param run whether the fonction is called or not (default = FALSE)
#' @param phen.code phenotype coding :
#'        - 'R' : 0:control ; 1:case ; NA:unknown (default)
#'        - 'plink' : 1:control ; 2:case ; 0/-9/NA:unknown
#' if 'plink' the function automatically convert it to 'R' to run logistic regression
#' 
#' @export

glm.HBD <- function( x, expl_var, covar_df, covar, test = c("bilateral", "right", "left"), run = FALSE, phen.code) {
  
  if(class(x)[1] != "atlas")
    stop("Need an atlas")
  
  if (run) { 
    if (expl_var == 'pHBD') {
      # Recovery pHBD
      hbd <- x@HBD_recap
    } else if (expl_var == 'FLOD') {
      # Recovery FLOD
      hbd <- x@FLOD_recap
    } else {
      stop("Explanatory variable must be 'pHBD' or 'FLOD'")
    }
    
    # Recovery phenotype
    id <- sub(".*:", "" , row.names(hbd))
    id.index <- match ( id, x@bedmatrix@ped$id )
    pheno <- x@bedmatrix@ped$pheno [id.index]
    if (phen.code == 'plink') {
      pheno <- ifelse(pheno == 1, 0, ifelse(pheno == 2, 1, NA))# Translate phenotype
    }
    
    # Recovery chr, snps, pos_cM and pos_Bp 
    final <- get.positions(x)	
    
    # unadjusted 
    if (missing(covar_df)) {
      message("No covariates given for the analysis = unadjusted data. To use covariates import a dataframe.")
      message(paste0("Call : glm(formula = pheno ~ ",expl_var,"[,i])"))
      x@logisticRegression$unadj <- cbind(final, glm.HBD.0(pheno, matrix(1, length(pheno)), hbd, test)) 
      message("-----------> GLM on UNADJUSTED data Done \n")
    }
    
    # adjusted 
    else {
      if(class(covar_df)[1] != "matrix" & class(covar_df)[1] != "data.frame") 
        stop("Need a matrix or a dataframe of covariates")
      
      if (class(covar_df)[1] == "data.frame")
        covar_df <- as.matrix(covar_df)
      
      if(missing(covar)) {
        message(paste0("No covariates specified - All covariates of the dataframe will be used : " , gsub(",", " +", toString(colnames(covar_df)))))
        message(paste0("Call : glm(formula = pheno ~ ",expl_var,"[,i] + ", gsub(",", " +", toString(colnames(covar_df))) ,")" ))
        df <- na.omit(covar_df[id,])				 # take all covar given in the dataframe
      } else {
        message(paste0("Covariates = ", gsub(",", " +", toString(covar))))
        message(paste0("Call : glm(formula = pheno ~ ",expl_var,"[,i] + ", gsub(",", " +", toString(covar)) ,")"))
        df <- na.omit(covar_df[ id , covar]) #rownames covar_df  = individual id 	
      }
      x@logisticRegression$adj <- cbind(final, glm.HBD.0(pheno, cbind(1,df), hbd, test)) 
      message("-----------> GLM on ADJUSTED data Done \n")
    }
  }
  x
} 
