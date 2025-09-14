
#' Define a psychometric function formula
#'
#' @param yVar The variable to predict in the psychometric function; probably a binary numeric (i.e., 0 or 1).
#' @param pf Character. The name of the psychometric function you want (i.e., "logistic" or "weibull"), or "identity" for no psychometric function
#' @param link_x The linking \code{x}
#' @param ... 
#'
#' @export
#'
#' @examples
#' 
#' c()
#' 
pfFormula <- function(yVar, pf = c('identity','logistic','weibull')
                      ,link_x = 'LINK_X_VARIABLE'
                      ,...){
  
  ## THIS DOES NOT SUCCESSFULLY USE EXISTS()
  
  
  ## also need to figure out a seemless way of augmenting the data with the bases
  
  curEnv <- environment()
 
  # print(curEnv)
  
    pfCharacter <- gsub('LINK_X_VARIABLE',link_x,pfChar(pf))
    
    pfForm <- list()
    
    pfForm[['pf']] <- formula(paste0(yVar, ' ~ ',pfCharacter))
    
    ## define threshold
    if(exists('threshold')){
      pfForm[['threshold']] <- threshold
    }else{
      pfForm[['threshold']] <- threshold ~ 1
    }
    
    if(pf == 'logistic'){
      ## define bias
      if(exists('bias')){
        pfForm[['bias']] <- bias
      }else{
        pfForm[['bias']] <- bias ~ 1
      }
      ## define lapse
     if(exists('lapseRate', inherits = F)){
      pfForm[['lapseRate']] <- lapseRate
    }else{
      pfForm[['lapseRate']] <- .01
    }
    }
  
    if(pf == 'weibull'){
      ## define shape
      if(exists('shape')){
        pfForm[['shape']] <- shape
      }else{
        pfForm[['shape']] <- shape ~ 1
      }
      ## define lapse
      if(exists('rhAsymptote')){
        pfForm[['rhAsymptote']] <- rhAsymptote
      }else{
        pfForm[['rhAsymptote']] <- .99
      }
    }
    return(pfForm)
}
