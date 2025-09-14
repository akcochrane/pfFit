
#' Get a psychometric function formula in character form
#'
#' @param pfType The type of psychometric function. Currently supported are \code{logistic} or \code{weibull}
#'
#'
#' @export
#'
#' @examples
#' 
#' c()
#' 
pfChar <- function(pfType){
  
  thresholdVal <- .75 # need to have the bases of the exponents be defined by user
  
  if(pfType == 'logistic'){
    pf <- paste0("lapseRate + (1-2*lapseRate)/(1+3^((bias-LINK_X_VARIABLE)/threshold))")
    attr(pf,'varNames') <- c('threshold','bias','lapseRateA','lapseRateB')
    
  }
  
  if(pfType == 'weibull'){
    pf <- "xIntercept+(rhAsymptote - xIntercept)*(1-2^(-(LINK_X_VARIABLE/threshold)^shape)) " 
    attr(pf,'varNames') <- c('threshold','shape','xIntercept','rhAsymptote')
  }
  
  attr(pf,'thresholdVal') <- thresholdVal
  return(pf)
}
