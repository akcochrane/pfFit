#' Define a b-spline term for a model formula
#'
#' @param basisVar Character. The name of the numeric variable to construct the basis over.
#' @param groupingVar Character. Optional. If provided, return a string formatted for a random effects structure with no intercept.
#' @param nbasis Number of bases. Note that, in most cases, the returned number will be smaller than this due to \code{dropLargest}.
#' @param moderators Character vector. Optional. If provided, names of variables to interact with bases.
#' @param dropLargest Number of bases to drop, to improve identifiability.
#'
#' @export
#'
#' @examples
#' 
#' splineTerm <- bSplineTerm(basisVar = 1:200)
#' 
bSplineTerm <- function(basisVar,nbasis = 5,groupingVar = '', moderators = c(),dropLargest = 1){
  require(fda)
  
  dfBases <- bSplineDF(basisVar, nbasis)
  
  if(length(moderators)>0){
    allCombs <- expand.grid(moderators,colnames(dfBases))
    basisTerm <- paste(allCombs[,1],allCombs[,2],sep='*',collapse=' + ')
    rm(allCombs)
  }else{
    basisTerm <- paste(moderators,colnames(dfBases),sep='',collapse=' + ')
  }
  
  if(nchar(groupingVar)>0){
    basisTerm <- paste0('(0 + ',basisTerm,' || ',groupingVar,')')
  }
  
  attr(basisTerm,'fda_package_object') <- curBases
  attr(basisTerm,'basisDF') <- dfBases
  attr(basisTerm,'basisVar') <- basisVar
  
  return(basisTerm)
}
