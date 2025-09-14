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
#' c()
#' 
bSplineTerm <- function(basisVar,nbasis = 5,groupingVar = '', moderators = c(),dropLargest = 1){
  require(fda)
  
  
  if(length(unique(basisVar))/nbasis < 4){
    nbasis <- floor(length(unique(basisVar))/4)
    warning(paste0('There are a large number of bases relative to unique data points. Automatically reducing the number of bases to ',nbasis,'.\n') )
  }
  
  curBases <- fda::create.bspline.basis(c(min(basisVar,na.rm=T)
                                          ,max(basisVar, na.rm=T))
                                        ,nbasis = nbasis)
  
  matBases <- fda::getbasismatrix(basisVar,curBases)
  
    lastBasis <- ncol(matBases)-dropLargest

  dfBases <- data.frame(matBases[,1:lastBasis])
  
  colnames(dfBases) <- paste0('bspl_',nbasis,'_',signif(basisVar[apply(dfBases,2,which.max)],2))
  
  
  
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
