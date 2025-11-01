
#' Make a data frame of b-spline bases
#'
#' @param returnDF If TRUE, return a data frame rather than matrix
#' @inheritParams bSplineTerm
#'
#' @export
#' 
#' @examples
#' 
#' d <- data.frame(trial_number = 1:100)
#' 
#' # get a variable with a b-spline with 5 bases, but dropping the one with the mode centered on the largest trial number
#' d$spl5dr1 <- bSplineMat(d$trial_number)
#' 
#'  # get a variable with a b-spline with 9 bases, but dropping the three with the modes centered on the largest trial numbers
#' d$spl9dr3 <- bSplineMat(d$trial_number,nbasis = 9, dropLargest = 3)
#' 
bSplineMat <- function(basisVar, nbasis = 5,dropLargest = 1,returnDF = F){
  
  require(fda)
  
  if(length(unique(basisVar))/nbasis < 4){
    nbasis <- floor(length(unique(basisVar))/4)
    warning(paste0('There are a large number of bases relative to unique data points. Automatically reducing the number of bases to ',nbasis,'.\n') )
  }
  
  curBases <- create.bspline.basis(c(min(basisVar,na.rm=T)
                                     ,max(basisVar, na.rm=T))
                                   ,nbasis = nbasis)
  
  matBases <- getbasismatrix(basisVar,curBases)
  
  lastBasis <- ncol(matBases)-dropLargest
  
  dfBases <- data.frame(matBases[,1:lastBasis])
  
  colnames(dfBases) <- paste0('bspl_',nbasis,'_',basisVar[apply(dfBases,2,which.max)])
  
  if(returnDF){
  return(dfBases)
  }else{
    return(as.matrix(dfBases))
  }
  
}