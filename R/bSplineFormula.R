
#' To a formula, add a set of terms corresponding to a set of b-spline basis functions
#' 
#' Construct a set of bases using \code{basisVar} and add the terms from these bases
#' to the formula \code{formIn}. Relies on \code{\link[fda]{create.bspline.basis}} and
#' \code{\link[fda]{getbasismatrix}} to construct bases.
#' 
#' Given the numeric vector \code{basisVar} and the number of bases \code{nbasis},
#' an attribute \code{attr(output,'basisDF')} is returned giving the density [weight]
#' for each of the bases at each of the observations.
#' 
#' To increase identifiability in common use cases, such as within conventional mixed-effects models,
#' by default the "largest" bases are removed (e.g., the basis with the mode at the largest
#' value of \code{basisVar} is deleted.). If identifiability is still posing difficulties,
#' \code{dropLargest} could be increased, in order to remove an increasing
#' number of the bases with modes nearest the largest values of \code{basisVar}. Alternatively, if the 
#' full coverage is desired, \code{dropLargest} should be lowered to 0 in order to not
#' remove any bases.
#'
#' @param formIn Formula,
#' @inheritParams bSplineTerm
#'
#' @export
#'
#' @examples
#' 
#' c()
#' 
bSplineFormula <- function(formIn,basisVar
                           ,groupingVar = '' 
                           ,nbasis = 5
                           ,nbasis_groupingVar = -1
                           ,dropLargest = 2
                           ,dropLargest_groupingVar = -1){
  
  formChar <- as.character(formIn)
  lhs <- formChar[2]
  rhs <- formChar[3]
  
  
  splineTerm <- bSplineTerm(basisVar = basisVar
                            ,nbasis = nbasis
                            ,dropLargest = dropLargest)
  if(groupingVar != ''){ # with a grouping variable
    
    if(nbasis_groupingVar<0){
      nbasis_groupingVar=nbasis
    }    
    if(dropLargest_groupingVar<0){
      dropLargest_groupingVar=dropLargest
    }
    
    splineTerm_group <-  bSplineTerm(basisVar = basisVar
                                     ,groupingVar = groupingVar
                                     ,nbasis = nbasis
                                     ,dropLargest = dropLargest_groupingVar)
    
    formOut <- formula(paste0(lhs,' ~ ',rhs, ' + ',splineTerm,' + ',splineTerm_group))
    
    attr(formOut,'fda_package_object_groupingVar') <- attr(splineTerm_group,'fda_package_object')
    
    basisDF_grouping <- attr(splineTerm_group,'basisDF')
    
    if(all(dim(basisDF_grouping)==dim(attr(splineTerm,'basisDF')))){
      attr(formOut,'basisDF') <- attr(splineTerm,'basisDF')
    }else{
      
      basisDF_grouping <- data.frame(attr(splineTerm,'basisDF'),basisDF_grouping)
      attr(formOut,'basisDF') <- basisDF_grouping
    }
  }else{ # without a grouping variable
    formOut <- formula(paste0(lhs,' ~ ',rhs, ' + ',splineTerm))
  }
  
  attr(formOut,'fda_package_object') <- attr(splineTerm,'fda_package_object')
  
  return(formOut)
}
