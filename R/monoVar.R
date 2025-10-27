
#' Create an ordered incremented variable
#' 
#' Creates an ordered factor from a numeric vector.
#'  
#' In situations where we expect a monotonic change over time with 
#' most change happening fairly rapidly, we could want to treat a continuous
#' variable as an ordered factor. This function puts the smallest \code{increment_size} values
#' into a single factor level, the next \code{increment_size*2} values into the next smallest
#' factor levels, the next \code{increment_size*3} into the third of the factor levels, the next
#' eight into the fourth, and so on. The last factor level is combined such 
#' that the largest values do not have a factor level with smaller members
#' than the second-smallest factor level.
#' 
#' Note that increment sizes larger than 1 are recommended. This is because, for a common application 
#' (i.e., a \code{brms} model), there should be fewer than 50 total levels.
#' 
#' 
#' 
#' 
#' @param continVar Numeric vector, such as trial number.
#' @param return_numeric The grouped incremented numeric variable can be returned. By default, however, this is turned into an ordered factor and then returned.
#' @param increment_size Size of the smallest level, and the incremental increase in size with each level. Defaults to \code{ceiling(length(unique(continVar))/500)+1}, which means that the increment size is 2 if there are fewer than 500 unique values in \code{continVar}.
#'
#' @export
#'
monoVar <- function(continVar, return_numeric = F, increment_size = 'auto'){
  
  # continVar <- c(1:5000,1:5000)

  if(increment_size == 'auto'){
    increment_size <- ceiling(length(unique(continVar))/500)+1
  }
    
  uObs <- sort(unique(continVar))
  
  nObs <- length(uObs)
  
  vect <- c()
  nVal <- increment_size
  
  while(length(vect) < nObs){
    vect <- c(vect,rep(length(unique(vect))+1 , nVal))
    nVal <- nVal + increment_size
  } ; rm(nVal)
  
  vect[vect == max(vect)] <- max(vect)-1
  vect <- vect[1:nObs]

  vect_native <- rep(NA,nObs)
  for(curVal in unique(vect)){
    vect_native[vect == curVal] <- floor(median(continVar[vect == curVal]))
  } ; rm(curVal)
  
  monoVar <- rep(NA,length(continVar))
  for(curInd in 1:length(continVar)){
    monoVar[curInd] <- vect_native[uObs==continVar[curInd]]
  } ; rm(curInd)
  
  if(return_numeric){
    return(monoVar)
  }
  
  nDigits <-  ceiling(log10(max(monoVar)))
  sprintCode <- paste0('%s%0',nDigits,'d')
  
  vect_ordered <- ordered(as.factor(sprintf(sprintCode,'m_',monoVar))
                          ,sprintf(sprintCode,'m_',sort(unique(monoVar),decreasing = T)))
  
  return(vect_ordered)
    
}