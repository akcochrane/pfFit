#' Fit a model to cumulative data
#' 
#' Unimplemented
#'
#' @param unimplemented 
#'
#' @export
#'
#' @examples
#' 
#' c()
#' 
fit_cumulative <- function(unimplemented){
  
  # # maybe:
  # Combining several ordinal measures in clinical studies 10.1002/sim.1778
  #
  # # ideally:
  # make a demo for multivariate `cumulative()` model, then using participant-level
  # ranefs to get a PC dominant dimension
  # OR
  # if all scales have the same number of options, have the data long, have 
  # scale be a factor, estimate the differences in the thresholds, and
  # have participant-level ranefs be extracted for the scale-agnostic scores
  
m <- brm(
    bForm
    ,data = dat
    ,family = cumulative()
    ,prior = priorIn
    ,...
  )
  
  try({
    m$fitted_param$xxxx <- fitted(m, ndraws = 200)
  })
  
  return(m)
  
}