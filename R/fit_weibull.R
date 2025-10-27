#' Fit a Weibull psychometric function generalized mixed-effects model
#' 
#' Fit a model to binary response data, where the link between positively-bounded \code{x}
#' and the binary \code{y} is a Weibull function. The Weibull function is defined
#' with \code{threshold} and \code{shape} parameters, which correspond broadly to
#' a "location" and "slope."
#' 
#' This function serves to set up a nonlinear \code{\link[brms]{brm}} model with a 
#' bernoulli response distribution.
#'
#' @param y The binary variable to be predicted. It should change as a monotonic function of \code{x}.
#' @param x Positive-bounded variable used as the "x-axis" of the psychometric function
#' @param thresholdForm A formula, either a linear model or linear mixed-effects model, specifying the predicted values of \code{scaleLogisThreshold}. This is the scaled inverse-logit-transformed threshold, which provides more stable estimation than sampling the distributions of the raw threshold.
#' @param shapeForm A formula, either a linear model or linear mixed-effects model, specifying the predicted values of \code{logShape}. This is the log-transformed shape, which provides more stable estimation than sampling the distributions of the raw shape.
#' @param dat Data frame with variables with which to fit the model.
#' @param xIntercept Scalar number: What should the value of the Weibull function be when \code{x} is zero? 
#' @param rhAsymptote Scalar number: What should the value of the Weibull function be when \code{x} is infinite? 
#' @param range_x Two numbers (e.g., \code{c(.1,5)}). Defines the plausible range of \code{threshold} values, which act to limit the possible values of the threshold. Defaults to the \code{min} and \code{max} of \code{x}, but if the user wants their threshold to have different boundaries, they can be user-specified.
#' @param ... 
#'
#' @export
#'
#' @examples
#' 
#'     
#'     d <- TEfits::anstrain
#'     d$absRat <- abs(d$ratio)
#'     
#'     d$spl <- bSplineMat(basisVar = d$trialNum)
#'                                    
#'     m_weibull_linearFixed <- fit_weibull(y='acc',x='absRat'
#'                                    , thresholdForm = scaleLogisThreshold ~ scale(trialNum) + (spl || subID)
#'                                    , shapeForm = logShape ~ (1 | subID)
#'                                    , dat = d
#'                                    , cores = 2 , chains = 2 # only run 2 chains for efficiency
#'                                    )
#'                                    
#'     d$trialM <- monoVar(d$trialNum)
#'     
#'     m_weibull_monoFixed <- fit_weibull(y='acc',x='absRat'
#'                                    , thresholdForm = scaleLogisThreshold ~ mo(trialM) + (spl || subID)
#'                                    , shapeForm = logShape ~ (1 | subID)
#'                                    , dat = d
#'                                    , cores = 2 , chains = 2 # only run 2 chains for efficiency
#'                                    )
#'                                    
#' 
fit_weibull <- function(y
                        , x
                        , thresholdForm = scaleLogisThreshold ~ 1
                        , shapeForm = logShape ~ 1
                        , dat
                        , xIntercept = .5
                        , rhAsymptote = .99
                        , range_x = 'auto'
                        , ...){
  require(brms)
  
  if(!is.null(attr(thresholdForm,'basisDF'))){
    dat <- data.frame(dat,attr(thresholdForm,'basisDF'))
  }
  
  if(length(range_x)==1){
    range_x <- signif(range(d[,x],na.rm = T),2)
  }
  
  threshForm_scale <- formula(paste0('threshold ~ ',range_x[1],' + (',range_x[2],'-',range_x[1],')*inv_logit(scaleLogisThreshold)'))
  
  bForm <- brmsformula(
    pfFormula(y,pf = 'weibull'
              ,link_x = x
    )$pf
    ,thresholdForm
    ,shapeForm
    ,nl = T
  ) +
    nlf(threshForm_scale) + 
    nlf(shape ~ exp(logShape)) +
    nlf(formula(paste0('xIntercept ~ ',xIntercept))) +
    nlf(formula(paste0('rhAsymptote ~ ',rhAsymptote)))
  
  m <- brm(bForm
           ,data = dat
           ,data2 = dat
           ,family = bernoulli(link = 'identity')
           , prior(normal(0,1.5), nlpar = 'logShape') +
             prior(normal(0,1.5),nlpar = 'scaleLogisThreshold')
           ,...
  )
  
  try({
    m$fitted_threshold <- fitted(m,nlpar = 'threshold',ndraws = 200)
  },silent=T)
  
  return(m)
  
}
