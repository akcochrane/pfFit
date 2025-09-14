#' Fit a logistic psychometric function generalized mixed-effects model
#' 
#' Fit a model to binary response data, where the link between positively-bounded \code{x}
#' and the binary \code{y} is a logistic function. The logistic function is defined
#' with \code{threshold} and \code{bias} parameters, which correspond broadly to
#' a "slope" and "location" respectively.
#' 
#' This function serves to set up a nonlinear \code{\link[brms]{brm}} model with a 
#' bernoulli response distribution.
#'
#' @param y The binary variable to be predicted. It should change as a monotonic function of \code{x}.
#' @param x Numeric variable used as the "x-axis" of the psychometric function
#' @param thresholdForm A formula, either a linear model or linear mixed-effects model, specifying the predicted values of \code{scaleLogisThreshold}. This is the scaled inverse-logit-transformed threshold, which provides more stable estimation than sampling the distributions of the raw threshold.
#' @param biasForm A formula, either a linear model or linear mixed-effects model, specifying the predicted values of \code{bias}.
#' @param dat Data frame with variables with which to fit the model.
#' @param lapseRate The value the logistic function takes at \code{x} values of \code{-Inf}, with the values approaching \code{1-lapseRate} as \code{x} approaches \code{Inf}.
#' @param range_x Two numbers (e.g., \code{c(.1,5)}). Defines the plausible range of \code{threshold} values, which act to limit the possible values of the threshold. Defaults to the \code{min} and \code{max} of \code{abs(x)}, but if the user wants their threshold to have different boundaries, they can be user-specified.
#' @param ... 
#'
#' @export
#'
#' @examples
#' 
#'   d <- TEfits::anstrain
#'   d$absRat <- abs(d$ratio)
#'   
#'   threshSpline <- bSplineFormula(scaleLogisThreshold ~ trialNum
#'   ,d$trialNum,groupingVar = 'subID')
#'   
#'   m_logistic <- fit_logistic(y='resp',x='ratio'
#'   ,thresholdForm = threshSpline
#'   ,bias = bias ~ (1 | subID)
#'   ,dat = d
#'    )
#' 
fit_logistic <- function(y
                         , x
                         , thresholdForm = scaleLogisThreshold ~ 1
                         , biasForm = bias ~ 1
                         , dat
                         , lapseRate = .01
                         , range_x = 'auto'
                         , ...){
  require(brms)
  
  if(!is.null(attr(thresholdForm,'basisDF'))){
    dat <- data.frame(dat,attr(thresholdForm,'basisDF'))
  }
  
  if(length(range_x)==1){
    range_x <- signif(range(abs(d[,x]),na.rm = T),2)
  }
  
  threshForm_scale <- formula(paste0('threshold ~ ',range_x[1],' + (',range_x[2],'-',range_x[1],')*inv_logit(scaleLogisThreshold)'))
  
  bForm <- brmsformula(
    pfFormula(y,pf = 'logistic'
              ,link_x = x
    )$pf
    ,thresholdForm
    ,biasForm
    ,nl = T
  ) +
    nlf(threshForm_scale) + 
    nlf(formula(paste0('lapseRate ~ ',lapseRate)))
  
  m <- brm(bForm
           ,data = dat
           ,data2 = dat
           ,family = bernoulli(link = 'identity')
           , prior_string(paste0('normal(0,',signif(sd(d[,x],na.rm=T),2),')'), nlpar = 'bias') +
             prior(normal(0,1.5),nlpar = 'scaleLogisThreshold')
           ,...
  )
  
  try({
    m$fitted_threshold <- fitted(m,nlpar = 'threshold',ndraws = 200)
  },silent=T)
  
  return(m)
  
}
