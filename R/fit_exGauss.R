#' Fit a response time model using an ex-Gaussian distribution
#' 
#' Fits a 3-parameter ex-Gaussian model using \code{\link[brms]{brm}},
#' with some adjusted defaults in order to make things somewhat easier for the user. 
#' 
#' By default, \code{brms} has the main model formula predict the sum of the exponential and
#' Gaussian components of the ex-Gaussian distribution. This function separates these so they
#' can be estimated separately from one another. This allows simpler partitioning of the model into
#' its constituent components. For instance, it may be theoretically meaningful to have one variable
#' of interest predict the exponential component and a different variable of interest predict the 
#' Gaussian component.
#' 
#' Default priors are defined as follows: The prior for the exponential component intercept,
#' which by \code{brms} default is estimated on a log scale, is defined as
#' a normal distribution with a mean of \code{mean(log(rt))} and a standard deviation of
#' \code{mad(log(rt))}; non-intercept coefficients have a prior of a normal distribution 
#' with a mean of 0 and a standard deviation of \code{mad(log(rt))}.
#' The prior for mean of the Gaussian component intercept is defined as
#' a normal distribution with a mean of \code{exp(mean(log(rt)))} and a standard deviation of
#' \code{mad(rt)}; non-intercept coefficients have a prior of a normal distribution 
#' with a mean of 0 and a standard deviation of \code{mad(rt)}.
#' 
#'
#' @param rtVar Character. The name of the response time variable to be predicted.
#' @param expoFormula linear [mixed-effects] model formula defining the exponential component of the ex-Gaussian distribution.
#' @param dat Data frame
#' @param gaussFormula linear [mixed-effects] model formula defining the mean of the Gaussian component of the ex-Gaussian distribution.
#' @param sigmaFormula linear [mixed-effects] model formula defining the dispersion of the Gaussian component of the ex-Gaussian distribution.
#' @param priorIn Priors to provide to the model. For defaults see \code{Description}.
#' @param ... 
#'
#' @export
#'
#' @examples
#' 
#' d <- TEfits::anstrain
#' 
#' ## generate example response time data
#' d$absRat <- abs(d$ratio)
#' d$rt <- brms::rwiener(nrow(d)
#' ,delta = .2 + d$absRat/5 + scale(d$trialNum)/20
#' , tau = .15
#' , alpha = 1 + sin(d$trialNum/20)/10
#' , beta = .5)[,'q']
#' 
#' expoFormula <- bSplineFormula(expoMean ~ (0 + scale(absRat) | subID) , basisVar = d$trialNum, groupingVar = 'subID')
#' 
#' m_EG <- fit_exGauss('rt'
#' ,expoFormula = expoFormula
#' ,gaussFormula = gaussMean ~ scale(sizeRat) + (1 | subID)
#' ,dat = d
#' , cores = 2
#' )
#' 
fit_exGauss <- function(rtVar
                        ,expoFormula = expoMean ~ 1
                        ,dat
                        ,gaussFormula = gaussMean ~ 1
                        ,sigmaFormula = sigma ~ 1
                        ,priorIn = 'empirical'
                        ,...){
  
  require(brms)
  
  # to note in docs:
  #
  # - how empirical priors are calculated
  #
  # - how we de-convolve relative to brms default gauss+expo
  #
  # - how we augment data for bases
  #
  # - the nature of the `fitted_param`
  #
  # - how to use brm_optimizing() or even better, algorithm='pathfinder'
  
  # get "empirical" uninformative priors
  if(is.character(priorIn)){
    if(priorIn == 'empirical'){
      
      rt_log_mean <- signif(mean(log(dat[,rtVar]),na.rm = T),3)
      rt_log_sd <- signif(mad(log(dat[,rtVar]),na.rm = T),3)
      rt_mean <- signif(exp(rt_log_mean),3) # geometric mean
      rt_sd <- signif(mad(dat[,rtVar],na.rm = T),3)
      
      priorIn <- prior_string(paste0('normal(0,',rt_log_sd,')'),nlpar = 'expoMean') +
        prior_string(paste0('normal(',rt_log_mean,',',rt_log_sd,')'),nlpar = 'expoMean', coef = 'Intercept') +
        prior_string(paste0('normal(0,',rt_sd,')'),nlpar = 'gaussMean') +
        prior_string(paste0('normal(',rt_mean,',',rt_sd,')'),nlpar = 'gaussMean', coef = 'Intercept') 
      
    }
  }
  
  # get the non-overlapping bases from the different formulas:
  bases <- data.frame(nullVar = rep(NA,nrow(d)))
  for(curForm in list(expoFormula, gaussFormula, sigmaFormula)){
    if(!is.null(attr(curForm,'basisDF'))){
      tmpBases <- attr(curForm,'basisDF')
      for(curCol in colnames(tmpBases)){
        if(is.null(bases[[curCol]])){ 
          bases[,curCol] <- tmpBases[,curCol]
        }
      }
    }
  }
  bases$nullVar <- NULL
  dat <- data.frame(dat,bases)
  
  bForm <- brmsformula(
    formula(paste0(rtVar,' ~ gaussMean + expoMean'))
    ,gaussFormula
    ,expoFormula
    ,sigmaFormula
    ,nl = T
  ) + 
    nlf(beta ~ expoMean)
  
  m <- brm(
    bForm
    ,data = dat
    ,family = exgaussian()
    ,prior = priorIn
    ,...
  )
  
  try({
    m$fitted_param$expoMean <- fitted(m, nlpar = 'expoMean' , ndraws = 200,scale = 'resp') 
    m$fitted_param$gaussMean <- fitted(m, nlpar = 'gaussMean' , ndraws = 200)
  })
  
  return(m)
  
  ##-##-##-##
  ##-##-##-##

}
