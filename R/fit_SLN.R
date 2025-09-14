
#' Fit a response time model using a shifted lognormal distribution
#' 
#' Fits a 3-parameter shifted lognormal model using \code{\link[brms]{brm}},
#' with some adjusted defaults in order to make things somewhat easier for the user. 
#'
#' @param rtFormula linear [mixed-effects] model formula defining model to predict RT (i.e., the mean of the lognormal distribution)
#' @param dat Data frame
#' @param sigmaFormula linear [mixed-effects] model formula defining sigma
#' @param ndtFormula linear [mixed-effects] model formula defining non-decision time (NDT). NDT is estimated as a scaled logistic function, which helps with sampling. 
#' @param ndtBound The scalar numeric upper bound of the scaled logistic function defining NDT. Defaults to an upper bound that is the 5th percentile of the RT distribution (i.e., \code{quantile(rt,.05)}). If the model is not initializing successfully, reducing this value may help.
#' @param priorIn Priors to provide to the model. Note the default value which provides a reasonable prior over the possible values of NDT; if providing user-defined priors, it is recommended to also provide this one as well.
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
#' rtFormula <- bSplineFormula(rt ~ (0 + scale(absRat) | subID) , basisVar = d$trialNum, groupingVar = 'subID')
#' 
#' m_SLN <- fit_SLN(rtFormula = rtFormula
#' ,ndtFormula = scaleLogisNDT ~ (1 | subID)
#' ,dat = d
#' , cores = 2
#' )
#' 
fit_SLN <- function(rtFormula
                    ,dat
                    ,sigmaFormula = sigma ~ 1
                    ,ndtFormula = scaleLogisNDT ~ 1
                    ,ndtBound = 'auto'
                    ,priorIn = prior(normal(0,1.5),nlpar = 'scaleLogisNDT')
                    ,...){
  
  require(brms)
  
  # to note in docs:
  #
  # - how empirical priors are calculated
  #
  # - how we augment data for bases
  #
  # - the nature of the `fitted_param`
  #
  # - how to use brm_optimizing() or even better, algorithm='pathfinder'
  
  # # get "empirical" uninformative priors
  # if(is.character(priorIn)){
  #   if(priorIn == 'empirical'){
  #     
  #     rt_log_mean <- signif(mean(log(dat[,rtVar]),na.rm = T),3)
  #     rt_log_sd <- signif(mad(log(dat[,rtVar]),na.rm = T),3)
  #     rt_mean <- signif(exp(rt_log_mean),3) # geometric mean
  #     rt_sd <- signif(mad(dat[,rtVar],na.rm = T),3)
  #     
  #     priorIn <- prior_string(paste0('normal(0,',rt_log_sd,')'),nlpar = 'expoMean') +
  #       prior_string(paste0('normal(',rt_log_mean,',',rt_log_sd,')'),nlpar = 'expoMean', coef = 'Intercept') +
  #       prior_string(paste0('normal(0,',rt_sd,')'),nlpar = 'gaussMean') +
  #       prior_string(paste0('normal(',rt_mean,',',rt_sd,')'),nlpar = 'gaussMean', coef = 'Intercept') 
  #     
  #   }
  # }
  
  # get the non-overlapping bases from the different formulas:
  bases <- data.frame(nullVar = rep(NA,nrow(dat)))
  for(curForm in list(rtFormula, ndtFormula, sigmaFormula)){
    if(!is.null(attr(curForm,'basisDF'))){
      tmpBases <- attr(curForm,'basisDF')
      # colnames(tmpBases) <- gsub('spl4',paste0('spl',ncol(tmpBases)),colnames(tmpBases)) # it might be necessary to try to do something like this, if different basis densities are requested
      for(curCol in colnames(tmpBases)){
        if(is.null(bases[[curCol]])){ 
          bases[,curCol] <- tmpBases[,curCol]
        }
      }
    }
  }
  bases$nullVar <- NULL
  dat <- data.frame(dat,bases)
  
  if(!is.numeric(ndtBound)){
    rtVar <- as.character(rtFormula)[2]
    ndtBound <- signif(quantile(dat[,rtVar],.05,na.rm=T),2)
  }
  ndtLogisForm <- nlf(formula(paste0(
    'ndt ~ ',ndtBound,'*inv_logit(scaleLogisNDT)'
  )))
  
  bForm <- brmsformula(
    rtFormula
    ,sigmaFormula
    ,ndtFormula
  ) + 
    ndtLogisForm
  
  m <- brm(
    bForm
    ,data = dat
    ,family = shifted_lognormal(link_ndt = 'identity')
    ,prior = priorIn
    ,...
  )
  
  try({
    m$fitted_param$LN_mean <- fitted(m, ndraws = 200)
    m$fitted_param$LN_sigma <- fitted(m, dpar = 'sigma', ndraws = 200)
    m$fitted_param$ndt <- fitted(m, dpar = 'ndt', ndraws = 200)
  })
  
  return(m)
  
}
