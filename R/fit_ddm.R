
#' Fit drift-diffusion model
#' 
#' Fits a 4-parameter drift-diffusion model using \code{\link[brms]{brm}},
#' with some adjusted defaults in order to make things somewhat easier for the user.
#' 
#' The primary change from the default \code{brms} implementation is estimating the NDT 
#' within a scaled logistic function. This provides natural boundaries, improves the 
#' chances that model initialization errors will not be encountered, and in preliminary tests
#' has helped restrict parameter trade-offs with other potential parameters of interest (e.g., drift rate).
#' If initialization errors are still encountered, it's recommended to decrease 
#' \code{init_r} (see \code{\link[brms]{brm}}).
#'
#' @param drFormula linear [mixed-effects] model formula defining drift rate
#' @param dat data frame
#' @param bsFormula linear [mixed-effects] model formula defining boundary separation
#' @param ndtFormula linear [mixed-effects] model formula defining non-decision time (NDT). NDT is estimated as a scaled logistic function, which helps with sampling. 
#' @param biasFormula linear [mixed-effects] model formula defining bias
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
#' d$drSpl <- bSplineMat(basisVar = d$trialNum, nbasis = 7)
#' d$bsSpl <- bSplineMat(basisVar = d$trialNum)
#' 
#' m_ddm <- fit_ddm(drFormula = rt | dec(acc) ~ scale(trialNum) +  (drSpl + scale(absRat) || subID)
#'  ,bsFormula = bs ~ (bsSpl || subID)
#'  ,dat = d
#'  , cores = 2 , chains = 2 # only run 2 chains for efficiency
#'  )
#' 
fit_ddm <- function(drFormula
                    ,dat
                    ,bsFormula = bs ~ 1
                    ,ndtFormula = scaleLogisNDT ~ 1
                    ,biasFormula = bias ~ 1
                    ,ndtBound = 'auto'
                    ,priorIn = prior(normal(0,1.5),nlpar = 'scaleLogisNDT')
                    ,...){
  
  require(brms)
  
  # get the non-overlapping bases from the different formulas:
  bases <- data.frame(nullVar = rep(NA,nrow(d)))
  for(curForm in list(drFormula,bsFormula,ndtFormula,biasFormula)){
    if(!is.null(attr(curForm,'basisDF'))){
      tmpBases <- attr(curForm,'basisDF') # colnames(tmpBases) <- gsub('spl4',paste0('spl',ncol(tmpBases)),colnames(tmpBases)) # it might be necessary to try to do something like this, if different basis densities are requested
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
    rtVar <- gsub(' ','',strsplit(as.character(drFormula)[2],split = '|',fixed = T)[[1]][1])
    ndtBound <- signif(quantile(dat[,rtVar],.05,na.rm=T),2)
  }
  
  ndtLogisForm <- nlf(formula(paste0(
    'ndt ~ ',ndtBound,'*inv_logit(scaleLogisNDT)'
  )))
  
  m <- brm(
    brmsformula(
      drFormula
      ,bsFormula
      ,biasFormula
      ,ndtFormula
    ) + ndtLogisForm
    ,data = dat
    ,family = wiener(link_ndt = 'identity')
    ,prior = priorIn
    ,...
  )
  
  try({
    m$fitted_param$dr <- fitted(m, scale = 'linear' , ndraws = 200)
    for(curParam in c('bs','ndt','bias')){
      m$fitted_param[[curParam]] <- fitted(m, dpar = curParam, ndraws = 200) # because `scale='response'` is used here we're returning raw bs and bias; if `scale='linear'` was used we would be returning log-bs and logit-bias
    }
  },silent = T)
  
  return(m)
  ##-##-##-##
  ##-##-##-##
  ##-##-##-##
  ##-##-##-##
  
  
}

