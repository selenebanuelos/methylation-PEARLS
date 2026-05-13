### Author: Alan E. Hubbard
### Note: comments below line 11 were added by Selene Banuelos

gee.post.estimate=function(glmob,comps,labs=NULL,rounded=3,exponentiate=FALSE) {
  # glmob is the gee object 
  # comps is the vector (or matrix if more than one) of constants to multiply by 
  ###  coefficients for the desired parameter estimate.
  # labs are labels (character vector) that you wish to call each comparison
  # rounded is the number of digits beyond the decimal.
  # exponentiate exponentiates the resulting estimate (e.g., if you used 
  ###  logistic regression).
  if(is.matrix(comps)==FALSE) {
    comps=t(as.matrix(comps))}

  # get robust variance from gee object
  vc=glmob$robust.variance
  
  # get estimates of coefficients
  ests = coef(glmob)
  
  # multiply estimates by constants given in comps argument
  linear.ests=as.vector(comps%*%ests)
  
  vcests=comps%*%vc%*%t(comps)
  ses=sqrt(diag(vcests))
  pvalue=(2*(1-pnorm(abs(linear.ests/ses))))
  if(exponentiate) {
    l95ci = exp(linear.ests - 1.96 * ses)
    exp_beta = exp(linear.ests)
    u95ci = exp(linear.ests + 1.96 * ses)
    summ = cbind(Ratio.est=format(round(exp_beta,rounded),nsmall=rounded),CI=paste(format(round(l95ci,rounded),nsmall=rounded),format(round(u95ci,rounded),nsmall=rounded),sep=", "), pvalue=format(round(pvalue,rounded),nsmall=rounded))     
    rownames(summ)=labs
  }
  if(exponentiate==FALSE) {
    l95ci = (linear.ests - 1.96 * ses)
    beta = (linear.ests)
    u95ci = (linear.ests + 1.96 * ses)
    summ = cbind(Est=format(round(beta,rounded),nsmall=rounded),CI=paste(format(round(l95ci,rounded),nsmall=rounded),format(round(u95ci,rounded),nsmall=rounded),sep=", "), pvalue=format(round(pvalue,rounded+1),nsmall=rounded))  
    rownames(summ)=labs
  }
  return(summ) }
