rm.all.na=function(df){
df = as.data.frame(Filter(function(x)!all(is.na(x)), df))
df = df[rowSums(is.na(df)) != ncol(df), ]
}

substrRight <- function(x, n){
  substr(x, nchar(x)-n+1, nchar(x))
}

rho.func = function(z){
  num =  exp(2*z)-1
  den = exp(2*z)+1
  num/den
}
covar.one = function(u,t,theta6,theta7,theta8){
  exp(2*theta7) + (rho.func(theta8)^(abs(t-u)))*(exp(2*theta6)/(1-(rho.func(theta8)^2)))
}

mu.cond = Vectorize(function(tt, theta){
  if(theta[4]==0){return(theta[2]+ (theta[1]-theta[2])*exp(-exp(theta[3])*tt) )}else{
    a = 1+(exp(theta[4]*(theta[2]-theta[1])) -1)
    theta[2] - 1/theta[4]*(log(a)-exp(theta[2]*theta[4]-theta[3])*tt)
  }
  
},vectorize.args = "tt")

slice = Vectorize(function(j, vcv.emp, mu){
  animal = unique(zerogy.long$Animal)[j]
  xj = rep(NA,length(days))
  names(xj)=days
  times.xj = zerogy.long$time[zerogy.long$Animal==animal]
  xj[as.character(times.xj)] = zerogy.long$LTV[zerogy.long$Animal==animal]
  vec = as.matrix(xj-mu)
  vec[is.na(vec)]=0
  t(vec)%*%solve(vcv.emp)%*%(vec)
}, vectorize.args = "j")

compute.neg2lnL = function(theta, tt, theta4){
  theta0=c(theta[1:3], theta4)
  theta1=theta[4:6]
  if(!(theta[2]>=theta[1])){return(Inf)}else{
    vcv.emp = outer(tt,tt,covar.one, theta6=theta1[1],theta7=theta1[2],theta8=theta1[3])
    mu = mu.cond(tt, theta=theta0)
    n*length(tt)*log(2*pi)+n*log(det(vcv.emp)) + sum(slice(1:n, vcv.emp,mu),na.rm=T)
    }
}

compute.lnL.for.gradient = function(theta0, t, theta4, theta1,j){
  theta0=c(theta0, theta4)
  if(!(theta[2]>=theta[1])){return(Inf)}else{
    vcv.emp = outer(t,t,covar.one, theta6=theta1[1],theta7=theta1[2],theta8=theta1[3])
    mu = mu.cond(t, theta=theta0)
    -0.5*(n*length(t)*log(2*pi)+n*log(det(vcv.emp)) + sum(slice(j, vcv.emp,mu),na.rm=T))
  }
}