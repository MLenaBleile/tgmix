rho.func = function(z){
  num =  exp(2*z)-1
  den = exp(2*z)+1
  num/den
}



covar.one = function(u,v,theta6,theta7,theta8){
  exp(2*theta7) + (rho.func(theta8)^(abs(u-v)))*  ( exp(2*theta6) / (1-rho.func(theta8)^2))
}

theta0 = c(5,11,-2,0)
theta6 = -1.2
theta7 = -2
theta8 = .7
theta=c(theta0, theta6, theta7, theta8)

time=seq(0,63,by=7)

#variance-covariance matrix for a single mouse
vcv.mat = outer(time,time,covar.one,theta6=theta6, theta7=theta7, theta8=theta8)

mu.cond = Vectorize(function(tt, theta, log.theta3=T){
  if(log.theta3==T){
  if(theta[4]==0){return(theta[2]+ (theta[1]-theta[2])*exp(-exp(theta[3])*tt) )}else{
    a = 1+(exp(theta[4]*(theta[2]-theta[1])) -1)
    theta[2] - 1/theta[4]*(log(a)-exp(theta[2]*theta[4]-theta[3])*tt)
  }}
  else{
    if(theta[4]==0){return(theta[2]+ (theta[1]-theta[2])*exp(-theta[3]*tt) )}else{
      a = 1+(exp(theta[4]*(theta[2]-theta[1])) -1)
      theta[2] - 1/theta[4]*(log(a)-exp(theta[2]*theta[4])*(-theta[3])*tt)
    }
  }
  
},vectorize.args = "tt")

v.pg = Vectorize(function(tt, theta, dt){
  t1=dt[1]
  t2=dt[2]
  if(tt<=t1){return(mu.cond(tt, theta[1:4]))}
  else if (tt<=t2){
    thet.t1 <- c(mu.cond(tt=t1, theta[1:4]), theta[2], exp(theta[3])*theta[5], theta[4])
    return(mu.cond(tt=tt-t1, thet.t1, log.theta3=F))
  }
  else{
    thet.t1 <- c(mu.cond(tt=t1, theta[1:4]), theta[2], exp(theta[3])*theta[5], theta[4])
    thet.t2 = c(mu.cond(tt=t2-t1, thet.t1, log.theta3 = F), theta[2], exp(theta[3])*theta[5]*theta[6], theta[4])
    return(mu.cond(tt=tt-t2, theta = thet.t2, log.theta3 = F))
  }
}, vectorize.args = "tt")

mouse.baseline.pg=numeric(10)
for(v in 1:length(time)){
  mouse.baseline.pg = v.pg(time, theta1, dt=c(20,30))
}


M=300
n=114

mouse.baseline = mu.cond(time, theta=c(5,11,-2,0,-1.2,-2,0))

mouse.sim = array(data=NA, dim=c(M,n, length(time)))

for(m in 1:M){
  for(i in 1:n){
  mouse.sim[m,i,] = mouse.baseline+MASS::mvrnorm(n=1,rep(0, length(time)),Sigma=vcv.mat)
  }
}

simulation.slice = Vectorize(function(j,m, mu, vcv.emp){
  t(as.matrix(mouse.sim[m,j,] - mu))%*%solve(vcv.emp)%*%(as.matrix(mouse.sim[m,j,]-mu))
}, vectorize.args = "j")

compute.2lnL = Vectorize(function(theta, m=1){
  if(!(theta[2]>=theta[1])){return(Inf)}else{
  vcv.emp = outer(time,time,covar.one, theta6=theta[6-1],theta7=theta[7-1],theta8=theta[8-1])
  mu = mu.cond(time, theta=theta[1:4])
 
  return(n*length(time)*log(2*pi)+n*log(det(vcv.emp)) + sum(simulation.slice(1:n,m,mu, vcv.emp)))
  }
},vectorize.args = "m")

starting.params = c(3, 10,-1.5,0.5,-1,-3,.5)

max.L = function(m){optim(par=starting.params,compute.2lnL, method="Nelder-Mead", m=m )}

params = matrix(NA, nrow=M,ncol=length(theta))
colnames(params) = paste("theta",c(1,2,3,4,6,7,8))



lnL.for.gradient.computation = function(thetaa,thetavar,unit, m=1){
  vcv.emp = outer(time,time,covar.one, theta6=thetavar[1],theta7=thetavar[2],theta8=thetavar[3])
  mu = mu.cond(time, theta=thetaa)
  -0.5*(length(time)*log(2*pi)+log(det(vcv.emp)) + t(as.matrix(mouse.sim[m,unit,] - mu))%*%solve(vcv.emp)%*%(as.matrix(mouse.sim[m,unit,]-mu)))
}


M=10
vcv.avg = array(NA, dim=c(M,4,4))
gradients=array(NA, dim=c(n,4,4))

m=1
for(m in 1:M){
    for(unit in 1:n){
    p = max.L(m)
    params[m,] = p$par
    stopifnot(p$convergence==0)
    #p
    g = rootSolve::gradient(lnL.for.gradient.computation, x=params[m,1:4],thetavar=params[m,5:7],unit=unit, m=m, pert=.01)
    gradients[unit,,] = as.matrix(t(g))%*%as.matrix(g)
    }
    vcv.avg[m,,] = 1/(apply(gradients,2:3, sum))
  }
  
params[1:10,]
