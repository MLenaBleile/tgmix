prob.drop = Vectorize(function(data.vec, beta, k){
  tt= length(data.vec)
  unlogit(beta[1]*k+beta[2]*data.vec[k] + beta[3]*data.vec[k-1])
}, vectorize.args = "k")

unlogit = function(x){return(exp(x)/(1+exp(x)))}

dropout.process = function(data.vec,beta){
    drop.probs = prob.drop(data.vec, beta, 4:length(data.vec))
    #print(drop.probs)
    drop = 0
    time=5
    while(drop==0&time<=length(data.vec)){
      drop = Rlab::rbern(1, prob=drop.probs[time-4])
      print(c(drop.probs[time],time))
      time=time+1
    }
    if(time<=length(data.vec)){data.vec[(time-1):length(data.vec)]=NA}
    data.vec
}

one.dataset=mouse.sim[1,,]
dd=t(apply(one.dataset,1, dropout.process, beta=c(-0.1,-0.1, -0.01)))
dd=as.data.frame(dd)
dd$Animal=1:nrow(dd)
colnames(dd)[-c(1)] = paste("Day", time, sep="")
dd.long = reshape(data=dd, 
                  direction="long",
                  varying = paste("Day",time,sep=""),
                  v.names="LTV",
                  idvar= "Animal",
                  times = seq(0,63, by=7))
write.csv(dd.long, file="testdat.csv", na=".")
samples=read.csv("samples_extra.csv", na.strings = ".")

df.imputed =impute(1,dd,samples)
df.imputed$imputation=1

for(imp in unique(samples$sample)[-1]){
  new= impute(imp,dd,samples)
  new$imputation=imp
  df.imputed = rbind(df.imputed, new)
}

df.imputed= df.imputed%>%dplyr::select(imputation,Animal, everything())
dd = dd%>%dplyr::select(Animal, everything())
make_betterplots(dd, df.imputed)

diffs.real = matrix(NA, nrow=nrow(dd), ncol=10)
diffs.imputed = matrix(NA, nrow=nrow(dd), ncol=10)
for(i in 1:nrow(diffs.real)){
  diffs.real[i,] = as.numeric(dd[i,-c(1)]) - mouse.baseline
  diffs.imputed[i,] = as.numeric(df.imputed[i,-c(1:2)])-mouse.baseline
}

one_group_orig=filter_all_na(dd[1:8,-c(1)])
one_group_imputed=df.imputed[df.imputed$Animal%in%1:8,1:(ncol(one_group_orig)+2)]
make_betterplots(one_group_orig, one_group_imputed, day=time[1:7])
lines(time[1:7], mouse.baseline[1:7])
