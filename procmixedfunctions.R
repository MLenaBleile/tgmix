#############################################
##procmixedfunctions.R                     ##
##Tumor Growth Data, Multiple Imputation   ##
##MaryLena Bleile, 12 December 2020        ## 
#############################################

###utility functions
substrRight <- function(x, n){
  substr(x, nchar(x)-n+1, nchar(x))
}
rm.na=function(x){if(is.na(mean(x))){return(x[-which(is.na(x))])}else{return(x)}}
filter_all_na = function(data.orig){
  data.orig = as.data.frame(Filter(function(x)!all(is.na(x)), data.orig))
  data.orig = data.orig[rowSums(is.na(data.orig)) != ncol(data.orig), ]
  data.orig
}


#'tgmix_imputation
#'takes a data frame with a categorical columns corresponding to
#'animal id and optionally a column corresponding to treatment group. can specify animal column with option animalcol
#'treatment group col specified by tgcol
#'optional vector of time points "day"; otherwise timepoints are assumed to be 
#' 1 through length(tumor observations)
#'produces filled-in dataset and optionally plots results (plots=T)
#' can specify number of samples (nsamp)
#' works for a single treatment group. need to apply by group for complete dataset.
#' 
#'
#' This function loads a file as a matrix. It assumes that the first column
#' contains the rownames and the subsequent columns are the sample identifiers.
#' Any rows with duplicated row names will be dropped with the first one being
#' kepted.
#'
#' @param one.group data frame with a categorical columns corresponding to animal id and optionally a column corresponding to treatment group.
#' @return options = nsamp versions of the impured dataset with optional plots
#' @export

tgmix_imputation = function(one.group, day=NA, plots=T,animalcol=1, nsamp=20, chains=4, cores=4,iter=2000,warmup=500){
  kcat=1
  
 
  
  ##single imputation function
  ## Function to perform a single round of sampled imputations for a dataset (df)
  ## impnum:imputation iteration number
  ##df : dataframe with observed/missing values, 
  ## samples: dataframe with least squares means(lsm) and covariance parameters (covp) from mixed model 
  #tgcol: treatment group column
  impute = function(impnum, data,samples,animalcol=1){
    for(a in data$Animal){
      #get index of desired animal
      aidx = which(data$Animal==a)
      #find which values are missing; obtain their indices
      naidx = which(is.na(data[which(data$Animal==a),-c(animalcol)]))
      #get a dataframe with just the least squares means
      midx=which(substr(colnames(samples), 1,3)=="lsm")
      meanss = samples[,midx]
      if(length(naidx>0)){
        #create diagonal and off-diagonal matrix elements of the partitioned Sigma matrix (corresponding to known vs unknown obs)
        sig.offdiag = matrix(data=rep(samples$covp1[impnum],length(naidx)*(ncol(meanss)-length(naidx))), ncol=length(naidx))
        p1 = ncol(meanss) - length(naidx)
        p2=length(naidx)
        sig.11 = matrix(data=rep(samples$covp1[impnum], p1*p1), ncol=p1)
        diag(sig.11) = rep(samples$covp2[impnum]+samples$covp1[impnum], p1)
        sig.22 = matrix(data=rep(samples$covp1[impnum], p2*p2), ncol=p2)
        diag(sig.22) = rep(samples$covp2[impnum]+samples$covp1[impnum], p2)
        #get the desired parameters for the distribtution
        muvec = meanss[impnum,naidx] +t(sig.offdiag)%*%solve(sig.11)%*%(as.matrix(as.numeric(meanss[impnum, -naidx])-as.numeric(data[aidx,-c(1, 1+naidx)])))
        
        sigmat = sig.22 - t(sig.offdiag)%*%solve(sig.11)%*%sig.offdiag
        #generate obs from this distribution and impute
        data[aidx,(naidx+1)] = MASS::mvrnorm(n=1, mu=as.numeric(muvec), Sigma=sigmat)
        
      }
    }
    data
  }
  
  ###Function to generate samples using a mixed model
  get_samples = function(one.group,days,animalcol=1){
    names(one.group)[animalcol]="Animal"
    data.long = reshape(one.group,
                        direction="long",
                        varying = list(names(one.group)[-c(animalcol)]),
                        v.names="LTV",
                        idvar= "Animal",
                        times = days)
    tgmix = data.long
    
    ##fit using time as a categorical like before
    tgmix_fit = brms::brm(data = tgmix,
                    family = gaussian,
                    formula = LTV ~ as.factor(time) + (1|Animal),
                    iter = iter, warmup = warmup, chains = chains, cores = cores,
                    control = list(adapt_delta = .975, max_treedepth = 20),
                    seed = 190831)
    samples_init = brms::posterior_samples(tgmix_fit, subset = sample(warmup:iter, nsamp))
    get.mus1 = function(colnum){samples_init$b_Intercept+ samples_init[,colnum]}
    post.mus = sapply(2:length(days), get.mus1)
    post.mus = cbind(as.matrix(samples_init[,1]), post.mus)
    colnames(post.mus) = paste("lsm",1:length(days),sep="")
    
    covpTime = (samples_init$sigma)^2
    
    
    covpAnimal = (samples_init$sd_Animal__Intercept)^2
    
    samples = cbind(as.matrix(covpTime),as.matrix(covpAnimal),post.mus)
    colnames(samples)[1:2] = c("covp2","covp1")
    
    samples = as.data.frame(samples)
    mono.inc.idx= apply(samples[,-c(1,2)], 1, function(x)all(x==cummax(x)))
    samples = samples[all()]
    samples$sample = 1:nsamp
    samples
  }
  
  #rename names to be standardized
  names(one.group)[animalcol]="Animal"
  
  ##get time indices
  if(is.na(day)){day = 1:(ncol(one.group)-kcat)}
  ##generate posterior samples
  samples = get_samples(one.group,days=day)
  #use function to do all of the imputation numbers
  df.imputed =impute(1,one.group,samples)
  df.imputed$imputation=1
  
  for(imp in unique(samples$sample)[-1]){
    new= impute(imp,one.group,samples)
    new$imputation=imp
    df.imputed = rbind(df.imputed, new)
  }
  cols = ncol(df.imputed)
  
  df.imputed= df.imputed[,c("imputation",names(df.imputed)[-c(cols)] )]
                        
  rownames(one.group)=one.group$Animal
  if(plots==T){
    
    suppressWarnings(make_betterplots(one.group[,-c(animalcol)],df.imputed, day=day, animalcol=which(names(df.imputed)=="Animal")))
  }
  df.imputed
}


#'plotting function.
#'takes two data frames: one "original" dataframe with animal names as rownames (and no extra columns), one imputed with two additional columns:
#'imputation number("sample number") and Animal,
#'The sample number column must be column number 1, but animal column is allowed to vary
#'If there are m multiple imputations and n animals, the first data frame (one_group_orig) should have n rows
#'and the imputed dataframe (one_group_imputed) should have n * m rows. e.g. 6 mice with 5 multiple imputations and d days 
#' means one_group_orig is 6 x (d+1) and one_group_imputed is 35 x (d+2)
#' @param one_group_orig original dataset with missing values with 
#' @param one_group_imputed set of imputed versions of one_group_orig, compiled into a dataset with the first column corresponding to imutation number.
make_betterplots = function(one_group_orig,one_group_imputed,day=NA, animalcol=2,plotname="Counterfactual Imputation Plot"){
  kcat=length(c(animalcol))
  if(suppressWarnings(is.na(day))){day=1:(ncol(one_group_imputed)-kcat-1)}
  day.numeric=as.numeric(day)
  cols = rainbow(nrow(one_group_orig))
  plot(day.numeric, one_group_orig[1,], pch=16,main=plotname,xlim=c(min(day.numeric), max(day.numeric)+6),
       ylim=c(min(one_group_orig, na.rm=T),max(apply(one_group_imputed[,-c(1,(animalcol))],2,as.numeric),na.rm = T)), xlab="time", ylab="LTV", type="n")
  rect.corners = numeric()
  for(i in 1:nrow(one_group_orig)){
    one_mouse_imputed=one_group_imputed[as.character(one_group_imputed[,(animalcol)])==rownames(one_group_orig)[i],-c(1, animalcol)]
    
    if(!(length(rm.na(as.numeric(one_mouse_imputed[1,])))==length(rm.na(as.numeric(one_group_orig[i,]))))){
      
      # print(day.numeric)
      points(day.numeric, one_group_orig[i,],pch=16)
      lines(day.numeric, one_group_orig[i,],col=cols[i])
      
      mouse.known=rm.na(one_group_orig[i,])
      mouseinit = one_group_orig[i,]
      uknownidx=which(is.na(mouseinit))
      mouse.unknown = as.matrix(one_mouse_imputed[,uknownidx[uknownidx<=ncol(one_mouse_imputed)]])
      day.unknown=day.numeric[uknownidx]
      med=apply(mouse.unknown, 2, median, na.rm=T)
      iqr = apply(mouse.unknown,2, IQR, na.rm=T)
      avg = apply(mouse.unknown,2,mean,na.rm=T)
      sd = apply(mouse.unknown,2,sd,na.rm=T) 
      imputeidx = which(!is.na(avg))
      for(j in 1:length(avg)){
        x.prop = day.unknown[j]
        while (x.prop %in% rect.corners) {
          x.prop=x.prop+.8
        }
        arrows(x.prop,y0=avg[j]-3*sd[j],y1=avg[j]+3*sd[j],angle=90, code=3, lty=2)
        rect(x.prop-0.4,med[j]-0.5*iqr[j],x.prop+0.4,med[j]+0.5*iqr[j], col=cols[i])
        rect.corners = c(rect.corners,x.prop)
        points(x.prop, avg[j], pch=8)
        
      }
      lines(day.numeric,c(mouse.known, avg), col=cols[i])
    }
  }
}  




