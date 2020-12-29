#'plotting function.
#'takes two data frames: one "original" dataframe with animal names as rownames (and no extra columns), one imputed with two additional columns:
#'imputation number("sample number") and Animal,
#'The sample number column must be column number 1, but animal column is allowed to vary
#'If there are m multiple imputations and n animals, the first data frame (one_group_orig) should have n rows
#'and the imputed dataframe (one_group_imputed) should have n * m rows. e.g. 6 mice with 5 multiple imputations and d days 
#' means one_group_orig is 6 x (d+1) and one_group_imputed is 35 x (d+2)
#' @param one_group_orig original dataset with missing values with 
#' @param one_group_imputed set of imputed versions of one_group_orig, compiled into a dataset with the first column corresponding to imutation number.
#' @export
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
