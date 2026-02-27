keeporder <- function(x){
  x <- as.character(x)
  x <- factor(x, levels=unique(x))
  x
}

# Functions from Rand Wilcox's Rallfun-v40.txt
# https://osf.io/xhe8u/
trimci <- function(x,tr=.2,alpha=.05,null.value=0){
  #  Compute a 1-alpha confidence interval for the trimmed mean
  #  The default amount of trimming is tr=.2
  # x<-elimna(x)
  se<-sqrt(winvar(x,tr))/((1-2*tr)*sqrt(length(x)))
  trimci<-vector(mode="numeric",length=2)
  df<-length(x)-2*floor(tr*length(x))-1
  trimci[1]<-mean(x,tr)-qt(1-alpha/2,df)*se
  trimci[2]<-mean(x,tr)+qt(1-alpha/2,df)*se
  test <- (mean(x,tr)-null.value)/se
  # sig<-2*(1-pt(abs(test),df))
  sig <- 1-pt(abs(test),df)
  list(tmean=mean(x,tr),
       ci=trimci,
       test.stat=test,
       se=se,
       p.value=sig,
       n=length(x),
       df=df)
}

winvar <- function(x,tr=.2){
  #
  #  Compute the gamma Winsorized variance for the data in the vector x.
  #  tr is the amount of Winsorization which defaults to .2.
  #
  remx=x
  # x<-x[!is.na(x)]
  y<-sort(x)
  n<-length(x)
  ibot<-floor(tr*n)+1
  itop<-length(x)-ibot+1
  xbot<-y[ibot]
  xtop<-y[itop]
  y<-ifelse(y<=xbot,xbot,y)
  y<-ifelse(y>=xtop,xtop,y)
  wv<-var(y)
  # if(!na.rm)if(sum(is.na(remx)>0))wv=NA
  wv
}

outbox<-function(x,mbox=FALSE,gval=NA,plotit=FALSE,STAND=FALSE){
  #
  # This function detects outliers using the
  # boxplot rule, but unlike the R function boxplot,
  # the ideal fourths are used to estimate the quartiles.
  #
  # Setting mbox=TRUE results in using the modification
  # of the boxplot rule suggested by Carling (2000).
  #
  x<-x[!is.na(x)] # Remove missing values
  if(plotit)boxplot(x)
  n<-length(x)
  temp<-idealf(x)
  if(mbox){
    if(is.na(gval))gval<-(17.63*n-23.64)/(7.74*n-3.71)
    cl<-median(x)-gval*(temp$qu-temp$ql)
    cu<-median(x)+gval*(temp$qu-temp$ql)
  }
  if(!mbox){
    if(is.na(gval))gval<-1.5
    cl<-temp$ql-gval*(temp$qu-temp$ql)
    cu<-temp$qu+gval*(temp$qu-temp$ql)
  }
  flag<-NA
  outid<-NA
  vec<-c(1:n)
  for(i in 1:n){
    flag[i]<-(x[i]< cl || x[i]> cu)
  }
  if(sum(flag)==0)outid<-NULL
  if(sum(flag)>0)outid<-vec[flag]
  keep<-vec[!flag]
  outval<-x[flag]
  n.out=sum(length(outid))
  list(out.val=outval,out.id=outid,keep=keep,n=n,n.out=n.out,cl=cl,cu=cu)
}

idealf<-function(x,na.rm=FALSE){
  #
  # Compute the ideal fourths for data in x
  #
  if(na.rm)x<-x[!is.na(x)]
  j<-floor(length(x)/4 + 5/12)
  y<-sort(x)
  g<-(length(x)/4)-j+(5/12)
  ql<-(1-g)*y[j]+g*y[j+1]
  k<-length(x)-j+1
  qu<-(1-g)*y[k]+g*y[k-1]
  list(ql=ql,qu=qu)
}


