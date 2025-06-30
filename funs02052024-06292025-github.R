Kern.FUN <- function(zz,zi,bw) 
{ 
  out = (VTM(zz,length(zi))- zi)/bw
  dnorm(out)/bw
  
}


VTM<-function(vc, dm){
  matrix(vc, ncol=length(vc), nrow=dm, byrow=T)
}


resam<- function(vv,t,t.0,tt,data,data1,data2,indexindex){
  n1=nrow(data1)
  n2=nrow(data2)
  # stept=0.05
  # tt=seq(t.0,5,stept)
  causal=rep(NA,length(tt)); causal2=rep(NA,length(tt))
  causals=rep(NA,length(tt)); causals2=rep(NA,length(tt))
  causalind=rep(NA,length(tt)); causalind2=rep(NA,length(tt))
  causalsind=rep(NA,length(tt)); causalsind2=rep(NA,length(tt))
  for (j in 1:length(tt)){
  t=tt[j] 
  ################ pte2 given data1  
  xob=data1[,1];deltaob=data1[,2];aob=data1[,3];sob=data1[,4];n=n1;v=vv[indexindex]
  
  from = min(sob[sob!=0],na.rm = T); to = quantile(sob[sob!=0],.95,na.rm = T); step=((to - from)/nn)
  s=seq(from, to, by = step)
  
  #### optimal function 
  {
  temp=WEIGHT.p(xob,deltaob,aob,n=n,v)
  wei.t0=temp[,1];wei.t=temp[,2]
  
  bw = 1.06*sd(sob,na.rm=T)*n^(-1/5)/(n^0.06)
  kern = Kern.FUN(zz=s,zi=sob,bw)
  
  f0.s.t0.t0.hat=apply(as.numeric(v)*(aob==0)*(xob>t.0)*(sob<=t.0)*kern*wei.t0,2,sum,na.rm=T)/sum(as.numeric(v)*(aob==0)*wei.t0)
  f0.s.t.t0.hat=apply(as.numeric(v)*(aob==0)*(xob>t)*(sob<=t.0)*kern*wei.t,2,sum,na.rm=T)/sum(as.numeric(v)*(aob==0)*wei.t)
  f1.s.t0.t0.hat=apply(as.numeric(v)*(aob==1)*(xob>t.0)*(sob<=t.0)*kern*wei.t0,2,sum,na.rm=T)/sum(as.numeric(v)*(aob==1)*wei.t0)
  f1.s.t.t0.hat=apply(as.numeric(v)*(aob==1)*(xob>t)*(sob<=t.0)*kern*wei.t,2,sum,na.rm=T)/sum(as.numeric(v)*(aob==1)*wei.t)
  temp=(sob>t.0); temp[is.na(temp)]=1
  p0.t0.t0.hat=sum(as.numeric(v)*(aob==0)*(xob>t.0)*temp*wei.t0, na.rm = T)/sum(as.numeric(v)*(aob==0)*wei.t0)
  p0.t.t0.hat=sum(as.numeric(v)*(aob==0)*(xob>t)*temp*wei.t, na.rm = T)/sum(as.numeric(v)*(aob==0)*wei.t)
  p1.t0.t0.hat=sum(as.numeric(v)*(aob==1)*(xob>t.0)*temp*wei.t0, na.rm = T)/sum(as.numeric(v)*(aob==1)*wei.t0)
  p1.t.t0.hat=sum(as.numeric(v)*(aob==1)*(xob>t)*temp*wei.t, na.rm = T)/sum(as.numeric(v)*(aob==1)*wei.t)
  
  integrand<-f0.s.t0.t0.hat^2/f1.s.t0.t0.hat
  temp=(integrand[1] + integrand[nn+1] + 2*sum(integrand[seq(2,nn,by=2)]) + 4 *sum(integrand[seq(3,nn-1, by=2)]) )*step/3
  de=temp+p0.t0.t0.hat^2/p1.t0.t0.hat
  
  mu0t=sum(as.numeric(v)*(aob==0)*(xob>t)*wei.t)/sum(as.numeric(v)*(aob==0)*wei.t)
  mu1t=sum(as.numeric(v)*(aob==1)*(xob>t)*wei.t)/sum(as.numeric(v)*(aob==1)*wei.t)
  
  integrand<-f0.s.t0.t0.hat*f1.s.t.t0.hat/f1.s.t0.t0.hat
  temp=(integrand[1] + integrand[nn+1] + 2*sum(integrand[seq(2,nn,by=2)]) + 4 *sum(integrand[seq(3,nn-1, by=2)]) )*step/3
  nu=mu0t-temp-p0.t0.t0.hat*p1.t.t0.hat/p1.t0.t0.hat
  
  lambda=nu/de
  
  gs.hat=(lambda*f0.s.t0.t0.hat+f1.s.t.t0.hat)/f1.s.t0.t0.hat
  g3.hat=(lambda*p0.t0.t0.hat+p1.t.t0.hat)/p1.t0.t0.hat
  }
  
  #### pte
  xob=data2[,1];deltaob=data2[,2];aob=data2[,3];sob=data2[,4];n=n2;v=vv[-indexindex]
  temp=WEIGHT.p(xob,deltaob,aob,n=n,v)
  wei.t0=temp[,1];wei.t=temp[,2]
  
  # causal=mu1t-mu0t
  causal.t0=sum(as.numeric(v)*(xob>t.0)*aob*wei.t0)/sum(as.numeric(v)*aob*wei.t0)-sum(as.numeric(v)*(xob>t.0)*(1-aob)*wei.t0)/sum(as.numeric(v)*(1-aob)*wei.t0)
  causal[j]=sum(as.numeric(v)*(xob>t)*aob*wei.t)/sum(as.numeric(v)*aob*wei.t)-sum(as.numeric(v)*(xob>t)*(1-aob)*wei.t)/sum(as.numeric(v)*(1-aob)*wei.t)
  tempind=c(sapply(1:n, function(kk){which.min(abs(sob[kk]-s))})); tempind=as.numeric(tempind)
  temp=(sob>t.0); temp[is.na(temp)]=1
  causals[j]=sum(as.numeric(v)*aob*wei.t0*(xob>t.0)*(sob<=t.0)*gs.hat[tempind],na.rm = T)/sum(as.numeric(v)*aob*wei.t0)+
    sum(as.numeric(v)*(xob>t.0)*temp*g3.hat*aob*wei.t0)/sum(as.numeric(v)*aob*wei.t0) -
    (sum(as.numeric(v)*(1-aob)*wei.t0*(xob>t.0)*(sob<=t.0)*gs.hat[tempind],na.rm = T)/sum(as.numeric(v)*(1-aob)*wei.t0)+
       sum(as.numeric(v)*(xob>t.0)*temp*g3.hat*(1-aob)*wei.t0)/sum(as.numeric(v)*(1-aob)*wei.t0))
  pte.es=causals[j]/causal[j]
  
  # causalind[j]=sum(as.numeric(v)*(xob>t)*aob*wei.t)/sum(as.numeric(v)*aob*wei.t)-sum(as.numeric(v)*(xob>t)*(1-aob)*wei.t)/sum(as.numeric(v)*(1-aob)*wei.t)
  # temp=(xob>t.0); temp[is.na(temp)]=1
  # causalsind[j]=sum(as.numeric(v)*(xob>t.0)*temp*g.hat*aob*wei.t0)/sum(as.numeric(v)*aob*wei.t0) -
  #   sum(as.numeric(v)*(xob>t.0)*temp*g.hat*(1-aob)*wei.t0)/sum(as.numeric(v)*(1-aob)*wei.t0)
  # pteind.es=causalsind[j]/causalind[j]
  
  ################ pte1 given data2  
  xob=data2[,1];deltaob=data2[,2];aob=data2[,3];sob=data2[,4];n=n2;v=vv[-indexindex]
  
  from = min(sob[sob!=0],na.rm = T); to = quantile(sob[sob!=0],.95,na.rm = T); step=((to - from)/nn)
  s=seq(from, to, by = step)
  
  #### optimal function 
  {
  temp=WEIGHT.p(xob,deltaob,aob,n=n,v)
  wei.t0=temp[,1];wei.t=temp[,2]
  
  bw = 1.06*sd(sob,na.rm=T)*n^(-1/5)/(n^0.06)
  kern = Kern.FUN(zz=s,zi=sob,bw)
  
  f0.s.t0.t0.hat=apply(as.numeric(v)*(aob==0)*(xob>t.0)*(sob<=t.0)*kern*wei.t0,2,sum,na.rm=T)/sum(as.numeric(v)*(aob==0)*wei.t0)
  f0.s.t.t0.hat=apply(as.numeric(v)*(aob==0)*(xob>t)*(sob<=t.0)*kern*wei.t,2,sum,na.rm=T)/sum(as.numeric(v)*(aob==0)*wei.t)
  f1.s.t0.t0.hat=apply(as.numeric(v)*(aob==1)*(xob>t.0)*(sob<=t.0)*kern*wei.t0,2,sum,na.rm=T)/sum(as.numeric(v)*(aob==1)*wei.t0)
  f1.s.t.t0.hat=apply(as.numeric(v)*(aob==1)*(xob>t)*(sob<=t.0)*kern*wei.t,2,sum,na.rm=T)/sum(as.numeric(v)*(aob==1)*wei.t)
  temp=(sob>t.0); temp[is.na(temp)]=1
  p0.t0.t0.hat=sum(as.numeric(v)*(aob==0)*(xob>t.0)*temp*wei.t0, na.rm = T)/sum(as.numeric(v)*(aob==0)*wei.t0)
  p0.t.t0.hat=sum(as.numeric(v)*(aob==0)*(xob>t)*temp*wei.t, na.rm = T)/sum(as.numeric(v)*(aob==0)*wei.t)
  p1.t0.t0.hat=sum(as.numeric(v)*(aob==1)*(xob>t.0)*temp*wei.t0, na.rm = T)/sum(as.numeric(v)*(aob==1)*wei.t0)
  p1.t.t0.hat=sum(as.numeric(v)*(aob==1)*(xob>t)*temp*wei.t, na.rm = T)/sum(as.numeric(v)*(aob==1)*wei.t)
  
  integrand<-f0.s.t0.t0.hat^2/f1.s.t0.t0.hat
  temp=(integrand[1] + integrand[nn+1] + 2*sum(integrand[seq(2,nn,by=2)]) + 4 *sum(integrand[seq(3,nn-1, by=2)]) )*step/3
  de=temp+p0.t0.t0.hat^2/p1.t0.t0.hat
  
  mu0t=sum(as.numeric(v)*(aob==0)*(xob>t)*wei.t)/sum(as.numeric(v)*(aob==0)*wei.t)
  mu1t=sum(as.numeric(v)*(aob==1)*(xob>t)*wei.t)/sum(as.numeric(v)*(aob==1)*wei.t)
  
  integrand<-f0.s.t0.t0.hat*f1.s.t.t0.hat/f1.s.t0.t0.hat
  temp=(integrand[1] + integrand[nn+1] + 2*sum(integrand[seq(2,nn,by=2)]) + 4 *sum(integrand[seq(3,nn-1, by=2)]) )*step/3
  nu=mu0t-temp-p0.t0.t0.hat*p1.t.t0.hat/p1.t0.t0.hat
  
  lambda=nu/de
  
  gs.hat2=(lambda*f0.s.t0.t0.hat+f1.s.t.t0.hat)/f1.s.t0.t0.hat
  g3.hat2=(lambda*p0.t0.t0.hat+p1.t.t0.hat)/p1.t0.t0.hat
  }
  
  #### pte
  xob=data1[,1];deltaob=data1[,2];aob=data1[,3];sob=data1[,4];n=n1;v=vv[indexindex]
  temp=WEIGHT.p(xob,deltaob,aob,n=n,v)
  wei.t0=temp[,1];wei.t=temp[,2]
  
  # causal=mu1t-mu0t
  causal.t0=(causal.t0+sum(as.numeric(v)*(xob>t.0)*aob*wei.t0)/sum(as.numeric(v)*aob*wei.t0)-
               sum(as.numeric(v)*(xob>t.0)*(1-aob)*wei.t0)/sum(as.numeric(v)*(1-aob)*wei.t0))/2
  causal2[j]=sum(as.numeric(v)*(xob>t)*aob*wei.t)/sum(as.numeric(v)*aob*wei.t)-sum(as.numeric(v)*(xob>t)*(1-aob)*wei.t)/sum(as.numeric(v)*(1-aob)*wei.t)
  tempind=c(sapply(1:n, function(kk){which.min(abs(sob[kk]-s))})); tempind=as.numeric(tempind)
  temp=(sob>t.0); temp[is.na(temp)]=1
  causals2[j]=sum(as.numeric(v)*aob*wei.t0*(xob>t.0)*(sob<=t.0)*gs.hat2[tempind],na.rm = T)/sum(as.numeric(v)*aob*wei.t0)+
    sum(as.numeric(v)*(xob>t.0)*temp*g3.hat2*aob*wei.t0)/sum(as.numeric(v)*aob*wei.t0) -
    (sum(as.numeric(v)*(1-aob)*wei.t0*(xob>t.0)*(sob<=t.0)*gs.hat2[tempind],na.rm = T)/sum(as.numeric(v)*(1-aob)*wei.t0)+
       sum(as.numeric(v)*(xob>t.0)*temp*g3.hat2*(1-aob)*wei.t0)/sum(as.numeric(v)*(1-aob)*wei.t0))
  pte.es=(pte.es+causals2[j]/causal2[j])/2
  
  # causalind2[j]=sum(as.numeric(v)*(xob>t)*aob*wei.t)/sum(as.numeric(v)*aob*wei.t)-sum(as.numeric(v)*(xob>t)*(1-aob)*wei.t)/sum(as.numeric(v)*(1-aob)*wei.t)
  # temp=(xob>t.0); temp[is.na(temp)]=1
  # causalsind2[j]=sum(as.numeric(v)*(xob>t.0)*temp*g.hat*aob*wei.t0)/sum(as.numeric(v)*aob*wei.t0) -
  #   sum(as.numeric(v)*(xob>t.0)*temp*g.hat*(1-aob)*wei.t0)/sum(as.numeric(v)*(1-aob)*wei.t0)
  # pteind.es=(pteind.es+causalsind2[j]/causalind2[j])/2
  
  g3.hat=(g3.hat+g3.hat2)/2
  gs.hat=(gs.hat+gs.hat2)/2
  }
  
  # mm=length(causal)-1
  # integrand=causal
  # rmst=(integrand[1] + integrand[mm+1] + 2*sum(integrand[seq(2,mm,by=2)]) + 4 *sum(integrand[seq(3,mm-1, by=2)]) )*stept/3
  # integrand=causals
  # rmsts=(integrand[1] + integrand[mm+1] + 2*sum(integrand[seq(2,mm,by=2)]) + 4 *sum(integrand[seq(3,mm-1, by=2)]) )*stept/3
  # ptermst=(rmsts+causal.t0)/(rmst+causal.t0)
  # 
  # integrand=causal2
  # rmst=(integrand[1] + integrand[mm+1] + 2*sum(integrand[seq(2,mm,by=2)]) + 4 *sum(integrand[seq(3,mm-1, by=2)]) )*stept/3
  # integrand=causals2
  # rmsts=(integrand[1] + integrand[mm+1] + 2*sum(integrand[seq(2,mm,by=2)]) + 4 *sum(integrand[seq(3,mm-1, by=2)]) )*stept/3
  # ptermst=(ptermst+(rmsts+causal.t0)/(rmst+causal.t0))/2
  # 
  # integrand=causalind
  # rmst=(integrand[1] + integrand[mm+1] + 2*sum(integrand[seq(2,mm,by=2)]) + 4 *sum(integrand[seq(3,mm-1, by=2)]) )*stept/3
  # integrand=causalsind
  # rmsts=(integrand[1] + integrand[mm+1] + 2*sum(integrand[seq(2,mm,by=2)]) + 4 *sum(integrand[seq(3,mm-1, by=2)]) )*stept/3
  # ptermstind=(rmsts+causal.t0)/(rmst+causal.t0)
  # 
  # integrand=causalind2
  # rmst=(integrand[1] + integrand[mm+1] + 2*sum(integrand[seq(2,mm,by=2)]) + 4 *sum(integrand[seq(3,mm-1, by=2)]) )*stept/3
  # integrand=causalsind2
  # rmsts=(integrand[1] + integrand[mm+1] + 2*sum(integrand[seq(2,mm,by=2)]) + 4 *sum(integrand[seq(3,mm-1, by=2)]) )*stept/3
  # ptermstind=(ptermstind+(rmsts+causal.t0)/(rmst+causal.t0))/2
  
#### output
  out=c(pte.es,g3.hat,gs.hat)#,pteind.es,ptermst,ptermstind
  
}


gen.perturb.weights=function(data.num, n, num.perturb=500){
  set.seed(data.num)
  #matrix(rexp(n*num.perturb, rate=1), nrow=n, ncol=num.perturb)
  index = sapply(1:num.perturb,function(x) sample(1:n,n,replace=T))
  apply(index,2,function(x) tabulate(x,nbins=n))
}


gen.bootstrap.weights=function(data.num, n, num.perturb=500){
  set.seed(data.num)
  sapply(1:num.perturb,function(x) sample(1:n,n,replace=T))
  #weights = apply(index,2,function(x) tabulate(x,nbins=n))
  #list(index=index,weights=weights)
}


WEIGHT<-function(xob,deltaob,aob,n){
  x=xob[aob==0]; delta=1-deltaob[aob==0]
  xsort=sort(x[delta==1])
  risk=VTM(x,length(xsort))>=xsort
  risk.n=apply(risk,1,sum)
  Lam=cumsum(1/risk.n)
  s0=exp(-Lam)
  sur0=data.frame("time"=xsort,"surv"=s0)
  x=xob[aob==1]; delta=1-deltaob[aob==1]
  xsort=sort(x[delta==1])
  risk=VTM(x,length(xsort))>=xsort
  risk.n=apply(risk,1,sum)
  Lam=cumsum(1/risk.n)
  s1=exp(-Lam)
  sur1=data.frame("time"=xsort,"surv"=s1)
  ind0=c(sapply(1:n, function(kk){which.min(abs(xob[kk]-sur0$time))}))
  G0=sur0$surv[ind0]
  ind1=c(sapply(1:n, function(kk){which.min(abs(xob[kk]-sur1$time))}))
  G1=sur1$surv[ind1]
  G=(1-aob)*G0+aob*G1
  G.t0=(1-aob)*G0[which.min(abs(t.0-xob))]+aob*G1[which.min(abs(t.0-xob))]
  G.t=(1-aob)*G0[which.min(abs(t-xob))]+aob*G1[which.min(abs(t-xob))]
  wei.t0=(xob<=t.0)*deltaob/G+(xob>t.0)/G.t0; #wei.t0[is.nan(wei.t0)]=max(wei.t0[!is.nan(wei.t0)])
  wei.t=(xob<=t)*deltaob/G+(xob>t)/G.t; #wei.t[is.nan(wei.t)]=max(wei.t[!is.nan(wei.t)])
  
  out=cbind(wei.t0,wei.t,G.t,G.t0,G1,G0,G)#,pte2
}


WEIGHT.p<-function(xob,deltaob,aob,n,v){
  x=xob[aob==0];delta=1-deltaob[aob==0]
  xsort=sort(x[delta==1])
  risk=VTM(x,length(xsort))>=xsort
  risk.n=apply(risk*VTM(v[aob==0],length(xsort)),1,sum)
  N=(VTM(x,length(xsort))<=xsort)*VTM(delta,length(xsort))
  dN=rbind( N[1,],N[-1,]-N[-length(xsort),] )
  nu=apply(VTM(v[aob==0],length(xsort))*dN, 1, sum)
  Lam=cumsum(nu/risk.n)
  s0=exp(-Lam)
  sur0=data.frame("time"=xsort,"surv"=s0)
  x=xob[aob==1];delta=1-deltaob[aob==1]
  xsort=sort(x[delta==1])
  risk=VTM(x,length(xsort))>=xsort
  risk.n=apply(risk*VTM(v[aob==1],length(xsort)),1,sum)
  N=(VTM(x,length(xsort))<=xsort)*VTM(delta,length(xsort))
  dN=rbind( N[1,],N[-1,]-N[-length(xsort),] )
  nu=apply(dN*VTM(v[aob==1],length(xsort)), 1, sum)
  Lam=cumsum(nu/risk.n)
  s1=exp(-Lam)
  sur1=data.frame("time"=xsort,"surv"=s1)
  ind0=c(sapply(1:n, function(kk){which.min(abs(xob[kk]-sur0$time))}))
  G0=sur0$surv[ind0]
  ind1=c(sapply(1:n, function(kk){which.min(abs(xob[kk]-sur1$time))}))
  G1=sur1$surv[ind1]
  G=(1-aob)*G0+aob*G1
  G.t0=(1-aob)*G0[which.min(abs(t.0-xob))]+aob*G1[which.min(abs(t.0-xob))]
  G.t=(1-aob)*G0[which.min(abs(t-xob))]+aob*G1[which.min(abs(t-xob))]
  wei.t0=(xob<=t.0)*deltaob/G+(xob>t.0)/G.t0; #wei.t0[is.nan(wei.t0)]=max(wei.t0[!is.nan(wei.t0)])
  wei.t=(xob<=t)*deltaob/G+(xob>t)/G.t; #wei.t[is.nan(wei.t)]=max(wei.t[!is.nan(wei.t)])

  out=cbind(wei.t0,wei.t,G.t,G.t0,G1,G0,G)#,pte2
}


est<- function(t,t.0,tt,nn=200,re=5,data){
  # set.seed(2023)
  n=nrow(data)
  indexindex=sample(n, n/2, replace = FALSE)
  data1=data[indexindex,]
  data2=data[-indexindex,]
  n1=nrow(data1)
  n2=nrow(data2)
  
  causal=rep(NA,length(tt)); causal2=rep(NA,length(tt))
  causals=rep(NA,length(tt)); causals2=rep(NA,length(tt))
  causalind=rep(NA,length(tt)); causalind2=rep(NA,length(tt))
  causalsind=rep(NA,length(tt)); causalsind2=rep(NA,length(tt))
  for (j in 1:length(tt)){
    t=tt[j]
    ################ pte2 given data1  
    xob=data1[,1];deltaob=data1[,2];aob=data1[,3];sob=data1[,4];n=n1
    
    ####
    from = min(sob[sob!=0],na.rm = T); to = quantile(sob[sob!=0],.95,na.rm = T); step=((to - from)/nn)
    s=seq(from, to, by = step)
    # s.ob=(s.ob-mean(s.ob))/sd(s.ob)
    # s.ob=pnorm(s.ob, mean = 0, sd = 1)
    # from = 0.01; to = .99; step=((to - from)/nn)
    # s=seq(from, to, by = step)
    
    #### optimal function 
    {
      temp=WEIGHT(xob,deltaob,aob,n=n)
      wei.t0=temp[,1];wei.t=temp[,2]
      
      bw = 1.06*sd(sob,na.rm=T)*n^(-1/5)/(n^0.06)
      kern = Kern.FUN(zz=s,zi=sob,bw)
      
      f0.s.t0.t0.hat=apply((aob==0)*(xob>t.0)*(sob<=t.0)*kern*wei.t0,2,sum,na.rm=T)/sum((aob==0)*wei.t0)
      f0.s.t.t0.hat=apply((aob==0)*(xob>t)*(sob<=t.0)*kern*wei.t,2,sum,na.rm=T)/sum((aob==0)*wei.t)
      f1.s.t0.t0.hat=apply((aob==1)*(xob>t.0)*(sob<=t.0)*kern*wei.t0,2,sum,na.rm=T)/sum((aob==1)*wei.t0)
      f1.s.t.t0.hat=apply((aob==1)*(xob>t)*(sob<=t.0)*kern*wei.t,2,sum,na.rm=T)/sum((aob==1)*wei.t)
      temp=(sob>t.0); temp[is.na(temp)]=1
      p0.t0.t0.hat=sum((aob==0)*(xob>t.0)*temp*wei.t0, na.rm = T)/sum((aob==0)*wei.t0)
      p0.t.t0.hat=sum((aob==0)*(xob>t)*temp*wei.t, na.rm = T)/sum((aob==0)*wei.t)
      p1.t0.t0.hat=sum((aob==1)*(xob>t.0)*temp*wei.t0, na.rm = T)/sum((aob==1)*wei.t0)
      p1.t.t0.hat=sum((aob==1)*(xob>t)*temp*wei.t, na.rm = T)/sum((aob==1)*wei.t)
      
      integrand<-f0.s.t0.t0.hat^2/f1.s.t0.t0.hat
      temp=(integrand[1] + integrand[nn+1] + 2*sum(integrand[seq(2,nn,by=2)]) + 4 *sum(integrand[seq(3,nn-1, by=2)]) )*step/3
      de=temp+p0.t0.t0.hat^2/p1.t0.t0.hat
      
      mu0t=sum((aob==0)*(xob>t)*wei.t)/sum((aob==0)*wei.t)
      mu1t=sum((aob==1)*(xob>t)*wei.t)/sum((aob==1)*wei.t)
      mu0t0=sum((aob==0)*(xob>t.0)*wei.t0)/sum((aob==0)*wei.t0)
      mu1t0=sum((aob==1)*(xob>t.0)*wei.t0)/sum((aob==1)*wei.t0)
      
      integrand<-f0.s.t0.t0.hat*f1.s.t.t0.hat/f1.s.t0.t0.hat
      temp=(integrand[1] + integrand[nn+1] + 2*sum(integrand[seq(2,nn,by=2)]) + 4 *sum(integrand[seq(3,nn-1, by=2)]) )*step/3
      nu=mu0t-temp-p0.t0.t0.hat*p1.t.t0.hat/p1.t0.t0.hat
      
      lambda=nu/de
      lambda2=(mu0t-mu0t0*mu1t/mu1t0)/(mu0t0^2/mu1t0)
      
      gs.hat=(lambda*f0.s.t0.t0.hat+f1.s.t.t0.hat)/f1.s.t0.t0.hat
      g3.hat=(lambda*p0.t0.t0.hat+p1.t.t0.hat)/p1.t0.t0.hat
      g.hat=(lambda2*mu0t0+mu1t)/mu1t0
    }
    
    #### pte
    xob=data2[,1];deltaob=data2[,2];aob=data2[,3];sob=data2[,4];n=n2
    temp=WEIGHT(xob,deltaob,aob,n=n)
    wei.t0=temp[,1];wei.t=temp[,2]
    
    # causal=mu1t-mu0t
    causal.t0=sum((xob>t.0)*aob*wei.t0)/sum(aob*wei.t0)-sum((xob>t.0)*(1-aob)*wei.t0)/sum((1-aob)*wei.t0)
    causal[j]=sum((xob>t)*aob*wei.t)/sum(aob*wei.t)-sum((xob>t)*(1-aob)*wei.t)/sum((1-aob)*wei.t)
    tempind=c(sapply(1:n, function(kk){which.min(abs(sob[kk]-s))})); tempind=as.numeric(tempind)
    temp=(sob>t.0); temp[is.na(temp)]=1
    causals[j]=sum(aob*wei.t0*(xob>t.0)*(sob<=t.0)*gs.hat[tempind],na.rm = T)/sum(aob*wei.t0)+
      sum((xob>t.0)*temp*g3.hat*aob*wei.t0)/sum(aob*wei.t0) -
      (sum((1-aob)*wei.t0*(xob>t.0)*(sob<=t.0)*gs.hat[tempind],na.rm = T)/sum((1-aob)*wei.t0)+
         sum((xob>t.0)*temp*g3.hat*(1-aob)*wei.t0)/sum((1-aob)*wei.t0))
    pte.es=causals[j]/causal[j]
    
    causalind[j]=sum((xob>t)*aob*wei.t)/sum(aob*wei.t)-sum((xob>t)*(1-aob)*wei.t)/sum((1-aob)*wei.t)
    temp=(xob>t.0); temp[is.na(temp)]=1
    causalsind[j]=sum((xob>t.0)*temp*g.hat*aob*wei.t0)/sum(aob*wei.t0) -
      sum((xob>t.0)*temp*g.hat*(1-aob)*wei.t0)/sum((1-aob)*wei.t0)
    pteind.es=causalsind[j]/causalind[j]
    
    ################ pte1 given data2  
    xob=data2[,1];deltaob=data2[,2];aob=data2[,3];sob=data2[,4];n=n2
    
    ####
    from = min(sob[sob!=0],na.rm = T); to = quantile(sob[sob!=0],.95,na.rm = T); step=((to - from)/nn)
    s=seq(from, to, by = step)
    # s.ob=(s.ob-mean(s.ob))/sd(s.ob)
    # s.ob=pnorm(s.ob, mean = 0, sd = 1)
    # from = 0.01; to = .99; step=((to - from)/nn)
    # s=seq(from, to, by = step)
    
    #### optimal function 
    {
      temp=WEIGHT(xob,deltaob,aob,n=n)
      wei.t0=temp[,1];wei.t=temp[,2]
      
      bw = 1.06*sd(sob,na.rm=T)*n^(-1/5)/(n^0.06)
      kern = Kern.FUN(zz=s,zi=sob,bw)
      
      f0.s.t0.t0.hat=apply((aob==0)*(xob>t.0)*(sob<=t.0)*kern*wei.t0,2,sum,na.rm=T)/sum((aob==0)*wei.t0)
      f0.s.t.t0.hat=apply((aob==0)*(xob>t)*(sob<=t.0)*kern*wei.t,2,sum,na.rm=T)/sum((aob==0)*wei.t)
      f1.s.t0.t0.hat=apply((aob==1)*(xob>t.0)*(sob<=t.0)*kern*wei.t0,2,sum,na.rm=T)/sum((aob==1)*wei.t0)
      f1.s.t.t0.hat=apply((aob==1)*(xob>t)*(sob<=t.0)*kern*wei.t,2,sum,na.rm=T)/sum((aob==1)*wei.t)
      temp=(sob>t.0); temp[is.na(temp)]=1
      p0.t0.t0.hat=sum((aob==0)*(xob>t.0)*temp*wei.t0, na.rm = T)/sum((aob==0)*wei.t0)
      p0.t.t0.hat=sum((aob==0)*(xob>t)*temp*wei.t, na.rm = T)/sum((aob==0)*wei.t)
      p1.t0.t0.hat=sum((aob==1)*(xob>t.0)*temp*wei.t0, na.rm = T)/sum((aob==1)*wei.t0)
      p1.t.t0.hat=sum((aob==1)*(xob>t)*temp*wei.t, na.rm = T)/sum((aob==1)*wei.t)
      
      integrand<-f0.s.t0.t0.hat^2/f1.s.t0.t0.hat
      temp=(integrand[1] + integrand[nn+1] + 2*sum(integrand[seq(2,nn,by=2)]) + 4 *sum(integrand[seq(3,nn-1, by=2)]) )*step/3
      de=temp+p0.t0.t0.hat^2/p1.t0.t0.hat
      
      mu0t=sum((aob==0)*(xob>t)*wei.t)/sum((aob==0)*wei.t)
      mu1t=sum((aob==1)*(xob>t)*wei.t)/sum((aob==1)*wei.t)
      mu0t0=sum((aob==0)*(xob>t.0)*wei.t0)/sum((aob==0)*wei.t0)
      mu1t0=sum((aob==1)*(xob>t.0)*wei.t0)/sum((aob==1)*wei.t0)
      
      integrand<-f0.s.t0.t0.hat*f1.s.t.t0.hat/f1.s.t0.t0.hat
      temp=(integrand[1] + integrand[nn+1] + 2*sum(integrand[seq(2,nn,by=2)]) + 4 *sum(integrand[seq(3,nn-1, by=2)]) )*step/3
      nu=mu0t-temp-p0.t0.t0.hat*p1.t.t0.hat/p1.t0.t0.hat
      
      lambda=nu/de
      lambda2=(mu0t-mu0t0*mu1t/mu1t0)/(mu0t0^2/mu1t0)
      
      gs.hat2=(lambda*f0.s.t0.t0.hat+f1.s.t.t0.hat)/f1.s.t0.t0.hat
      g3.hat2=(lambda*p0.t0.t0.hat+p1.t.t0.hat)/p1.t0.t0.hat
      g.hat2=(lambda2*mu0t0+mu1t)/mu1t0
      g3.es=(g3.hat+g3.hat2)/2
      gs.es=(gs.hat+gs.hat2)/2
      g.es=(g.hat+g.hat2)/2
    }
    
    #### pte
    xob=data1[,1];deltaob=data1[,2];aob=data1[,3];sob=data1[,4];n=n1
    temp=WEIGHT(xob,deltaob,aob,n=n)
    wei.t0=temp[,1];wei.t=temp[,2]
    
    # causal=mu1t-mu0t
    causal.t0=(causal.t0+sum((xob>t.0)*aob*wei.t0)/sum(aob*wei.t0)-sum((xob>t.0)*(1-aob)*wei.t0)/sum((1-aob)*wei.t0))/2
    causal2[j]=sum((xob>t)*aob*wei.t)/sum(aob*wei.t)-sum((xob>t)*(1-aob)*wei.t)/sum((1-aob)*wei.t)
    tempind=c(sapply(1:n, function(kk){which.min(abs(sob[kk]-s))})); tempind=as.numeric(tempind)
    temp=(sob>t.0); temp[is.na(temp)]=1
    causals2[j]=sum(aob*wei.t0*(xob>t.0)*(sob<=t.0)*gs.hat2[tempind],na.rm = T)/sum(aob*wei.t0)+
      sum((xob>t.0)*temp*g3.hat2*aob*wei.t0)/sum(aob*wei.t0) -
      (sum((1-aob)*wei.t0*(xob>t.0)*(sob<=t.0)*gs.hat2[tempind],na.rm = T)/sum((1-aob)*wei.t0)+
         sum((xob>t.0)*temp*g3.hat2*(1-aob)*wei.t0)/sum((1-aob)*wei.t0))
    pte.es=(pte.es+causals2[j]/causal2[j])/2
    
    causalind2[j]=sum((xob>t)*aob*wei.t)/sum(aob*wei.t)-sum((xob>t)*(1-aob)*wei.t)/sum((1-aob)*wei.t)
    temp=(xob>t.0); temp[is.na(temp)]=1
    causalsind2[j]=sum((xob>t.0)*temp*g.hat2*aob*wei.t0)/sum(aob*wei.t0) -
      sum((xob>t.0)*temp*g.hat2*(1-aob)*wei.t0)/sum((1-aob)*wei.t0)
    pteind.es=(pteind.es+causalsind2[j]/causalind2[j])/2
  }
  
  mm=length(causal)-1
  integrand=causal
  rmst=(integrand[1] + integrand[mm+1] + 2*sum(integrand[seq(2,mm,by=2)]) + 4 *sum(integrand[seq(3,mm-1, by=2)]) )*stept/3
  integrand=causals
  rmsts=(integrand[1] + integrand[mm+1] + 2*sum(integrand[seq(2,mm,by=2)]) + 4 *sum(integrand[seq(3,mm-1, by=2)]) )*stept/3
  ptermst=(rmsts+causal.t0)/(rmst+causal.t0)
  
  integrand=causal2
  rmst=(integrand[1] + integrand[mm+1] + 2*sum(integrand[seq(2,mm,by=2)]) + 4 *sum(integrand[seq(3,mm-1, by=2)]) )*stept/3
  integrand=causals2
  rmsts=(integrand[1] + integrand[mm+1] + 2*sum(integrand[seq(2,mm,by=2)]) + 4 *sum(integrand[seq(3,mm-1, by=2)]) )*stept/3
  ptermst=(ptermst+(rmsts+causal.t0)/(rmst+causal.t0))/2
  
  integrand=causalind
  rmst=(integrand[1] + integrand[mm+1] + 2*sum(integrand[seq(2,mm,by=2)]) + 4 *sum(integrand[seq(3,mm-1, by=2)]) )*stept/3
  integrand=causalsind
  rmsts=(integrand[1] + integrand[mm+1] + 2*sum(integrand[seq(2,mm,by=2)]) + 4 *sum(integrand[seq(3,mm-1, by=2)]) )*stept/3
  ptermstind=(rmsts+causal.t0)/(rmst+causal.t0)
  
  integrand=causalind2
  rmst=(integrand[1] + integrand[mm+1] + 2*sum(integrand[seq(2,mm,by=2)]) + 4 *sum(integrand[seq(3,mm-1, by=2)]) )*stept/3
  integrand=causalsind2
  rmsts=(integrand[1] + integrand[mm+1] + 2*sum(integrand[seq(2,mm,by=2)]) + 4 *sum(integrand[seq(3,mm-1, by=2)]) )*stept/3
  ptermstind=(ptermstind+(rmsts+causal.t0)/(rmst+causal.t0))/2
  
  causal.est=(causal[length(causal)]+causal2[length(causal2)])/2
  causals.est=(causals[length(causals)]+causals2[length(causals2)])/2
  
  # ################ layla
  # temp=R.q.event(xone = data$xob[data$aob==1], xzero = data$xob[data$aob==0], deltaone = data$deltaob[data$aob==1],
  #           deltazero = data$deltaob[data$aob==0], sone = data$sob[data$aob==1], szero = data$sob[data$aob==0], t = 5,
  #           landmark=t.0, type = "np")
  # ptela[i]=temp$R.q
  
  
  ######## variance of proposed estimators
  vv=matrix(rexp((n1+n2)*re),nrow=n1+n2)
  temp=apply(vv,2,resam,t,t.0,tt=5,data,data1,data2,indexindex)
  aa=which(temp[1,]>1 | temp[1,]<0)
  if (length(aa)>= re-1){
    pte.se=NA
    # pteind.se=NA
    # ptermst.se=NA
    # ptermstind.se=NA
    g3.se=NA
    gs.se=rep(NA,nn+1)
  } else {
    pte.se=sd(temp[1,temp[1,]<=1 & temp[1,]>=0])
    # pteind.se=sd(temp[2,temp[1,]<=1 & temp[1,]>=0])
    # ptermst.se=sd(temp[3,temp[3,]<=1 & temp[3,]>=0])
    # ptermstind.se=sd(temp[4,temp[3,]<=1 & temp[3,]>=0])
    g3.se=sd(temp[2,temp[1,]<=1 & temp[1,]>=0])
    gs.se=apply( temp[-(1:2),temp[1,]<=1 & temp[1,]>=0],1,sd)
  }
  
  # result=data.frame( matrix(c(pte.es[i],pteind.es[i],ptermst[i],ptermstind[i],g3.es[i],
  #                             pte.se[i],pteind.se[i],ptermst.se[i],ptermstind.se[i],g3.se[i],ptela[i]),nrow=1) )
  # write.table(result, paste(getwd(),"/set",set,'_t0_',t.0,".txt", sep =""),
  #             append=TRUE,col.names = FALSE, row.names= FALSE)
  out=list('pte.es'=pte.es,'pteind.es'=pteind.es,'ptermst.es'=ptermst,'ptermstind.es'=ptermstind,
           'g3.es'=g3.es,'gs.es'=gs.es,
           'pte.se'=pte.se,#'pteind.se'=pteind.se,'ptermst.se'=ptermst.se,'ptermstind.se'=ptermstind.se,
           'g3.se'=g3.se,'gs.se'=gs.se)
}
