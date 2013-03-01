
lesum<-function(lx) { # finds log of sum of exponentials
        max.l=which.max(lx);
        out=lx[max.l]+log(sum(exp(-abs(lx-lx[max.l]))));
        return(out);
     }

mleNB <- function(x,eps.nb=1e-8){
require(MASS);
if ((length(x)>1)&(mean(x)>0)) out=fitdistr(x,"Negative Binomial",control=list(reltol=eps.nb))$estimate else out=c(100,x[1])+eps.nb;
return(out);
}

postCPcrit<-function(data,seg=numeric(),model,eps.nb=1e-8, prior=0.5, prior.type="n"){

  k=length(seg)+1;
  print(paste("Number of segments:",k));
  n=length(data);
  bestcp=numeric();
  if (!is.element(model,1:3)) stop("Please enter model= 1 (Poisson), 2 (normal), 3 (negative binomial)");
  if (k==1){
    if (model==2){
      mBIC=0;
      BIC=-sum(dnorm(data,mean(data),sd(data)*sqrt((n-1)/n),TRUE))+(2)*log(n);
      entropy=0;
      AIC=-sum(dnorm(data,mean(data),sd(data)*sqrt((n-1)/n),TRUE))+(2)*2;
    } 
    if (model==1){
      mBIC=-n*mean(data)*log(mean(data))+log(n);
      BIC=-sum(dpois(data,mean(data),TRUE))+log(n);
      AIC=-sum(dpois(data,mean(data),TRUE))+2;
      entropy=0;
    } 
    if (model==3){
      negbin.est=mleNB(data,eps.nb=eps.nb);
      BIC=-sum(dnbinom(data,mu=negbin.est[2],size=negbin.est[1],log=TRUE))+1*log(n);
      AIC=-sum(dnbinom(data,mu=negbin.est[2],size=negbin.est[1],log=TRUE))+1*2;
      entropy=0;
    } 
  }else{
   if (model==3){
    state.temp=rep(1:k,diff(c(0,seg,n)));
    out.mle=lapply(split(data,state.temp),mleNB,eps.nb=eps.nb);
    out.sizes=sapply(out.mle,function(x) x[1]);
    out.means=sapply(out.mle,function(x) x[2]);
    # log-density of neg binomial, 1 column for each possible state
    lprob.matrix=matrix(dnbinom(rep(data,k),size=out.sizes[rep(1:k,each=n)],mu=out.means[rep(1:k,each=n)],log=TRUE),ncol=k);
    # initial run of postCP
    out=viterbi(lprob=lprob.matrix,verbose=FALSE,prior=prior,prior.type=prior.type);
    bestcp=out$bestcp;
    state.temp=rep(1:k,diff(c(0,out$bestcp,n)));
    out.mle=lapply(split(data,state.temp),mleNB,eps.nb=eps.nb);
    out.sizes=sapply(out.mle,function(x) x[1]);
    out.means=sapply(out.mle,function(x) x[2]);
    lprob.matrix=matrix(dnbinom(rep(data,k),size=out.sizes[rep(1:k,each=n)],mu=out.means[rep(1:k,each=n)],log=TRUE),ncol=k);
    # re-run postCP on aposteriori
    rm(out);
    out=postCP(lprob=lprob.matrix,keep=TRUE,verbose=FALSE,initsegci=FALSE,ci=0,prior=prior,prior.type=prior.type);
   
    BIC=-sum(dnbinom(data,size=out.sizes[state.temp],mu=out.means[state.temp],log=TRUE))+(k*1)*log(n);
    AIC=-sum(dnbinom(data,size=out.sizes[state.temp],mu=out.means[state.temp],log=TRUE))+(k*1)*2;
   }else{
      out=viterbi(data,seg,model,verbose=FALSE,prior=prior,prior.type=prior.type);
      bestcp=out$bestcp;
     # redo postCP using a posteriori bestCP (viterbi) from previous postCP run
      state.temp=rep(1:k,c(bestcp[1],diff(bestcp),n-max(bestcp)));
      out.means=aggregate(data,list(state.temp),mean)$x;
      out=postCP(data,bestcp,model,keep=TRUE,verbose=FALSE,debug=FALSE,ci=0,prior=prior,prior.type=prior.type);
      if (model==2) {
        # use for normal data
	SSB=sum(table(state.temp)*(out.means-mean(data))^2);
        SSA=sum((data-mean(data))^2);
        mBIC=-(0.5*(n-k+2)*log(1+SSB/(SSA-SSB))+lgamma(0.5*(n-k+2))-lgamma(0.5*(n+1))+0.5*(k-1)*log(SSA)-0.5*sum(log(table(state.temp)))+(1.5-k)*log(n));
        BIC=-sum(dnorm(data,out.means[state.temp],out$sds,TRUE))+(k+1)*log(n);
        AIC=-sum(dnorm(data,out.means[state.temp],out$sds,TRUE))+(k+1)*2;
      }
      if (model==1){
    # use for Poisson data
        mBIC=-(sum(table(state.temp)*out.means*log(out.means))-0.5*sum(log(table(state.temp)))+(0.5-k)*log(n));
        BIC=-sum(dpois(data,out.means[state.temp],TRUE))+(k)*log(n);
        AIC=-sum(dpois(data,out.means[state.temp],TRUE))+(k)*2;
      }
    }
    if (model==1){
       ldensmatrix=matrix(dpois(rep(data,k),out$means[rep(1:k,each=n)],log=TRUE),ncol=k);
    }
    if (model==2){
       ldensmatrix=matrix(dnorm(rep(data,k),out$means[rep(1:k,each=n)],out$sds,log=TRUE),ncol=k)+log(out$sds);
    } 
    if (model==3){
       ldensmatrix=matrix(dnbinom(rep(data,k),size=out.sizes[rep(1:k,each=n)],mu=out.means[rep(1:k,each=n)],log=TRUE),ncol=k);
    }    
   if (out$prior.type==1) prior=out$prior;
   if (out$prior.type==2) prior=out$prior;
   if (out$prior.type==3) prior=matrix(rep(out$prior,each=n-1),ncol=k-1);
   lpost.trans=ldensmatrix[-1,-1]+out$lbackward[-1,-1]-out$lbackward[-n,-k]+log(prior);   
   entropy=-sum(out$post.state[1,1]*log(out$post.state[1,1]),na.rm=TRUE)-sum(out$post.cp*lpost.trans,na.rm=TRUE);
  }
   ICL=BIC+entropy;
  if (model!=3) { scores=c(ICL,AIC,BIC,mBIC); names(scores)=c("ICL","AIC","BIC","mBIC")} else {scores=c(ICL,AIC,BIC); names(scores)=c("ICL","AIC","BIC")};
  results<-list(scores=scores,cp.loc=bestcp);
  results;
}

postCPmodelsel <- function(data,K.range,model=1,greedy=FALSE,eps.nb=1e-8){

  n=length(data);
  K.min=K.range[1];
  K.max=K.range[length(K.range)];
  entropy=numeric();
  BIC=numeric();
  mBIC=numeric();
  AIC=numeric();
  bestcp=list();
  k=K.min;
  if (length(K.range)<2) stop("Please choose models with at least two different values of K.");
  if (!greedy) {
     if (model==2) seg.matrix=cghseg:::segmeanCO(data, Kmax=K.max)$t.est else seg.matrix=getBreaks(Segmentor(data,model=model,Kmax=K.max));     
  }  else seg.matrix=numeric();

  modelsel1<- function(k,data,seg.matrix, model, eps.nb = eps.nb){
  	 if (k>1) {if (length(seg.matrix)==0) seg=GreedySegmente(data,k)[2:k]-1 else seg=seg.matrix[k,1:(k-1)];} else seg=numeric();
	out=postCPcrit(data,seg,model,eps.nb);
	return(out);
  }
  out=sapply(K.range,modelsel1,data=data,seg.matrix=seg.matrix,model=model,eps.nb=eps.nb);
  scores=unlist(out[1,]);
  ICL=scores[names(scores)=="ICL"];
  BIC=scores[names(scores)=="BIC"];
  AIC=scores[names(scores)=="AIC"];
  names(ICL)=K.range;
  names(AIC)=K.range;
  names(BIC)=K.range;
  if (model!=3) {
    mBIC=scores[names(scores)=="mBIC"];    names(mBIC)=K.range;
  }
  best.K.ICL <- which.min(ICL);
  best.K.AIC <- which.min(AIC);
  best.K.BIC <- which.min(BIC);
  bestcp=out[2,];
  cp.loc=list();
  cp.loc$ICL = bestcp[[best.K.ICL]];
  cp.loc$BIC = bestcp[[best.K.BIC]];
  cp.loc$AIC = bestcp[[best.K.AIC]];
  if (model!=3) { 
	names(mBIC)=K.range;
	best.K.mBIC=which.min(mBIC); 
	} else best.K.mBIC=NA;
  if (model!=3) cp.loc$mBIC = bestcp[[best.K.mBIC]];
   
   if (model!=3) { scores=cbind(ICL,AIC,BIC,mBIC); colnames(scores)=c("ICL","AIC","BIC","mBIC");rownames=K.range;} else {scores=cbind(ICL,AIC,BIC); colnames(scores)=c("ICL","AIC","BIC");rownames=K.range;};
  if (model!=3) results<-list(ICL=K.range[best.K.ICL],AIC=K.range[best.K.AIC],BIC=K.range[best.K.BIC],mBIC=K.range[best.K.mBIC],scores=scores,cp.loc=cp.loc)else{
    results<-list(ICL=K.range[best.K.ICL],AIC=K.range[best.K.AIC],BIC=K.range[best.K.BIC],scores=scores,cp.loc=cp.loc,eps.nb=eps.nb)
  }
  results;
}
