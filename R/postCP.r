postCPsample <- function(postCP.res, nsamples=100, gen.data=FALSE, prior=0.5, prior.type="n", verbose=TRUE,debug=FALSE) UseMethod("postCPsample")

postCPsample.default <- function(postCP.res, nsamples=100, gen.data=FALSE, prior=0.5, prior.type="n",verbose=TRUE,debug=FALSE){
# postCP.res: results from function postCP
# nsamples: number of sets of generated changepoints
# gen.data: generate a matrix of data, one row of length n for each of nsamples replicates
  if(missing(postCP.res)){
     stop("Need to use results of postCP first")
  }
  data=as.double(postCP.res$data)
  mu=as.double(postCP.res$means)
  sigma=as.double(postCP.res$sds)
  prior=as.double(postCP.res$prior)
  ptype=as.integer(postCP.res$prior.type)
  lprob=numeric()
  if (!is.element("lbackward",names(postCP.res))) {
    stop("\nRe-run postCP with keep=TRUE");
  }
  lbackward=as.double(as.vector(postCP.res$lbackward))
  if (!is.element("lprob",names(postCP.res))){
  n=postCP.res$n
  J=length(mu)
  model=1
  if (postCP.res$model=="normal") model=2} else{
    n=nrow(postCP.res$lprob);
    J=ncol(postCP.res$lprob);
    lprob=as.vector(postCP.res$lprob);
   gen.data=FALSE;
  }

 if (prior.type=="n") ptype=1 else
 if (prior.type=="o") ptype=2 else
 if (prior.type=="s") ptype=3 
  cpsize=nsamples*(J-1)
  cpvector=rep(0,cpsize)

  if (!is.element("lbackward",names(postCP.res))) {
  .C("postCPsample",rdata=data,rmu=mu,rsigma=sigma,rlbackward=as.double(lbackward),as.double(cpvector),rnsamples=as.integer(nsamples),nn=as.integer(n),JJ=as.integer(J),rmodel=as.integer(model),rprior=as.double(prior),rpriortype=as.integer(ptype),rverbose=as.integer(verbose),rdebug=as.integer(debug), DUP=FALSE,PACKAGE="postCP")
  cpmatrix=matrix(cpvector,ncol=J-1)} else{
  .C("postCPsampleext",lprob=as.double(lprob),rlbackward=as.double(lbackward),as.double(cpvector),rnsamples=as.integer(nsamples),nn=as.integer(n),JJ=as.integer(J),rprior=as.double(prior),rpriortype=as.integer(ptype),rverbose=as.integer(verbose),rdebug=as.integer(debug), DUP=FALSE,PACKAGE="postCP")
  cpmatrix=matrix(cpvector,ncol=J-1)}

  fbsample <- list(changepoints=cpmatrix)
  if (gen.data){
     getSegment <- function(cpv,nsamples) rowSums(outer((1:nsamples),cpv,">"))+1
     segments <- c(apply(cpmatrix,1,getSegment,nsamples=n))
     if (model==1) x=matrix(rpois(length(segments),mu[segments]),ncol=n,byrow=T)
     if (model==2) x=matrix(rnorm(length(segments),mu[segments],sigma),ncol=n,byrow=T)
     fbsample$data=x
  }
  class(fbsample) <- "postCPsample"
  fbsample
}

print.postCPsample <- function(x,...)
{ cat("\nGenerated change-points ($changepoints):\n")
  print(x[[1]])
  if (length(x)==2){
    cat("\nMatrix of generated data ($data):\n")
    str(x[[2]])
  }
}

postCP <- function(data=numeric(),seg=integer(),model=1,lprob=numeric(),keep=FALSE,ci=0.9,initsegci=TRUE,nsamples=0,gen.data=FALSE,prior=0.5,prior.type="n",verbose=TRUE,debug=FALSE) UseMethod("postCP")

postCP.default <- function(data=numeric(),seg=integer(),model=1,lprob=numeric(),keep=FALSE,ci=0.9,initsegci=TRUE,nsamples=0,gen.data=FALSE,prior=0.5,prior.type="n",verbose=TRUE,debug=FALSE) {
  if ((model!=1)&(model!=2)){
    stop("Choose model=1 (Poisson) or 2 (normal)")
  }
  n=length(data);
  J=length(seg)+1
  table.seg=table(seg)
  if (sum(table.seg>1)>0) {    
    stop(paste("Change-point",names(table.seg)[table.seg>1],"is repeated at least twice, please include only once"));
  }
  seg=sort(seg)
  if (length(lprob)>0 & length(data)>0){
    stop("Both data/segmentation/model and log-densities entered, specify only one");
  }
  if (length(lprob)==0){
  if (n<2) {
    stop("Not enough data to segment");
  }
  if (J>=n-1) {
    stop("Number of segments larger than data");
  }
  if (max(seg)>=n) {
    stop("At least one change-point initialized to larger than n");
  }
  
  if ((sum(data<0)>0)&model==1) {
    stop("Negative numbers with Poisson distribution specified, choose model=2");
  }
  if ((sum(round(data)!=data)>0)&model==1) {
    stop("Non-integer numbers with Poisson distribution specified, choose model=2");
  }
  } else{
  if (!is.matrix(lprob)) {
    stop("Enter lprob as a matrix, with n rows (one for each observation) and J columns (one for each hidden state)");
  }

    n=nrow(lprob);
    J=ncol(lprob);
    lprob=as.vector(lprob);
  }
  if (sum(prior<0|prior>1)>0) {
    stop("Choose all priors between 0 and 1");
  }

  if ((length(prior)!=n)&(prior.type=="o")) {
    stop("n transitions (1 for each observation) need to be chosen for prior vector");
  }
  if ((length(prior)!=J)&(prior.type=="s")) {
    stop("J transitions (1 for each state) need to be chosen for prior vector");
  }

 if (prior.type=="n") ptype=1 else
 if (prior.type=="o") ptype=2 else
 if (prior.type=="s") ptype=3 else
   stop("if prior.type specified, must be by (o)bservation or (s)tate");
 if (!keep){
   lforward=0
   lbackward=0
   post.cp=0
 } else{
   lforward=rep(0,n*J)
   lbackward=rep(0,n*J)
   post.cp=rep(0,(n-1)*(J-1))
 } 
 rmu=rep(0,J)
 rsigma=0

 cpconfint=rep(0,(J-1)*3)

 cpvector=rep(0,nsamples*(J-1))
  if (length(lprob)==0){ 
 .C("postCP",data=as.double(data),rseg=as.integer(seg),nn=as.integer(n),JJ=as.integer(J),lforward=as.double(lforward),lbackward=as.double(lbackward),cp=as.double(post.cp),rmu=as.double(rmu),rsigma=as.double(rsigma),cpconfint=as.double(cpconfint),cpvector=as.double(cpvector),rci=as.double(ci),rinitsegci=as.integer(initsegci),rnsamples=as.integer(nsamples),rmodel=as.integer(model),rprior=as.double(prior),rpriortype=as.integer(ptype),rmout=as.integer(keep),rverbose=as.integer(verbose),rdebug=as.integer(debug),DUP=FALSE,PACKAGE="postCP")} else{
 .C("postCPext",lprob=as.double(lprob),nn=as.integer(n),JJ=as.integer(J),lforward=as.double(lforward),lbackward=as.double(lbackward),cp=as.double(post.cp),cpconfint=as.double(cpconfint),cpvector=as.double(cpvector),rci=as.double(ci),rnsamples=as.integer(nsamples),rprior=as.double(prior),rpriortype=as.integer(ptype),rmout=as.integer(keep),rverbose=as.integer(verbose),rdebug=as.integer(debug),DUP=FALSE,PACKAGE="postCP")
   }

 
   if (nsamples>0) {
     cpmatrix=matrix(cpvector,ncol=J-1)
     fbsample <- list(changepoints=cpmatrix)
     names(fbsample) <- "changepoints"
     if (gen.data){
       getSegment <- function(cpv,nsamples) rowSums(outer((1:nsamples),cpv,">"))+1
       segments <- c(apply(cpmatrix,1,getSegment,nsamples=n))
       if (model==1) x=matrix(rpois(length(segments),rmu[segments]),ncol=n,byrow=T)
       if (model==2) x=matrix(rnorm(length(segments),rmu[segments],rsigma),ncol=n,byrow=T)
       fbsample$x=x
       names(fbsample)[2] <- "data"
     }
     cp.out=matrix(cpconfint,ncol=3)
     colnames(cp.out)=c("est",paste("lo.",ci,sep=""),paste("hi.",ci,sep=""))
     if (model==1) {
        model.dist="Poisson"
        postCP.res=list(model=model.dist,n=n,cp.est=cp.out,means=rmu,prior=prior,prior.type=ptype,fbsample=fbsample)
     }
     if (model==2) {
        model.dist="normal"
        postCP.res=list(model=model.dist,n=n,cp.est=cp.out,means=rmu,sds=rsigma,prior=prior,prior.type=ptype,fbsample=fbsample)
     }
   } else{ 
     cp.out=matrix(cpconfint,ncol=3)
     colnames(cp.out)=c("est",paste("lo.",ci,sep=""),paste("hi.",ci,sep=""))
     if (model==1)  {
       model.dist="Poisson"
       postCP.res=list(model=model.dist,n=n,cp.est=cp.out,means=rmu,prior=prior,prior.type=ptype)
     }
     if (model==2) {
       model.dist="normal"
      postCP.res=list(model=model.dist,n=n,cp.est=cp.out,means=rmu,sds=rsigma,prior=prior,prior.type=ptype)
     }     
  } 
  if (keep){
     if (length(lprob)==0) postCP.res$data=data else{
     postCP.res$lprob=matrix(lprob,ncol=J)
     }
     postCP.res$lforward=matrix(lforward,ncol=J)
     postCP.res$lbackward=matrix(lbackward,ncol=J)
     postCP.res$post.cp=matrix(post.cp,ncol=J-1)
     postCP.res$post.state=exp(postCP.res$lforward+postCP.res$lbackward-postCP.res$lforward[1,1]-postCP.res$lbackward[1,1])
  }
  if (length(lprob)>0) postCP.res$model="blank"
  postCP.res$call <- match.call()
  class(postCP.res) <- "postCP"
  postCP.res
}

print.postCP <- function(x,...)
{ cat("Call:\n")
  print(x$call)
  data.ent=x$model!="blank" # T if data entered instead of log-densities
  if (data.ent) out.table=data.frame(seg=1:length(x$means),start=c(1,x$cp.est[,1]+1),end=c(x$cp.est[,1],x$n),size=c(x$cp.est[,1],x$n)-c(1,x$cp.est[,1]+1)+1,mean=x$means) else out.table=data.frame(seg=1:length(x$means),start=c(1,x$cp.est[,1]+1),end=c(x$cp.est[,1],x$n),size=c(x$cp.est[,1],x$n)-c(1,x$cp.est[,1]+1)+1)

  if (data.ent) {cat("\nDistribution:")
  print(x$model)}
  cat("\nSegment information:\n")
  print(out.table)
  if (!is.element("fbsample",names(x))) print(x[!is.element(names(x),c("call","data","lforward","lbackward","lprob","means","post.cp","means","post.state","model","prior","prior.type"))])
  else {
    cat("\nEstimated change points (with ci):\n")
    print(x[!is.element(names(x),c("fbsample","data","lprob","means","lforward","lbackward","post.cp","post.state","model","call","prior","prior.type"))])
    cat("Generated change-points: ($fbsample$changepoints)\n")
    print(x$fbsample[[1]])
    if (length(x$fbsample)==2){
      cat("Generated data: ($fbsample$data)\n")
      str(x$fbsample[[2]])
    }
  }
  if (is.element("post.cp",names(x))){
    if (data.ent){
      cat("\nData: ($data) ")
      str(x$data)
    } else{
      cat("\nLog-densities ($lprob) ")
      str(x$lprob)
    }    
    cat("\nBackward matrix: ($lforward)")
    str(x$lbackward)
    cat("\nForward matrix: ($lbackward)")
    str(x$lforward)
    cat("\nChange-point matrix: ($post.cp)")
    str(x$post.cp)
    cat("\nState probability matrix: ($post.state)")
    str(x$post.state)
  }
  cat("\nPrior:\n")
  str(x$prior)
  cat("\nPrior type: ")
  if (x$prior.type==1) cat("homogeneous HMM\n")
  if (x$prior.type==2) cat("by observation\n")
  if (x$prior.type==3) cat("by hidden state\n")
}

