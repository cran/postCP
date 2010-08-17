#define INF 1e300
//#include "postCP.h"

#include <R.h>
#include <Rmath.h>
#include <stddef.h>
#include <vector>

#include <iostream>
#include <iomanip>
#include <fstream>
#include <string>
#include <vector>
#include <cstdlib>
#include <cmath>
#include <ctime>

#define INF 1e300


extern "C"{
using std::vector;
using std::string;
using std::cout;
using std::setprecision;
using std::endl;

// #include <R.h>

using namespace std;

//vector<double> data;
vector<int> seg;
vector<double> mu;


int n;
int J;
int mout;

int normal;
double sigma;
double ci;

int priortype;
int nsamples;
int verbose,debug;
double priortemp;
double diffclock(clock_t clock1,clock_t clock2) {
	double diffticks=clock1-clock2;
	double diffs=(diffticks)/CLOCKS_PER_SEC;
	return diffs;
}

// returns log(exp(la)+exp(lb))
double lsum(double la, double lb) {
  if (la>lb){
    if (lb>(-INF)) return la+log1p(exp(lb-la));
    else return la;
    }
  else{
    if (la>(-INF)) return lb+log1p(exp(la-lb));    
    else return lb;
  }
}
double lfactorial(int x) {
double lfact=0;
if ((x==0)|(x==1)) lfact=0;
else {
  if (x<=23) {
    for (int i=1; i<=x;i++){
      lfact+=log(i);
     }
  }
 else lfact=0.5*log(2*x+1/3)+0.5*log(M_PI)+x*log(x)-x;

}
return(lfact);
}

void forward(double* data, double* lforward, double rprior[]) {
  vector<double> ldist(J);

  double aux1=-0.5*log(2*M_PI);
  double aux2=-0.5;
  //double aux3=log(0.5);
 
 // bool lower;
if (debug) Rprintf("normal? %i\n",normal);
  if (normal) lforward[0]=aux1+aux2*(data[0]-mu[0])*(data[0]-mu[0])/sigma/sigma;
  else lforward[0]=data[0]*log(mu[0])-mu[0]-lfactorial(data[0]);
  for (int j=1; j<J; j++) lforward[n*j]=-INF;
  if (verbose) Rprintf("\nCalculating forward matrix...\n");     
//  int a=0;
  for (int i=1; i<n; i++) {
    if (verbose) {
      if ((i+1)%10000==0)   Rprintf(".");
    }
    // precompute ldist for data[i]
    for (int j=0; j<J; j++) {
      if (normal) ldist[j]=aux1+aux2*(data[i]-mu[j])*(data[i]-mu[j])/sigma/sigma;
      else ldist[j]=data[i]*log(mu[j])-mu[j]-lfactorial(data[i]);
    }
    // recursion step
    if (priortype==2) priortemp=rprior[i];
    if (priortype==3) priortemp=rprior[0];

    lforward[i]=lforward[i-1]+log(1-priortemp)+ldist[0];

    for (int j=1; j<J; j++) {
    if (priortype==2) priortemp=rprior[i];
    if (priortype==3) priortemp=rprior[j];

      lforward[n*j+i]=lsum(lforward[n*(j-1)+i-1]+log(priortemp),lforward[(n*j)+i-1]+log(1-priortemp))+ldist[j];
    }
  }  
};

void backward(double* data, double* lbackward, double* lforward, double* cp, double rprior[]) {
  vector<double> ldist(J);

  double aux1=-0.5*log(2*M_PI);
  double aux2=-0.5;

 // double aux3=log(0.5);

 // bool upper;    
  if (debug) Rprintf("normal? %i\n",normal);
//  int a=J-2;
  lbackward[n*J-1]=0.0;
  for (int j=0; j<J-1; j++) lbackward[n*j+n-1]=-INF;
  if (verbose) Rprintf("\nCalculating backward matrix...\n");      
  for (int i=n-2; i>=0; i--) {
   if (verbose) {
      if ((i+1)%10000==0)  Rprintf("."); 
    }
    // precompute ldist for data[i+1]
    for (int j=0; j<J; j++) {
      if (normal) ldist[j]=aux1+aux2*(data[i+1]-mu[j])*(data[i+1]-mu[j])/sigma/sigma;
      else ldist[j]=data[i+1]*log(mu[j])-mu[j]-lfactorial(data[i+1]);
    }
    // recursion step
    if (priortype==2) priortemp=rprior[i];
    if (priortype==3) priortemp=rprior[J-1];

     lbackward[n*(J-1)+i]=lbackward[n*(J-1)+i+1]+log(1-priortemp)+ldist[J-1];
    for (int j=J-2; j>=0; j--) {
    if (priortype==2) priortemp=rprior[i];
    if (priortype==3) priortemp=rprior[j];

    lbackward[n*j+i]=lsum(lbackward[n*j+i+1]+ldist[j]+log(1-priortemp),lbackward[n*(j+1)+i+1]+ldist[j+1]+log(priortemp)); 
    cp[(n-1)*j+i]=exp(lforward[n*j+i]+lbackward[n*(j+1)+i+1]+ldist[j+1]+log(priortemp)-lforward[n*J-1]-lbackward[n*J-1]);
    }
  }  
};


void postCPsample(double* data, double* rmu, double *rsigma, double* rlbackward, double* cpvector, int *rnsamples, int *nn, int *JJ, int *rmodel, double rprior[], int *rpriortype, int *rverbose, int *rdebug) {
int nsamples=*rnsamples;

int n=*nn;
int J=*JJ;
int verbose=*rverbose;
int debug=*rdebug;
int priortype=*rpriortype;

normal=0;
if (*rmodel==2) normal=1;
double sigma=*rsigma;
 GetRNGstate(); 

    if (priortype==1) priortemp=rprior[0];

for (int r=0; r<nsamples; r++) {
  // int cp=0;
   int i=0;
  
   double rand;
   double aux;
  
   double aux1=-0.5*log(2*M_PI);
   double aux2=-0.5;
//   double aux3=log(0.5);
  
   if (verbose&((r+1)%500==0)) Rprintf("\nGenerating replicate %i\n",r+1);
  
   for (int j=0;j<J-1;j++){
      while (true) {
      if (debug)  Rprintf("i:%i,j:%i,n:%i,J:%i",i,j,n,J);   
      if (priortype==2) priortemp=rprior[i];
      if (priortype==3) priortemp=rprior[j];
       if (normal) aux=rlbackward[n*(j+1)+i+1]-rlbackward[n*j+i]+log(priortemp)+aux1+aux2*(data[i+1]-rmu[j+1])*(data[i+1]-rmu[j+1])/sigma;
       else aux=rlbackward[n*(j+1)+i+1]-rlbackward[n*j+i]+log(priortemp)+data[i+1]*log(rmu[j+1])-rmu[j+1]-lfactorial(data[i+1]);
       rand=log(runif(0,1.0));
       i++;
       if (rand<aux)
     break;
     }
     cpvector[nsamples*j+r]=i;
   }
 }; // end r loop
PutRNGstate();
//   return;
};


void changepointci(double* cpconfint, double* cp, int initsegci)
  { // obtain confidence interval of changepoint
    // start with cp with highest changepoint post prob, add adjacent points, whichever has higher post prob until total post prob within range is higher than entered ci

    double conf=0; // temporary cumulative confidence interval
    int ind=0,ind1,ind2; // ind: obs with highest prob of change-point, ind1: lower bound, ind2: higher bound
    
    for (int j=0;j<J-1; j++) {
      double bestll=0;  
      if(initsegci){
          cpconfint[j]=seg[j+1];
          ind=seg[j+1]-1;
      }
      else{
          for (int i=0;i<n-1;i++) if (cp[j*(n-1)+i]>bestll) {
            ind=i;
   	    bestll=cp[j*(n-1)+i]; 
             if (debug)      Rprintf("%i,%.8f..",ind,bestll);
          }
      cpconfint[j]=ind+1;
       if (debug)    Rprintf("%i\n",cpconfint[j]);
      }
      ind1=ind;
      ind2=ind;
      conf=cp[j*(n-1)+ind];

      while (conf<ci) {
	if ((ind1==1)|(cp[j*(n-1)+ind2+1]>cp[j*(n-1)+ind1-1])){
	  // increase upper bound 
	  ind2+=1;
	  conf+=cp[j*(n-1)+ind2];
	} else if ((ind2==n-1)|(cp[j*(n-1)+ind1-1]>cp[j*(n-1)+ind2+1])) {
	  // increase lower bound 
	  ind1-=1;
	  conf+=cp[j*(n-1)+ind1];
	} else if (cp[j*(n-1)+ind1-1]==cp[j*(n-1)+ind2+1]) {
	  // increase both if tied
          ind1-=1;
          ind2+=1;
          conf+=cp[j*(n-1)+ind1]+cp[j*(n-1)+ind2];
	}
      }
      cpconfint[J-1+j]=ind1+1;
      cpconfint[2*(J-1)+j]=ind2+1;
   //   cout<<cpconfint[J-1+j]<<" "<<cpconfint[2*(J-1)+j]<<endl;
    }
  }


void postCP (double* data, int* rseg, int *nn, int *JJ, double* lforward, double* lbackward, double* cp, double* rmu, double *rsigma, double* cpconfint, double* cpvector, double *rci, int *rinitsegci, int *rnsamples, int *rmodel, double rprior[], int *rpriortype, int *rmout, int *rverbose, int *rdebug) {

int initsegci=*rinitsegci;
verbose=*rverbose;
debug=*rdebug;
mout=*rmout;
n=*nn; 
J=*JJ;

ci=*rci;
priortype=*rpriortype;
nsamples=*rnsamples;
normal=0;
if (priortype==1) priortemp=rprior[0];
if (*rmodel==2) normal=1;


if (debug)  Rprintf("normal? %i verbose %i\n",verbose);

 if (debug) Rprintf("n=%i,J=%i\n",n,J);
//data.resize(n);
//for (int i=0;i<n;i++){
//  data[i]=data[i];
//}
seg.resize(J);
seg[0]=0;
for (int j=1;j<J;j++){
  seg[j]=rseg[j-1];
 if (debug) Rprintf("%i %i\n",j,seg[j]);
}

seg.push_back(n);

  
//  J=seg.size()-1;
//  n=data.size();
 
  // compute the parameters from seg
  {
    mu.resize(J);
    for (int j=0; j<J; j++) {
      double sum1=0.0,sum0=0.0;
      if (debug) Rprintf("%i, seg=%i, seg[j+1]=%i\n",j,seg[j],seg[j+1]);
//      if (debug) cout<<j<<"seg="<<seg[j]<<"\tseg[j+1]="<<seg[j+1]<<endl;
      for (int i=seg[j]; i<seg[j+1]; i++) {
        sum1+=data[i];
 	sum0+=1.0;
     }
    mu[j]=sum1/sum0;
    rmu[j]=sum1/sum0;
     if (debug)  Rprintf("%i %.8f\n",j,rmu[j]);
 }

    double sum2=0.0,sum0=0.0;
   for (int j=0; j<J; j++) {
     for (int i=seg[j]; i<seg[j+1]; i++) {
	sum2+=(data[i]-mu[j])*(data[i]-mu[j]);
	sum0+=1.0;
      }
    }
    sigma=sqrt(sum2/sum0);
   *rsigma=sqrt(sum2/sum0);
     if (debug) Rprintf("%.8f\n",*rsigma); 
  }

  // print the parameters
  if (debug) Rprintf("sigma=%.8f\n",sigma);
  for (int j=0; j<J; j++)
    if (debug) Rprintf("mu[%i]=%.8f\n",j,mu[j]); 
  if (debug) Rprintf("\n"); 


  clock_t begin=clock();
  
  if (!mout) { 
    if (debug) Rprintf("creating matrices");
    double* lforward1=NULL;
    lforward1=new double[n*J];
    double* lbackward1=NULL;
    lbackward1=new double[n*J];
    if (debug) Rprintf("%i x %i\n",n,J); 
    double* cp1=NULL;
    cp1=new double[(n-1)*(J-1)];
    if (debug) Rprintf("matrices created");
    forward(data,lforward1,rprior);
    backward(data,lbackward1,lforward1,cp1,rprior);
    if (ci>0) changepointci(cpconfint,cp1,initsegci);

    if (verbose) Rprintf("\n log(pevidence)=%.8f=%.8f\n",lforward1[0]+lbackward1[0],lforward1[n*J-1]+lbackward1[n*J-1]);

   //cout<<setprecision(12)<<"\nlog(pevidence)="<<lforward1[0]+lbackward1[0]<<"="<<lforward1[n*J-1]+lbackward1[n*J-1]<<endl;

    if (nsamples>0) postCPsample(data, rmu, rsigma, lbackward1, cpvector, rnsamples, nn, JJ, rmodel,rprior,rpriortype,rverbose,rdebug);
    delete [] lforward1;
    delete [] lbackward1;
    delete [] cp1;
    lforward1=NULL;
    lbackward1=NULL;
    cp1=NULL;
  }
  else {
    forward(data,lforward,rprior);
    backward(data,lbackward,lforward,cp,rprior);
    if (ci>0) changepointci(cpconfint,cp,initsegci);

   if (verbose) Rprintf("\n log(pevidence)=%.8f=%.8f\n",lforward[0]+lbackward[0],lforward[n*J-1]+lbackward[n*J-1]);

//cout<<setprecision(12)<<"\nlog(pevidence)="<<lforward1[0]+lbackward1[0]<<"="<<lforward1[n*J-1]+lbackward1[n*J-1]<<endl;
    if (nsamples>0) postCPsample(data, rmu, rsigma, lbackward, cpvector, rnsamples, nn, JJ, rmodel,rprior,rpriortype,rverbose,rdebug);
  }

  //forward2();
  //backward2();
  
  clock_t end=clock();
   if (debug)    Rprintf("Time elapsed: %.3f s\n", double(diffclock(end,begin)));




//  return 0;
}

void forwardext(double* lprob, double* lforward, double rprior[]) {
//void forward(double* data, double* lforward, double rprior[]) {
  lforward[0]=lprob[0];

  for (int j=1; j<J; j++) lforward[n*j]=-INF;
  if (verbose) Rprintf("\nCalculating forward matrix...\n");      
  //  int a=0;
  for (int i=1; i<n; i++) {
    if (verbose) {
      if ((i+1)%10000==0)   Rprintf("."); 
    }
    // recursion step
    if (priortype==2) priortemp=rprior[i];
    if (priortype==3) priortemp=rprior[0];

    lforward[i]=lforward[i-1]+log(1-priortemp)+lprob[i];

    for (int j=1; j<J; j++) {
      if (priortype==2) priortemp=rprior[i];
      if (priortype==3) priortemp=rprior[j];
      lforward[n*j+i]=lsum(lforward[n*(j-1)+i-1]+log(priortemp),lforward[(n*j)+i-1]+log(1-priortemp))+lprob[n*j+i];
    }
  }  
};

void backwardext(double* lprob, double* lbackward, double* lforward, double* cp, double rprior[]) {

  lbackward[n*J-1]=0.0;
  for (int j=0; j<J-1; j++) lbackward[n*j+n-1]=-INF;
  if (verbose) Rprintf("\nCalculating backward matrix...\n");      
  for (int i=n-2; i>=0; i--) {
   if (verbose) {
      if ((i+1)%10000==0)   Rprintf("."); 
    }
    // precompute ldist for data[i+1]
    // recursion step
    if (priortype==2) priortemp=rprior[i];
    if (priortype==3) priortemp=rprior[J-1];

    lbackward[n*(J-1)+i]=lbackward[n*(J-1)+i+1]+log(1-priortemp)+lprob[n*(J-1)+i+1];

    for (int j=J-2; j>=0; j--) {
    if (priortype==2) priortemp=rprior[i];
    if (priortype==3) priortemp=rprior[j];

    lbackward[n*j+i]=lsum(lbackward[n*j+i+1]+lprob[n*j+i+1]+log(1-priortemp),lbackward[n*(j+1)+i+1]+lprob[n*(j+1)+i+1]+log(priortemp)); 
    cp[(n-1)*j+i]=exp(lforward[n*j+i]+lbackward[n*(j+1)+i+1]+lprob[n*(j+1)+i+1]+log(priortemp)-lforward[n*J-1]-lbackward[n*J-1]);
    }
  }  
};



void postCPsampleext(double* lprob, double* rlbackward, double* cpvector, int *rnsamples, int *nn, int *JJ, double rprior[], int *rpriortype, int *rverbose, int *rdebug) {
int nsamples=*rnsamples;

int n=*nn;
int J=*JJ;
int verbose=*rverbose;
int debug=*rdebug;
int priortype=*rpriortype;

 GetRNGstate(); 

   if (priortype==1) priortemp=rprior[0];

   for (int r=0; r<nsamples; r++) {
  // int cp=0;
   int i=0;
  
   double rand;
   double aux;
  

//   double aux3=log(0.5);
  
   if (verbose&((r+1)%500==0)) Rprintf("\nGenerating replicate %i\n",r+1);

   for (int j=0;j<J-1;j++){
     while (true) {
      if (debug)     Rprintf("i:%i,j:%i,n:%i,J:%i",i,j,n,J);   
      if (priortype==2) priortemp=rprior[i];
      if (priortype==3) priortemp=rprior[j];
       aux=rlbackward[n*(j+1)+i+1]-rlbackward[n*j+i]+log(priortemp)+lprob[n*(j+1)+i+1];       
   rand=log(runif(0,1.0));
       i++;
       if (rand<aux)
     break;
     }
     cpvector[nsamples*j+r]=i;
   }
 }; // end r loop
PutRNGstate();
};


void postCPext (double* lprob, int *nn, int *JJ, double* lforward, double* lbackward, double* cp, double* cpconfint, double* cpvector, double *rci, int *rnsamples, double rprior[], int *rpriortype, int *rmout, int *rverbose, int *rdebug) {

  verbose=*rverbose;
  debug=*rdebug;
  mout=*rmout;
  n=*nn;
  J=*JJ;

  ci=*rci;
  priortype=*rpriortype;
  nsamples=*rnsamples;
  if (priortype==1) priortemp=rprior[0];


if (debug)  Rprintf("normal? %i verbose %i\n",verbose);

 if (debug) Rprintf("n=%i,J=%i\n",n,J);

  clock_t begin=clock();
  
  if (!mout) { 
    if (debug) Rprintf("creating matrices");
    double* lforward1=NULL;
    lforward1=new double[n*J];
    double* lbackward1=NULL;
    lbackward1=new double[n*J];
    if (debug) Rprintf("%i x %i\n",n,J); 
    double* cp1=NULL;
    cp1=new double[(n-1)*(J-1)];
    if (debug) Rprintf("matrices created");
    forwardext(lprob,lforward1,rprior);
    backwardext(lprob,lbackward1,lforward1,cp1,rprior);
    if (ci>0) changepointci(cpconfint,cp1,0);
     if (verbose) Rprintf("\n log(pevidence)=%.8f=%.8f\n",lforward1[0]+lbackward1[0],lforward1[n*J-1]+lbackward1[n*J-1]);

//cout<<setprecision(12)<<"\nlog(pevidence)="<<lforward1[0]+lbackward1[0]<<"="<<lforward1[n*J-1]+lbackward1[n*J-1]<<endl;

    if (nsamples>0) postCPsampleext(lprob,  lbackward1, cpvector, rnsamples, nn, JJ, rprior,rpriortype,rverbose,rdebug);
    delete [] lforward1;
    delete [] lbackward1;
    delete [] cp1;
    lforward1=NULL;
    lbackward1=NULL;
    cp1=NULL;
  }
  else {
    forwardext(lprob,lforward,rprior);
    backwardext(lprob,lbackward,lforward,cp,rprior);
    if (ci>0) changepointci(cpconfint,cp,0);

     if (verbose) Rprintf("\n log(pevidence)=%.8f=%.8f\n",lforward[0]+lbackward[0],lforward[n*J-1]+lbackward[n*J-1]);

//cout<<setprecision(12)<<"\nlog(pevidence)="<<lforward1[0]+lbackward1[0]<<"="<<lforward1[n*J-1]+lbackward1[n*J-1]<<endl;
    if (nsamples>0) postCPsampleext(lprob, lbackward, cpvector, rnsamples, nn, JJ, rprior,rpriortype,rverbose,rdebug);
  }

  //forward2();
  //backward2();
  
  clock_t end=clock();

    if (debug)    Rprintf("Time elapsed: %.3f s\n", double(diffclock(end,begin)));



//  return 0;
};


}


