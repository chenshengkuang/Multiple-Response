multi_outcome_sgid=function(Y,X,Trt,P,resp_type,maxiter,mu,lossf=c("logistic","square"),base=c("spline","linear","stump")){
  n=dim(X)[1];p=dim(X)[2];q=dim(Y)[2]
  
  IntX=matrix(0,nrow=n,ncol=p)
  for (i in 1:n)
    IntX[i,]=X[i,]*(2*Trt[i]-1)
  CbX=cbind(X,(2*Trt-1),IntX) 
  #AIPWE
  shiftY=matrix(0,nrow=n,ncol=q)
  eff_notrt=matrix(0,nrow=n,ncol=q)
  eff_trt=matrix(0,nrow=n,ncol=q)
  CbX_notrt=cbind(X,rep(-1,n),-X)
  CbX_trt=cbind(X,rep(1,n),X) 

  for (j in 1:q){
    if (resp_type[j]=="binary"){
      lasmod=cv.glmnet(y=Y[,j],x=CbX,family="binomial")
      eff_notrt[,j]=predict(lasmod,newx=CbX_notrt,s="lambda.min",type="response")
      eff_trt[,j]=predict(lasmod,newx=CbX_trt,s="lambda.min",type="response")
    }
    if (resp_type[j]=="continuous"){
      lasmod=cv.glmnet(y=Y[,j],x=CbX)
      eff_notrt[,j]=predict(lasmod,newx=CbX_notrt,s="lambda.min")
      eff_trt[,j]=predict(lasmod,newx=CbX_trt,s="lambda.min")
    }
    shiftY[,j]=Y[,j]-(1-P)*eff_trt[,j]-P*eff_notrt[,j]
  }

  Con=matrix(0,nrow=n,ncol=q) #Contrast function
  for (i in 1:n)
    for (j in 1:q)
      Con[i,j]=Trt[i]*shiftY[i,j]/P-(1-Trt[i])*shiftY[i,j]/(1-P)
  
  sign_Con=matrix(0,nrow=n,ncol=q) #Sign of contrast function
  sign_Con[Con>0]=1
  W=abs(Con)
  multiboost_mod=Multiboost(Y=sign_Con,W=W,X=X,maxiter=maxiter,mu=mu,lossf=lossf,base=base)
  
  return(multiboost_mod)
}