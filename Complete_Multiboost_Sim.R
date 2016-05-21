complete_multiple_response_sim=function(resp_type,non_linear,partial_effect,n,p,sigma){ 
  q=6 #p: # of covariates, q:# of responses, n: # of obs
  set.seed(2)
  W1=matrix(sign(rnorm(q*(p+1),0,1)),nrow=p+1,ncol=q)
  Q1=matrix(rep(c(1,0),c(11*q,(p-10)*q)),nrow=p+1,byrow=TRUE)
  B1=W1*Q1
  W2=matrix(2*sign(rnorm(q*(p+1),0,1)),nrow=p+1,ncol=q)
  Q2=matrix(rep(c(1,0),c(6*q,(p-5)*q)),nrow=p+1,byrow=TRUE)
  B2=W2*Q2
  B=rbind(B1,B2)
  if (partial_effect==1){
    for (i in (p+3):(p+7)){
      B[i,sample(1:q,3,replace=FALSE)]=0
    }
  }
  set.seed(NULL) #unset seed to allow random data
  P=0.5 #propensity score
  
  B[1,]=B[1,]/2;B[p+2,]=B[p+2,]/2
  
  ###########################
  #Generate Training Set
  #########################
  X=matrix(rnorm(n*p,0,1),nrow=n)
  Trt=rbinom(n,1,P) #Assignment
  IntX=matrix(0,nrow=n,ncol=p)
  for (i in 1:n)
    IntX[i,]=X[i,]*Trt[i]
  CbX=cbind(X,Trt,IntX) #Combined design matrix including interaction term
  
  Y=matrix(0,nrow=n,ncol=q)
  
  #nonlienar setting put some nonlienar term into main effect but keep main effect not too big
  if (non_linear==1){
    if (resp_type=="binary"){
      for (i in 1:n)
        for (j in 1:q)
          Y[i,j]=rbinom(1,1,1/(1+exp(-0.4*(B[1,j]+CbX[i,1:p]%*%B[2:(p+1),j]+CbX[i,1]*CbX[i,2]+CbX[i,3]*CbX[i,4]+CbX[i,(p+1):(2*p+1)]%*%B[(p+2):(2*p+2),j]+2*CbX[i,p+2]^2))))
    }
    if (resp_type=="continuous"){
      for (i in 1:n)
        for (j in 1:q)
          Y[i,j]=B[1,j]+CbX[i,1:p]%*%B[2:(p+1),j]+CbX[i,1]*CbX[i,2]+CbX[i,3]*CbX[i,4]+CbX[i,(p+1):(2*p+1)]%*%B[(p+2):(2*p+2),j]+2*CbX[i,p+2]^2+rnorm(1,0,sigma)
    }
  }
  
  #linear 
  if (non_linear==0){
    if (resp_type=="binary"){
      for (i in 1:n)
        for (j in 1:q)
          Y[i,j]=rbinom(1,1,1/(1+exp(-0.4*(B[1,j]+CbX[i,]%*%B[2:(2*p+2),j]))))
    }
    if (resp_type=="continuous"){
      for (i in 1:n)
        for (j in 1:q)
          Y[i,j]=B[1,j]+CbX[i,]%*%B[2:(2*p+2),j]+rnorm(1,0,sigma)
    }
  }
  
  resp_IntX=matrix(0,nrow=n,ncol=p)
  for (i in 1:n)
    resp_IntX[i,]=X[i,]*(2*Trt[i]-1)
  resp_CbX=cbind(X,(2*Trt-1),resp_IntX) 
  
  #AIPWE
  shiftY=matrix(0,nrow=n,ncol=q)
  eff_notrt=matrix(0,nrow=n,ncol=q)
  eff_trt=matrix(0,nrow=n,ncol=q)
  CbX_notrt=cbind(X,rep(-1,n),-X)
  CbX_trt=cbind(X,rep(1,n),X) 
  
  for (j in 1:q){
    if (resp_type=="binary"){
      lasmod=cv.glmnet(y=Y[,j],x=resp_CbX,family="binomial")
      eff_notrt[,j]=predict(lasmod,newx=CbX_notrt,s="lambda.min",type="response")
      eff_trt[,j]=predict(lasmod,newx=CbX_trt,s="lambda.min",type="response")
    }
    if (resp_type=="continuous"){
      lasmod=cv.glmnet(y=Y[,j],x=resp_CbX)
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
  
  ###########################
  #Generate validation Set
  #########################
  vali_X=matrix(rnorm(n*p,0,1),nrow=n)
  vali_Trt=rbinom(n,1,P) #Assignment
  vali_IntX=matrix(0,nrow=n,ncol=p)
  for (i in 1:n)
    vali_IntX[i,]=vali_X[i,]*vali_Trt[i]
  vali_CbX=cbind(vali_X,vali_Trt,vali_IntX) #Combined design matrix including interaction term
  
  vali_Y=matrix(0,nrow=n,ncol=q)
  
  #nonlienar setting put some nonlienar term into main effect but keep main effect not too big
  if (non_linear==1){
    if (resp_type=="binary"){
      for (i in 1:n)
        for (j in 1:6)
          vali_Y[i,j]=rbinom(1,1,1/(1+exp(-0.4*(B[1,j]+vali_CbX[i,1:p]%*%B[2:(p+1),j]+vali_CbX[i,1]*vali_CbX[i,2]+
                                                  vali_CbX[i,3]*vali_CbX[i,4]+vali_CbX[i,(p+1):(2*p+1)]%*%B[(p+2):(2*p+2),j]+
                                                  2*vali_CbX[i,p+2]^2))))
    }
    if (resp_type=="continuous"){
      for (i in 1:n)
        for (j in 1:6)
          vali_Y[i,j]=B[1,j]+vali_CbX[i,1:p]%*%B[2:(p+1),j]+vali_CbX[i,1]*vali_CbX[i,2]+
            vali_CbX[i,3]*vali_CbX[i,4]+vali_CbX[i,(p+1):(2*p+1)]%*%B[(p+2):(2*p+2),j]+2*vali_CbX[i,p+2]^2+rnorm(1,0,sigma)
    }
  }
  
  #linear 
  if (non_linear==0){
    if (resp_type=="binary"){
      for (i in 1:n)
        for (j in 1:6)
          vali_Y[i,j]=rbinom(1,1,1/(1+exp(-0.4*(B[1,j]+vali_CbX[i,]%*%B[2:(2*p+2),j]))))
    }
    if (resp_type=="continuous"){
      for (i in 1:n)
        for (j in 1:6)
          vali_Y[i,j]=B[1,j]+vali_CbX[i,]%*%B[2:(2*p+2),j]+rnorm(1,0,sigma)
    }
  }
  
  vali_resp_IntX=matrix(0,nrow=n,ncol=p)
  for (i in 1:n)
    vali_resp_IntX[i,]=vali_X[i,]*(2*vali_Trt[i]-1)
  vali_resp_CbX=cbind(vali_X,(2*vali_Trt-1),vali_resp_IntX) 
  
  #AIPWE
  vali_shiftY=matrix(0,nrow=n,ncol=q)
  vali_eff_notrt=matrix(0,nrow=n,ncol=q)
  vali_eff_trt=matrix(0,nrow=n,ncol=q)
  vali_CbX_notrt=cbind(vali_X,rep(-1,n),-vali_X)
  vali_CbX_trt=cbind(vali_X,rep(1,n),vali_X) 
  
  for (j in 1:q){
    if (resp_type=="binary"){
      vali_lasmod=cv.glmnet(y=vali_Y[,j],x=vali_resp_CbX,family="binomial")
      vali_eff_notrt[,j]=predict(vali_lasmod,newx=vali_CbX_notrt,s="lambda.min",type="response")
      vali_eff_trt[,j]=predict(vali_lasmod,newx=vali_CbX_trt,s="lambda.min",type="response")
    }
    if (resp_type=="continuous"){
      vali_lasmod=cv.glmnet(y=vali_Y[,j],x=vali_resp_CbX)
      vali_eff_notrt[,j]=predict(vali_lasmod,newx=vali_CbX_notrt,s="lambda.min")
      vali_eff_trt[,j]=predict(vali_lasmod,newx=vali_CbX_trt,s="lambda.min")
    }
    vali_shiftY[,j]=vali_Y[,j]-(1-P)*vali_eff_trt[,j]-P*vali_eff_notrt[,j]
  }
  
  vali_Con=matrix(0,nrow=n,ncol=q) #Contrast function
  for (i in 1:n)
    for (j in 1:q)
      vali_Con[i,j]=vali_Trt[i]*vali_shiftY[i,j]/P-(1-vali_Trt[i])*vali_shiftY[i,j]/(1-P)
  
  vali_sign_Con=matrix(0,nrow=n,ncol=q) #Sign of contrast function
  vali_sign_Con[vali_Con>0]=1
  vali_W=abs(vali_Con)
  
  ######################################
  #Generate Test Dataset
  ######################################
  ntest=100000
  Xtest=matrix(rnorm(ntest*p,0,1),nrow=ntest)
  Xtest=cbind(rep(1,ntest),Xtest) #Add Intercept
  
  if (non_linear==1 ){
    opt_dec=sign(Xtest%*%B[(p+2):(2*p+2),]+2*Xtest[,2]^2)
    opt_value=Xtest%*%B[(p+2):(2*p+2),]+2*Xtest[,2]^2
  }
  if (non_linear==0 ){
    opt_dec=sign(Xtest%*%B[(p+2):(2*p+2),])
    opt_value=Xtest%*%B[(p+2):(2*p+2),]
  }
  
  ##############################################
  #Response separate Lasso
  #############################################
  resp_lasso_beta=matrix(0,nrow=2*p+2,ncol=q)
  resp_lasso_lambda=rep(0,q)
  resp_lasso_mod_list=rep(list(list()),q)
  for (j in 1:q){
    resp_lasso_mod_list[[j]]=glmnet(x=resp_CbX,y=Y[,j])
  }
  #Tune lambda
  for (j in 1:q){
    resp_lasso_predict=predict(resp_lasso_mod_list[[j]],newx=vali_resp_CbX)
    resp_lasso_loss=rep(0,length(resp_lasso_mod_list[[j]]$lambda))
    for (l in 1:length(resp_lasso_mod_list[[j]]$lambda)){
      resp_lasso_loss[l]=sum((vali_Y[,j]-resp_lasso_predict[,l])^2)
    }
    resp_lasso_lambda[j]=resp_lasso_mod_list[[j]]$lambda[which.min(resp_lasso_loss)]
    resp_lasso_beta[,j]=as.vector(predict(resp_lasso_mod_list[[j]],type="coefficients",s=resp_lasso_lambda[j]))
  }
  #Misclassification rate for response lasso
  resp_lasso_dec=sign(Xtest%*%resp_lasso_beta[(p+2):(2*p+2),])
  resp_lasso_rate=(colSums(opt_dec==resp_lasso_dec)+0.5*colSums(resp_lasso_dec==0))/ntest
  resp_lasso_value=Xtest%*%resp_lasso_beta[(p+2):(2*p+2),]
  resp_lasso_cor=rep(0,q)
  for (j in 1:q){
    resp_lasso_cor[j]=cor(resp_lasso_value[,j],opt_value[,j],method="spearman")
  }
  
  ################################
  #Lasso
  ################################
  lasso_beta=matrix(0,nrow=p+1,ncol=q)
  lasso_list=rep(list(list()),q)
  lasso_lambda=rep(0,q)
  for (j in 1:q){
    lasso_list[[j]]=glmnet(x=X,y=sign_Con[,j],weights=W[,j])
    lasso_predict=predict(lasso_list[[j]],newx=vali_X)
    lasso_loss=rep(0,length(lasso_list[[j]]$lambda))
    for (l in 1:length(lasso_list[[j]]$lambda)){
      lasso_loss[l]=sum(vali_W[,j]*(vali_sign_Con[,j]-lasso_predict[,l])^2)
    }
    lasso_lambda[j]=lasso_list[[j]]$lambda[which.min(lasso_loss)]
    lasso_beta[,j]=as.vector(predict(lasso_list[[j]],type="coefficients",s=lasso_lambda[j]))
  }
  #misclassification rate for contrast lasso
  lasso_dec=sign(Xtest%*%lasso_beta-0.5)
  lasso_rate=(colSums(opt_dec==lasso_dec)+0.5*colSums(lasso_dec==0))/ntest
  lasso_value=Xtest%*%lasso_beta-0.5
  lasso_cor=rep(0,q)
  for (j in 1:q){
    lasso_cor[j]=cor(lasso_value[,j],opt_value[,j],method="spearman")
  }
  
  ###############################
  #SGL
  ###############################
  SGL_y=numeric(n*q)
  SGL_X=matrix(0,nrow=n*q,ncol=p*q)
  index=rep(seq(1,p),q)
  #Center each covariate and response
  for (j in 1:q){
    SGL_y[((j-1)*n+1):(j*n)]=sign_Con[,j]-sum(W[,j]*sign_Con[,j])/sum(W[,j])
    for (k in 1:p){
      SGL_X[((j-1)*n+1):(j*n),(p*(j-1)+k)]=X[,k]-sum(W[,j]*X[,k])/sum(W[,j])
    }
  }
    
  #Times weight to each data point
  for (j in 1:q){
    Con_mean=sum(W[,j]*sign_Con[,j])/sum(W[,j])
    Con_ss=sum(W[,j]*(sign_Con[,j]-Con_mean)^2)
    SGL_y[((j-1)*n+1):(j*n)]=SGL_y[((j-1)*n+1):(j*n)]*sqrt(W[,j]/Con_ss)
    for (k in 1:p){
      SGL_X[((j-1)*n+1):(j*n),(p*(j-1)+k)]=SGL_X[((j-1)*n+1):(j*n),(p*(j-1)+k)]*sqrt(W[,j]/Con_ss)
    }
  }
    
  SGL_data=list(x=SGL_X,y=SGL_y)
  index=rep(c(1:p),q)
  fit=SGL(SGL_data,index,type="linear",standardize=FALSE,nlam=100)
    
  SGL_loss=rep(0,100)
  for (l in 1:100){
    SGL_beta=matrix(0,nrow=p+1,ncol=q)
    for (j in 1:q){
      SGL_beta[2:(p+1),j]=fit$beta[((j-1)*p+1):(j*p),l]
      SGL_beta[1,j]=sum(W[,j]*sign_Con[,j])/sum(W[,j])
      for (k in 1:p){
        SGL_beta[1,j]=SGL_beta[1,j]-SGL_beta[k+1,j]*sum(W[,j]*X[,k])/sum(W[,j])
      }
      vali_Con_mean=sum(vali_W[,j]*vali_sign_Con[,j])/sum(vali_W[,j])
      vali_Con_ss=sum(vali_W[,j]*(vali_sign_Con[,j]-vali_Con_mean)^2)
      SGL_loss[l]=SGL_loss[l]+sum(vali_W[,j]*(vali_sign_Con[,j]-vali_X%*%SGL_beta[2:(p+1),j]-SGL_beta[1,j])^2)/vali_Con_ss
    }
  }
    
  SGL_lambda=which.min(SGL_loss)
  SGL_beta=matrix(0,nrow=p+1,ncol=q)
  for (j in 1:q){
    SGL_beta[2:(p+1),j]=fit$beta[((j-1)*p+1):(j*p),SGL_lambda]
    SGL_beta[1,j]=sum(W[,j]*sign_Con[,j])/sum(W[,j])
    for (k in 1:p){
      SGL_beta[1,j]=SGL_beta[1,j]-SGL_beta[k+1,j]*sum(W[,j]*X[,k])/sum(W[,j])
    }
  }
  #Misclassification of contrast SGL
  SGL_dec=sign(Xtest%*%SGL_beta-0.5)
  SGL_rate=(colSums(opt_dec==SGL_dec)+0.5*colSums(SGL_dec==0))/ntest
  SGL_value=Xtest%*%SGL_beta-0.5
  SGL_cor=rep(0,q)
  for (j in 1:q){
    SGL_cor[j]=cor(SGL_value[,j],opt_value[,j],method="spearman")
  }
  
  
  
  ##############################################################
  #Boosting
  ##############################################################
  multiboost_rate=uniboost_rate=array(NA,dim=c(2,3,q))
  dimnames(multiboost_rate)[[1]]=c("square","logistic")
  dimnames(multiboost_rate)[[2]]=c("linear","spline","stump")
  dimnames(multiboost_rate)[[3]]=c("Outcome1","Outcome 2","Outcome 3","Outcome 4","Outcome 5","Outcome 6")
  dimnames(uniboost_rate)[[1]]=c("square","logistic")
  dimnames(uniboost_rate)[[2]]=c("linear","spline","stump")
  dimnames(uniboost_rate)[[3]]=c("Outcome1","Outcome 2","Outcome 3","Outcome 4","Outcome 5","Outcome 6")
  multiboost_cor=uniboost_cor=array(NA,dim=c(2,3,q))
  dimnames(multiboost_cor)[[1]]=c("square","logistic")
  dimnames(multiboost_cor)[[2]]=c("linear","spline","stump")
  dimnames(multiboost_cor)[[3]]=c("Outcome1","Outcome 2","Outcome 3","Outcome 4","Outcome 5","Outcome 6")
  dimnames(uniboost_cor)[[1]]=c("square","logistic")
  dimnames(uniboost_cor)[[2]]=c("linear","spline","stump")
  dimnames(uniboost_cor)[[3]]=c("Outcome1","Outcome 2","Outcome 3","Outcome 4","Outcome 5","Outcome 6")
  
  
  for (lossf_choice in 1:2)
    for (base_choice in 1:3){
      lossf=c("square","logistic")[lossf_choice]
      base=c("linear","spline","stump")[base_choice]
      
      #############################
      #Mboost
      #############################
      Mb_maxiter=2000;mu=0.1
      Mb_model=Multiboost(Y=sign_Con,W=W,X=X,maxiter=Mb_maxiter,mu=mu,lossf=lossf,base=base)
      Mb_loss_red_ctr=Mb_model[[5]]
      Mb_loss=rep(0,Mb_maxiter+1)
      current_fit=matrix(0,nrow=n,ncol=q)
      vali_Con_mean=numeric(q)
      vali_Con_ss=numeric(q)
      if (lossf=="square"){
        for (j in 1:q){
          vali_Con_mean[j]=sum(vali_W[,j]*vali_sign_Con[,j])/sum(vali_W[,j])
          vali_Con_ss[j]=sum(vali_W[,j]*(vali_sign_Con[,j]-vali_Con_mean[j])^2)
          current_fit[,j]=sum(sign_Con[,j]*W[,j])/sum(W[,j])
          Mb_loss[1]=Mb_loss[1]+sum(vali_W[,j]*(vali_sign_Con[,j]-current_fit[,j])^2)/vali_Con_ss[j]
        }
      }
      if (lossf=="logistic"){
        for (j in 1:q){
          current_fit[,j]=optim(0.5,line_schlog,current_fit=rep(0,n),step_fit=rep(1,n),Y=sign_Con[,j],weight=W[,j],method="BFGS")$par
          vali_Con_ss[j]=optim(0.5,line_schlog,current_fit=rep(0,n),step_fit=rep(1,n),Y=vali_sign_Con[,j],weight=vali_W[,j],method="BFGS")$value
          Mb_loss[1]=Mb_loss[1]+sum(vali_W[,j]*log(1+exp(-(2*vali_sign_Con[,j]-1)*(2*current_fit[,j]-1))))/vali_Con_ss[j]
        }
      }
      
      for (k in 2:(Mb_maxiter+1))
      {
        mod_sel=Mb_model[[1]][[k-1]]
        cov_sel=Mb_model[[2]][k-1]
        res_sel=Mb_model[[3]][k-1]
        stepsize_sel=Mb_model[[4]][k-1]
        
        tempdata=as.data.frame(list(tempx=vali_X[,cov_sel]))
        current_fit[,res_sel]=current_fit[,res_sel]+predict(mod_sel,newdata=tempdata)*mu*stepsize_sel
        
        if (lossf=="logistic"){
          for (j in 1:q){
            Mb_loss[k]=Mb_loss[k]+sum(vali_W[,j]*log(1+exp(-(2*vali_sign_Con[,j]-1)*(2*current_fit[,j]-1))))/vali_Con_ss[j]
          }
        }
        if (lossf=="square"){
          for (j in 1:q){
            Mb_loss[k]=Mb_loss[k]+sum(vali_W[,j]*(vali_sign_Con[,j]-current_fit[,j])^2)/vali_Con_ss[j]
          }
        }
      } 
      Mb_opt_iter=which.min(Mb_loss)
      
      
      ###########################
      #Uboost
      ###########################
      Ub_maxiter=400
      Ub_opt_iter=rep(0,q)
      Ub_model=rep(list(list()),q)
      Ub_loss_red_ctr=matrix(0,nrow=p,ncol=q)
      current_fit=matrix(0,nrow=n,ncol=q)
      for (j in 1:q){
        if (lossf=="square")
          current_fit[,j]=sum(sign_Con[,j]*W[,j])/sum(W[,j])
        if (lossf=="logistic")
          current_fit[,j]=optim(0.5,line_schlog,current_fit=rep(0,n),step_fit=rep(1,n),Y=sign_Con[,j],weight=W[,j],method="BFGS")$par
      }
      for (j in 1:q)
      {
        Ub_model[[j]]=Multiboost(Y=as.matrix(sign_Con[,j]),W=as.matrix(W[,j]),X=X,maxiter=Ub_maxiter,mu=mu,lossf=lossf,base=base)
        Ub_loss=rep(0,Ub_maxiter+1)
        Ub_loss_red_ctr[,j]=Ub_model[[j]][[5]]
        if (lossf=="square")
          Ub_loss[1]=sum(vali_W[,j]*(vali_sign_Con[,j]-current_fit[,j])^2)
        if (lossf=="logistic")
          Ub_loss[1]=sum(vali_W[,j]*log(1+exp(-(2*vali_sign_Con[,j]-1)*(2*current_fit[,j]-1))))
        
        for (k in 2:(Ub_maxiter+1))
        {
          mod_sel=Ub_model[[j]][[1]][[k-1]]
          cov_sel=Ub_model[[j]][[2]][k-1]
          stepsize_sel=Ub_model[[j]][[4]][k-1]
          
          tempdata=as.data.frame(list(tempx=vali_X[,cov_sel]))
          current_fit[,j]=current_fit[,j]+predict(mod_sel,newdata=tempdata)*mu*stepsize_sel
          
          if (lossf=="logistic")
            Ub_loss[k]=sum(vali_W[,j]*log(1+exp(-(2*vali_sign_Con[,j]-1)*(2*current_fit[,j]-1))))
          if (lossf=="square"){
            Ub_loss[k]=sum(vali_W[,j]*(vali_sign_Con[,j]-current_fit[,j])^2)
          }
        } 
        Ub_opt_iter[j]=which.min(Ub_loss)
      }
      
      #Decision of Mboost
      current_fit=matrix(0,nrow=ntest,ncol=q)
      if (lossf=="square"){
        for (j in 1:q)
          current_fit[,j]=sum(sign_Con[,j]*W[,j])/sum(W[,j])
      }
      if(lossf=="logistic"){
        for (j in 1:q)
          current_fit[,j]=optim(0.5,line_schlog,current_fit=rep(0,n),step_fit=rep(1,n),Y=sign_Con[,j],weight=W[,j],method="BFGS")$par
      } 
      if (Mb_opt_iter>1){
        for (k in 2:Mb_opt_iter){
          mod_sel=Mb_model[[1]][[k-1]]
          cov_sel=Mb_model[[2]][k-1]
          res_sel=Mb_model[[3]][k-1]
          stepsize_sel=Mb_model[[4]][k-1]
          
          tempdata=as.data.frame(list(tempx=Xtest[,cov_sel+1]))
          current_fit[,res_sel]=current_fit[,res_sel]+predict(mod_sel,newdata=tempdata)*mu*stepsize_sel
        }
      }
      Mboost_dec=sign(current_fit-0.5)
      Mboost_value=current_fit-0.5
      
      #Decision of Uboost
      current_fit=matrix(0,nrow=ntest,ncol=q)
      for (j in 1:q){
        if (lossf=="square")
          current_fit[,j]=sum(sign_Con[,j]*W[,j])/sum(W[,j])
        if (lossf=="logistic")
          current_fit[,j]=optim(0.5,line_schlog,current_fit=rep(0,n),step_fit=rep(1,n),Y=sign_Con[,j],weight=W[,j],method="BFGS")$par
        if (Ub_opt_iter[j]>1){
          for (k in 2:Ub_opt_iter[j]){
            mod_sel=Ub_model[[j]][[1]][[k-1]]
            cov_sel=Ub_model[[j]][[2]][k-1]
            stepsize_sel=Ub_model[[j]][[4]][k-1]
            
            tempdata=as.data.frame(list(tempx=Xtest[,cov_sel+1]))
            current_fit[,j]=current_fit[,j]+predict(mod_sel,newdata=tempdata)*mu*stepsize_sel
          }
        }
      }
      Uboost_dec=sign(current_fit-0.5)
      Uboost_value=current_fit-0.5
      
      Mboost_rate=(colSums(opt_dec==Mboost_dec)+0.5*colSums(Mboost_dec==0))/ntest
      Uboost_rate=(colSums(opt_dec==Uboost_dec)+0.5*colSums(Uboost_dec==0))/ntest
      Mboost_cor=rep(0,q);Uboost_cor=rep(0,q)
      for (j in 1:q){
        Mboost_cor[j]=cor(Mboost_value[,j],opt_value[,j],method="spearman")
        Uboost_cor[j]=cor(Uboost_value[,j],opt_value[,j],method="spearman")
      }
      
      multiboost_rate[lossf_choice,base_choice,]=Mboost_rate
      uniboost_rate[lossf_choice,base_choice,]=Uboost_rate
      multiboost_cor[lossf_choice,base_choice,]=Mboost_cor
      uniboost_cor[lossf_choice,base_choice,]=Uboost_cor
      
      assign(paste(lossf,base,"Mb_iter",sep="_"),Mb_opt_iter)
      assign(paste(lossf,base,"Ub_iter",sep="_"),Ub_opt_iter)
      assign(paste(lossf,base,"Mb_cov_sel",sep="_"),Mb_model[[2]])
      assign(paste(lossf,base,"Mb_res_sel",sep="_"),Mb_model[[3]])     
      Ub_cov_sel=matrix(NA,nrow=q,ncol=Ub_maxiter)
      for (j in 1:q)
        Ub_cov_sel[j,]=Ub_model[[j]][[2]]
      assign(paste(lossf,base,"Ub_cov_sel",sep="_"),Ub_cov_sel)
    }

  
  result=list(resp_lasso_beta,resp_lasso_rate,resp_lasso_cor,
              lasso_beta,lasso_rate,lasso_cor,
              SGL_beta,SGL_rate,SGL_cor,
              multiboost_rate,uniboost_rate,multiboost_cor,uniboost_cor,
              square_linear_Mb_iter,square_linear_Ub_iter,square_linear_Mb_cov_sel,square_linear_Ub_cov_sel,square_linear_Mb_res_sel,
              square_spline_Mb_iter,square_spline_Ub_iter,square_spline_Mb_cov_sel,square_spline_Ub_cov_sel,square_spline_Mb_res_sel,
              square_stump_Mb_iter,square_stump_Ub_iter,square_stump_Mb_cov_sel,square_stump_Ub_cov_sel,square_stump_Mb_res_sel,
              logistic_linear_Mb_iter,logistic_linear_Ub_iter,logistic_linear_Mb_cov_sel,logistic_linear_Ub_cov_sel,logistic_linear_Mb_res_sel,
              logistic_spline_Mb_iter,logistic_spline_Ub_iter,logistic_spline_Mb_cov_sel,logistic_spline_Ub_cov_sel,logistic_spline_Mb_res_sel,
              logistic_stump_Mb_iter,logistic_stump_Ub_iter,logistic_stump_Mb_cov_sel,logistic_stump_Ub_cov_sel,logistic_stump_Mb_res_sel)
 
   return(result)
}