line_schlog=function(rho,current_fit,step_fit,Y,weight){
  return(sum(weight*log(1+exp(-(2*Y-1)*(2*(current_fit+rho*step_fit)-1)))))
}

Multiboost<-function(Y,W,X,maxiter,mu,lossf=c("logistic","square"),base=c("spline","linear","stump")){
  n=dim(X)[1];p=dim(X)[2];q=dim(Y)[2]
  mod_list=rep(list(list()),maxiter)
  cov_list=numeric(maxiter)
  res_list=numeric(maxiter)
  loss_reduce=numeric(maxiter)
  stepsize_list=numeric(maxiter)
  loss_red_ctr=matrix(0,nrow=p,ncol=q)
  current_fit<-matrix(0,nrow=n,ncol=q)
  if (lossf=="logistic")
  {
    for (j in 1:q)
      current_fit[,j]=optim(0.5,line_schlog,current_fit=rep(0,n),step_fit=rep(1,n),Y=Y[,j],weight=W[,j],method="BFGS")$par
    current_loss=sum(W*log(1+exp(-(2*Y-1)*(2*current_fit-1))))
    M=1
    while(M<=maxiter){
      presid=(2*Y-1)/(exp((2*Y-1)*(2*current_fit-1))+1)
      loss=matrix(0,nrow=p,ncol=q)
      for (i in 1:p)
        for (j in 1:q){
          if (base=="spline"){
            tempdata=as.data.frame(list(tempx=X[,i],tempy=presid[,j]))
            mod=gamboost(tempy~tempx,data=tempdata,baselearner="bbs",weights=W[,j],control=boost_control(mstop=1,nu=1),dfbase=4)
            mod_fit=predict(mod)
          }
          if (base=="linear"){
            tempdata=as.data.frame(list(tempx=X[,i],tempy=presid[,j]))
            mod=gamboost(tempy~tempx,data=tempdata,baselearner="bols",weights=W[,j],control=boost_control(mstop=1,nu=1))
            mod_fit=predict(mod)
          }
          if (base=="stump"){
            tempdata=as.data.frame(list(tempx=X[,i],tempy=presid[,j]))
            mod=gamboost(tempy~tempx,data=tempdata,baselearner="btree",weights=W[,j],control=boost_control(mstop=1,nu=1))
            mod_fit=predict(mod)
          }
          loss[i,j]=sum(W[,j]*(presid[,j]-mod_fit)^2)/sum(W[,j]*presid[,j]^2)
        }
      
      cov_sel=which.min(rowSums(loss))
      res_sel=which.min(loss[cov_sel,])
      
      if (base=="spline"){
        tempdata=as.data.frame(list(tempx=X[,cov_sel],tempy=presid[,res_sel]))
        mod_sel=gamboost(tempy~tempx,data=tempdata,baselearner="bbs",weights=W[,res_sel],control=boost_control(mstop=1,nu=1),dfbase=4)
        step_fit=predict(mod_sel)
      }
      if (base=="linear"){
        tempdata=as.data.frame(list(tempx=X[,cov_sel],tempy=presid[,res_sel]))
        mod_sel=gamboost(tempy~tempx,data=tempdata,baselearner="bols",weights=W[,res_sel],control=boost_control(mstop=1,nu=1))
        step_fit=predict(mod_sel)
      }
      if (base=="stump"){
        tempdata=as.data.frame(list(tempx=X[,cov_sel],tempy=presid[,res_sel]))
        mod_sel=gamboost(tempy~tempx,data=tempdata,baselearner="btree",weights=W[,res_sel],control=boost_control(mstop=1,nu=1))
        step_fit=predict(mod_sel)
      }
      opt_stepsize=optim(0,line_schlog,current_fit=current_fit[,res_sel],step_fit=step_fit,Y=Y[,res_sel],weight=W[,res_sel],method="BFGS")$par
      current_fit[,res_sel]=current_fit[,res_sel]+step_fit*opt_stepsize*mu
      loss_reduce[M]=current_loss-sum(W*log(1+exp(-(2*Y-1)*(2*current_fit-1))))
      loss_red_ctr[cov_sel,res_sel]=loss_red_ctr[cov_sel,res_sel]+loss_reduce[M]
      current_loss=sum(W*log(1+exp(-(2*Y-1)*(2*current_fit-1))))
      mod_list[[M]]=mod_sel
      cov_list[M]=cov_sel
      res_list[M]=res_sel
      stepsize_list[M]=opt_stepsize
      M=M+1
    }
  }
  if (lossf=="square")
  {
    for (j in 1:q)
      current_fit[,j]=sum(Y[,j]*W[,j])/sum(W[,j])
    current_loss=sum(W*(Y-current_fit)^2)
    M=1
    while(M<=maxiter){
      presid=Y-current_fit
      loss=matrix(0,nrow=p,ncol=q)
      for (i in 1:p)
        for (j in 1:q){
          if (base=="spline"){
            tempdata=as.data.frame(list(tempx=X[,i],tempy=presid[,j]))
            mod=gamboost(tempy~tempx,data=tempdata,baselearner="bbs",weights=W[,j],control=boost_control(mstop=1,nu=1),dfbase=4)
            mod_fit=predict(mod)
          }
          if (base=="linear"){
            tempdata=as.data.frame(list(tempx=X[,i],tempy=presid[,j]))
            mod=gamboost(tempy~tempx,data=tempdata,baselearner="bols",weights=W[,j],control=boost_control(mstop=1,nu=1))
            mod_fit=predict(mod)
          }
          if (base=="stump"){
            tempdata=as.data.frame(list(tempx=X[,i],tempy=presid[,j]))
            mod=gamboost(tempy~tempx,data=tempdata,baselearner="btree",weights=W[,j],control=boost_control(mstop=1,nu=1))
            mod_fit=predict(mod)
          }
          loss[i,j]=sum(W[,j]*(presid[,j]-mod_fit)^2)/sum(W[,j]*presid[,j]^2)
        }
      
      cov_sel=which.min(rowSums(loss))
      res_sel=which.min(loss[cov_sel,])
      
      if (base=="spline"){
        tempdata=as.data.frame(list(tempx=X[,cov_sel],tempy=presid[,res_sel]))
        mod_sel=gamboost(tempy~tempx,data=tempdata,baselearner="bbs",weights=W[,res_sel],control=boost_control(mstop=1,nu=1),dfbase=4)
        step_fit=predict(mod_sel)
      }
      if (base=="linear"){
        tempdata=as.data.frame(list(tempx=X[,cov_sel],tempy=presid[,res_sel]))
        mod_sel=gamboost(tempy~tempx,data=tempdata,baselearner="bols",weights=W[,res_sel],control=boost_control(mstop=1,nu=1))
        step_fit=predict(mod_sel)
      }
      if (base=="stump"){
        tempdata=as.data.frame(list(tempx=X[,cov_sel],tempy=presid[,res_sel]))
        mod_sel=gamboost(tempy~tempx,data=tempdata,baselearner="btree",weights=W[,res_sel],control=boost_control(mstop=1,nu=1))
        step_fit=predict(mod_sel)
      }
      current_fit[,res_sel]=current_fit[,res_sel]+step_fit*mu
      loss_reduce[M]=current_loss-sum(W*(Y-current_fit)^2)
      loss_red_ctr[cov_sel,res_sel]=loss_red_ctr[cov_sel,res_sel]+loss_reduce[M]
      current_loss=sum(W*(Y-current_fit)^2)
      mod_list[[M]]=mod_sel
      cov_list[M]=cov_sel
      res_list[M]=res_sel
      stepsize_list[M]=1
      M=M+1
    }
  }
  return(list(mod_list,cov_list,res_list,stepsize_list,loss_red_ctr,loss_reduce))   
}
