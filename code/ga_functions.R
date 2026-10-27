######################functions#########################

#tukey winsorization approach for truncating extreme outliers
tukey<-function(x){
  x[x>(quantile(x,probs=.75,names=FALSE,na.rm=TRUE)+3*IQR(x,na.rm=TRUE))]<-quantile(x,probs=.75,names=FALSE,na.rm=TRUE)+3*IQR(x,na.rm=TRUE)
  return(x)
}
#tukey winsorization, pareto standardization and log-transformation
pareto_log<-function(x){
  x[x<=0]<-mean(min(x[x>0],na.rm=TRUE),0)
  x<-log(x)
  x[x>(quantile(x,probs=.75,names=FALSE,na.rm=TRUE)+3*IQR(x,na.rm=TRUE))]<-quantile(x,probs=.75,names=FALSE,na.rm=TRUE)+3*IQR(x,na.rm=TRUE)
  x<-(x-mean(x,na.rm=TRUE))/sqrt(sd(x,na.rm=TRUE))
  x<-jitter(x,amount=0.00001)
  return(x)
}
#tukey winsorization and pareto standardization without log transformation
pareto_nolog<-function(x){
  x[x<=0]<-mean(min(x[x>0],na.rm=TRUE),0)
  x[x>(quantile(x,probs=.75,names=FALSE,na.rm=TRUE)+3*IQR(x,na.rm=TRUE))]<-quantile(x,probs=.75,names=FALSE,na.rm=TRUE)+3*IQR(x,na.rm=TRUE)
  x<-(x-mean(x,na.rm=TRUE))/sqrt(sd(x,na.rm=TRUE))
  x<-jitter(x,amount=0.00001)
  return(x)
}
#tukey winsorization with log transformation no standardization
tukey_log<-function(x){
  x[x<=0]<-mean(min(x[x>0],na.rm=TRUE),0)
  x<-log(x)
  x[x>(quantile(x,probs=.75,names=FALSE,na.rm=TRUE)+3*IQR(x,na.rm=TRUE))]<-quantile(x,probs=.75,names=FALSE,na.rm=TRUE)+3*IQR(x,na.rm=TRUE)
  return(x)
}

#tukey winsorization without log transformation and without standardization
tukey_nolog<-function(x){
  x[x<=0]<-mean(min(x[x>0],na.rm=TRUE),0)
  x[x>(quantile(x,probs=.75,names=FALSE,na.rm=TRUE)+3*IQR(x,na.rm=TRUE))]<-quantile(x,probs=.75,names=FALSE,na.rm=TRUE)+3*IQR(x,na.rm=TRUE)
  return(x)
}



prep.data<-function (data) {
  test<-select(data, study_id, gestage, gest_wk,sga3,sga10,lga90,premie,preterm,lbw, bweight,bweight_raw,Sex, multiple, N_17OHP, TREC,af_ratio, C4DC, TYR, C5, C5DC, Ala, Arg, bio, C0, C2, C3, C4, C5, C6, C8, C10, C12, C14, C16, C18, C10_1, C12_1, C14OH, C14_1, C14_2, C16OH, C16_1OH, C18OH, C18_1, C18_2, C18_1OH, C3DC, C4OH, C5OH, C5_1, C6DC, C8_1, Cit, GALT, Gly, IRT, LEU, MET, Orn, PHE, TSH, Val)
  test_covars<-select(test,study_id,gestage,gest_wk,sga3,sga10,lga90,premie,preterm,lbw,Sex, multiple,bweight_raw)
  
  test_trans <- select(test, N_17OHP, TREC, af_ratio, C4DC, TYR, C5DC, Ala, bweight,Arg, bio, C0, C2, C3, C4, C5, C6, C8, C10, C12, C14, C16, C18, C10_1, C12_1, C14OH, C14_1, C14_2, C16OH, C16_1OH, C18OH, C18_1, C18_2, C18_1OH, C3DC, C4OH, C5OH, C5_1, C6DC, C8_1, Cit, GALT, Gly, IRT, LEU, MET, Orn, PHE, TSH, Val)
  test_log <- select(test_trans, N_17OHP, TREC, C4DC, TYR, C5DC, Ala, Arg, bio, C0, C2, C3, C4, C5, C6, C8, C10, C12, C14, C16, C18, C10_1, C12_1, C14OH, C14_1, C14_2, C16OH, C16_1OH, C18OH, C18_1, C18_2, C18_1OH, C3DC, C4OH, C5OH, C5_1, C6DC, C8_1, Cit, GALT, Gly, IRT, LEU, MET, Orn, PHE, TSH, Val)
  test_nolog<-select(test_trans, af_ratio, bweight)
  
  nolog_test2<-apply(test_nolog,2,pareto_nolog)
  log_test2<-apply(test_log,2,pareto_log)
  test_final<-data.frame(test_covars,nolog_test2,log_test2)
}

mae_cv<-function(data,model,alpha,title){
  pred<-predict(object=model,newdata=data,alpha=alpha)
  resid<-pred-data$gestage
  avg_ga<-(data$gestage+pred)/2
  plotdata<-data.frame(data$gestage,data$bweight,data$bweight_raw,pred,resid,avg_ga)
  names(plotdata)<-c("gestage","bweight","bweight_raw","pred","resid","avg_ga")
  
  mae<-round(mean(abs(resid),na.rm=TRUE),3)
  rmse<-round(sqrt(mean(resid^2)),3)
  one_wk<-round(mean(abs(resid)<=1)*100,1)
  n<-length(resid)
  
  mae.term<-round(mean(abs(resid[data$gestage>=37]),na.rm=TRUE),3)
  rmse.term<-round(sqrt(mean(resid[data$gestage>=37]^2)),3)
  one_wk.term<-round(mean(abs(resid[data$gestage>=37])<=1)*100,1)
  n.term<-length(resid[data$gestage>=37])
  
  mae.preterm<-round(mean(abs(resid[data$gestage<37]),na.rm=TRUE),3)
  rmse.preterm<-round(sqrt(mean(resid[data$gestage<37]^2)),3)
  one_wk.preterm<-round(mean(abs(resid[data$gestage<37])<=1)*100,1)
  n.preterm<-length(resid[data$gestage<37])
  
  mae.sga10<-round(mean(abs(resid[data$sga10==1]),na.rm=TRUE),3)
  rmse.sga10<-round(sqrt(mean(resid[data$sga10==1]^2)),3)
  one_wk.sga10<-round(mean(abs(resid[data$sga10==1])<=1)*100,1)
  n.sga10<-length(resid[data$sga10==1])
  
  mae.sga3<-round(mean(abs(resid[data$sga3==1]),na.rm=TRUE),3)
  rmse.sga3<-round(sqrt(mean(resid[data$sga3==1]^2)),3)
  one_wk.sga3<-round(mean(abs(resid[data$sga3==1])<=1)*100,1)
  n.sga3<-length(resid[data$sga3==1])
  
  mae.lbw<-round(mean(abs(resid[data$lbw=="<2500 g"]),na.rm=TRUE),3)
  rmse.lbw<-round(sqrt(mean(resid[data$lbw=="<2500 g"]^2)),3)
  one_wk.lbw<-round(mean(abs(resid[data$lbw=="<2500 g"])<=1)*100,1)
  n.lbw<-length(resid[data$lbw=="<2500 g"])
  
  
  mae.multiple<-round(mean(abs(resid[data$multiple=="multiple"]),na.rm=TRUE),3)
  rmse.multiple<-round(sqrt(mean(resid[data$multiple=="multiple"]^2)),3)
  one_wk.multiple<-round(mean(abs(resid[data$multiple=="multiple"])<=1)*100,1)
  n.multiple<-length(resid[data$multiple=="multiple"])
  
 
  ba.plot1<-ggplot(data=plotdata,mapping=aes(y=resid,x=avg_ga),na.rm=TRUE) +
    geom_point(size=1,shape=1,alpha=.2) +
    geom_abline(intercept=-1,slope=0,linetype=2,linewidth=1) +
    geom_abline(intercept=0,slope=0,linetype=1,linewidth=1) + 
    geom_abline(intercept=1,slope=0,linetype=2,linewidth=1) +
    geom_smooth(se=FALSE) +
    geom_jitter(alpha=.2) +
    labs(title=paste(title))
  ba.plot2<-ggplot(data=plotdata,mapping=aes(y=resid,x=gestage),na.rm=TRUE) +
    geom_point(size=1,shape=1,alpha=.2) +
    geom_abline(intercept=-1,slope=0,linetype=2,linewidth=1) +
    geom_abline(intercept=0,slope=0,linetype=1,linewidth=1) + 
    geom_abline(intercept=1,slope=0,linetype=2,linewidth=1) +
    geom_smooth(se=FALSE) +
    geom_jitter(alpha=.2)+
    labs(title=paste(title))
  ba.plot3<-ggplot(data=plotdata,mapping=aes(y=resid,x=pred),na.rm=TRUE) +
    geom_point(size=1,shape=1,alpha=.2) +
    geom_abline(intercept=-1,slope=0,linetype=2,linewidth=1) +
    geom_abline(intercept=0,slope=0,linetype=1,linewidth=1) + 
    geom_abline(intercept=1,slope=0,linetype=2,linewidth=1) +
    geom_smooth(se=FALSE) +
    geom_jitter(alpha=.2)+
    labs(title=paste(title))
  ba.plot4<-ggplot(data=plotdata,mapping=aes(y=resid,x=bweight_raw),na.rm=TRUE) +
    geom_point(size=1,shape=1,alpha=.2) +
    geom_abline(intercept=-1,slope=0,linetype=2,linewidth=1) +
    geom_abline(intercept=0,slope=0,linetype=1,linewidth=1) + 
    geom_abline(intercept=1,slope=0,linetype=2,linewidth=1) +
    geom_smooth(se=FALSE) +
    geom_jitter(alpha=.2)+
    labs(title=paste(title))
  ba.plot.combined<-grid.arrange(ba.plot1,ba.plot2,ba.plot3,ba.plot4,nrow=2)
  return(list(paste("MAE (Overall) N=",n,"     MAE=",mae,"     RMSE=",rmse,"     +/- 1wk=",one_wk,"%"),
              paste("MAE (Term) N=",n.term,"     MAE=",mae.term,"     RMSE=",rmse.term,"     +/- 1wk=",one_wk.term,"%"),
              paste("MAE (Preterm) N=",n.preterm,"     MAE=",mae.preterm,"     RMSE=",rmse.preterm,"     +/- 1wk=",one_wk.preterm,"%"),
              paste("MAE (sga10) N=",n.sga10,"     MAE=",mae.sga10,"     RMSE=",rmse.sga10,"     +/- 1wk=",one_wk.sga10,"%"),
              paste("MAE (sga3) N=",n.sga3,"     MAE=",mae.sga3,"     RMSE=",rmse.sga3,"     +/- 1wk=",one_wk.sga3,"%"),
              paste("MAE (lbw) N=",n.lbw,"     MAE=",mae.lbw,"     RMSE=",rmse.lbw,"     +/- 1wk=",one_wk.lbw,"%"),
              paste("MAE (multiple) N=",n.multiple,"     MAE=",mae.multiple,"     RMSE=",rmse.multiple,"     +/- 1wk=",one_wk.multiple,"%"),
              ba.plot.combined))
  
}

mae_cv.boot<-function(data,model,alpha) {
  pred<-predict(newdata=data,object=model,alpha=alpha)
  resid<-pred-data$gestage
  mae<-round(mean(abs(resid),na.rm=TRUE),3)
  rmse<-round(sqrt(mean(resid^2,na.rm=TRUE)),3)
  one_wk<-round(mean(abs(resid)<=1,na.rm=TRUE)*100,1)
  
  mae.term<-round(mean(abs(resid[data$gestage>=37]),na.rm=TRUE),3)
  rmse.term<-round(sqrt(mean(resid[data$gestage>=37]^2,na.rm=TRUE)),3)
  one_wk.term<-round(mean(abs(resid[data$gestage>=37])<=1,na.rm=TRUE)*100,1)
  
  
  mae.preterm<-round(mean(abs(resid[data$gestage<37]),na.rm=TRUE),3)
  rmse.preterm<-round(sqrt(mean(resid[data$gestage<37]^2,na.rm=TRUE)),3)
  one_wk.preterm<-round(mean(abs(resid[data$gestage<37])<=1,na.rm=TRUE)*100,1)
 
  
  mae.sga10<-round(mean(abs(resid[data$sga10==1]),na.rm=TRUE),3)
  rmse.sga10<-round(sqrt(mean(resid[data$sga10==1]^2,na.rm=TRUE)),3)
  one_wk.sga10<-round(mean(abs(resid[data$sga10==1])<=1,na.rm=TRUE)*100,1)
  
  
  mae.sga3<-round(mean(abs(resid[data$sga3==1]),na.rm=TRUE),3)
  rmse.sga3<-round(sqrt(mean(resid[data$sga3==1]^2,na.rm=TRUE)),3)
  one_wk.sga3<-round(mean(abs(resid[data$sga3==1])<=1,na.rm=TRUE)*100,1)
 
  
  mae.lbw<-round(mean(abs(resid[data$lbw=="<2500 g"]),na.rm=TRUE),3)
  rmse.lbw<-round(sqrt(mean(resid[data$lbw=="<2500 g"]^2,na.rm=TRUE)),3)
  one_wk.lbw<-round(mean(abs(resid[data$lbw=="<2500 g"])<=1,na.rm=TRUE)*100,1)

  
  
  mae.multiple<-round(mean(abs(resid[data$multiple=="multiple"]),na.rm=TRUE),3)
  rmse.multiple<-round(sqrt(mean(resid[data$multiple=="multiple"]^2,na.rm=TRUE)),3)
  one_wk.multiple<-round(mean(abs(resid[data$multiple=="multiple"])<=1,na.rm=TRUE)*100,1)
 
  
valid_stats<-data.frame(mae,rmse,one_wk,mae.term,rmse.term,one_wk.term,mae.preterm,rmse.preterm,
one_wk.preterm,mae.sga10,rmse.sga10,one_wk.sga10,mae.sga3,rmse.sga3,one_wk.sga3,mae.lbw,rmse.lbw,
one_wk.lbw,mae.multiple,rmse.multiple,one_wk.multiple)

}


##bootstrap function for generating bootstrap samples
run_mae_boot<-function(data,reps,seed,model,alpha,title,digits=2) {
set.seed(seed) # Setting the seed for replication purposes
model<-model
alpha<-alpha
title<-title
digits<-digits
bootstrap.results <- data.frame(matrix(nrow=reps,ncol=22)) # Creating an empty dataframe
names(bootstrap.results)<-c("rep","mae","rmse","one_wk","mae.term","rmse.term","one_wk.term","mae.preterm","rmse.preterm",
                            "one_wk.preterm","mae.sga10","rmse.sga10","one_wk.sga10","mae.sga3","rmse.sga3","one_wk.sga3","mae.lbw","rmse.lbw",
                            "one_wk.lbw","mae.multiple","rmse.multiple","one_wk.multiple")
for (i in 1:reps)
{
  obs <- sample(1:nrow(data),size=nrow(data),replace=TRUE)
  bootstrap.results[i,1]<-i
  bootstrap.results[i,2:22] <- mae_cv.boot(data[obs,],alpha=alpha,model=model) # stats of the bootstrap sample
}
bootstrap.mean<-apply(bootstrap.results[,-1],2,mean,na.rm=TRUE)
bootstrap.median<-apply(bootstrap.results[,-1],2,quantile,probs=0.5,na.rm=TRUE)
bootstrap.lcl<-apply(bootstrap.results[,-1],2,quantile,probs=0.025,na.rm=TRUE)
bootstrap.ucl<-apply(bootstrap.results[,-1],2,quantile,probs=0.975,na.rm=TRUE)
n.overall<-length(data$gestage)
n.term<-length(data$gestage[data$gestage>=37])
n.preterm<-length(data$gestage[data$gestage<37])
n.sga10<-length(data$sga10[data$sga10==1])
n.sga3<-length(data$sga3[data$sga3==1])
n.lbw<-length(data$lbw[data$lbw=="<2500 g"])
n.multiple<-length(data$multiple[data$multiple=="multiple"])

n.overall<-paste0("n= ",n.overall)
n.term<-paste0("n.term= ",n.term)
n.preterm<-paste0("n.preterm= ",n.preterm)
n.sga10<-paste0("n.sga10= ",n.sga10)
n.sga3<-paste0("n.sga3= ",n.sga3)
n.lbw<-paste0("n.lbw= ",n.lbw)
n.multiple<-paste0("n.multiple= ",n.multiple)

mae<-paste0("mae: ",round(bootstrap.median["mae"],digits)," (",round(bootstrap.lcl["mae"],digits),",",round(bootstrap.ucl["mae"],digits),")")
rmse<-paste0("rmse: ",round(bootstrap.median["rmse"],digits)," (",round(bootstrap.lcl["rmse"],digits),",",round(bootstrap.ucl["rmse"],digits),")")
one_wk<-paste0("one_wk: ",round(bootstrap.median["one_wk"],digits)," (",round(bootstrap.lcl["one_wk"],digits),",",round(bootstrap.ucl["one_wk"],digits),")")
mae.term<-paste0("mae.term: ",round(bootstrap.median["mae.term"],digits)," (",round(bootstrap.lcl["mae.term"],digits),",",round(bootstrap.ucl["mae.term"],digits),")")
rmse.term<-paste0("rmse.term: ",round(bootstrap.median["rmse.term"],digits)," (",round(bootstrap.lcl["rmse.term"],digits),",",round(bootstrap.ucl["rmse.term"],digits),")")
one_wk.term<-paste0("one_wk.term: ",round(bootstrap.median["one_wk.term"],digits)," (",round(bootstrap.lcl["one_wk.term"],digits),",",round(bootstrap.ucl["one_wk.term"],digits),")")
mae.preterm<-paste0("mae.preterm: ",round(bootstrap.median["mae.preterm"],digits)," (",round(bootstrap.lcl["mae.preterm"],digits),",",round(bootstrap.ucl["mae.preterm"],digits),")")
rmse.preterm<-paste0("rmse.preterm: ",round(bootstrap.median["rmse.preterm"],digits)," (",round(bootstrap.lcl["rmse.preterm"],digits),",",round(bootstrap.ucl["rmse.preterm"],digits),")")
one_wk.preterm<-paste0("one_wk.preterm: ",round(bootstrap.median["one_wk.preterm"],digits)," (",round(bootstrap.lcl["one_wk.preterm"],digits),",",round(bootstrap.ucl["one_wk.preterm"],digits),")")
mae.sga10<-paste0("mae.sga10: ",round(bootstrap.median["mae.sga10"],digits)," (",round(bootstrap.lcl["mae.sga10"],digits),",",round(bootstrap.ucl["mae.sga10"],digits),")")
rmse.sga10<-paste0("rmse.sga10: ",round(bootstrap.median["rmse.sga10"],digits)," (",round(bootstrap.lcl["rmse.sga10"],digits),",",round(bootstrap.ucl["rmse.sga10"],digits),")")
one_wk.sga10<-paste0("one_wk.sga10: ",round(bootstrap.median["one_wk.sga10"],digits)," (",round(bootstrap.lcl["one_wk.sga10"],digits),",",round(bootstrap.ucl["one_wk.sga10"],digits),")")
mae.sga3<-paste0("mae.sga3: ",round(bootstrap.median["mae.sga3"],digits)," (",round(bootstrap.lcl["mae.sga3"],digits),",",round(bootstrap.ucl["mae.sga3"],digits),")")
rmse.sga3<-paste0("rmse.sga3: ",round(bootstrap.median["rmse.sga3"],digits)," (",round(bootstrap.lcl["rmse.sga3"],digits),",",round(bootstrap.ucl["rmse.sga3"],digits),")")
one_wk.sga3<-paste0("one_wk.sga3: ",round(bootstrap.median["one_wk.sga3"],digits)," (",round(bootstrap.lcl["one_wk.sga3"],digits),",",round(bootstrap.ucl["one_wk.sga3"],digits),")")
mae.lbw<-paste0("mae.lbw: ",round(bootstrap.median["mae.lbw"],digits)," (",round(bootstrap.lcl["mae.lbw"],digits),",",round(bootstrap.ucl["mae.lbw"],digits),")")
rmse.lbw<-paste0("rmse.lbw: ",round(bootstrap.median["rmse.lbw"],digits)," (",round(bootstrap.lcl["rmse.lbw"],digits),",",round(bootstrap.ucl["rmse.lbw"],digits),")")
one_wk.lbw<-paste0("one_wk.lbw: ",round(bootstrap.median["one_wk.lbw"],digits)," (",round(bootstrap.lcl["one_wk.lbw"],digits),",",round(bootstrap.ucl["one_wk.lbw"],digits),")")
mae.multiple<-paste0("mae.multiple: ",round(bootstrap.median["mae.multiple"],digits)," (",round(bootstrap.lcl["mae.multiple"],digits),",",round(bootstrap.ucl["mae.multiple"],digits),")")
rmse.multiple<-paste0("rmse.multiple: ",round(bootstrap.median["rmse.multiple"],digits)," (",round(bootstrap.lcl["rmse.multiple"],digits),",",round(bootstrap.ucl["rmse.multiple"],digits),")")
one_wk.multiple<-paste0("one_wk.multiple: ",round(bootstrap.median["one_wk.multiple"],digits)," (",round(bootstrap.lcl["one_wk.multiple"],digits),",",round(bootstrap.ucl["one_wk.multiple"],digits),")")

cat("___________________________",title,"___________________________",n.overall,mae,rmse,one_wk,n.term,mae.term,rmse.term,one_wk.term,n.preterm,mae.preterm,rmse.preterm,
one_wk.preterm,n.sga10,mae.sga10,rmse.sga10,one_wk.sga10,n.sga3,mae.sga3,rmse.sga3,one_wk.sga3,n.lbw,mae.lbw,rmse.lbw,
one_wk.lbw,n.multiple,mae.multiple,rmse.multiple,one_wk.multiple,"___________________________",sep="\n")
}

##bootstrap function for existing bootstrap samples
run_mae_boot2<-function(data,reps,model,alpha,title,digits=2) {
  model<-model
  alpha<-alpha
  title<-title
  digits<-digits
  bootstrap.results <- data.frame(matrix(nrow=reps,ncol=22)) # Creating an empty dataframe
  names(bootstrap.results)<-c("rep","mae","rmse","one_wk","mae.term","rmse.term","one_wk.term","mae.preterm","rmse.preterm",
                              "one_wk.preterm","mae.sga10","rmse.sga10","one_wk.sga10","mae.sga3","rmse.sga3","one_wk.sga3","mae.lbw","rmse.lbw",
                              "one_wk.lbw","mae.multiple","rmse.multiple","one_wk.multiple")
  for (i in 1:reps)
  {
    obs <- filter(data,imp==i)
    bootstrap.results[i,1]<-i
    bootstrap.results[i,2:22] <- mae_cv.boot(obs,alpha=alpha,model=model) # stats of the bootstrap sample
  }
  bootstrap.mean<-apply(bootstrap.results[,-1],2,mean,na.rm=TRUE)
  bootstrap.median<-apply(bootstrap.results[,-1],2,quantile,probs=0.5,na.rm=TRUE)
  bootstrap.lcl<-apply(bootstrap.results[,-1],2,quantile,probs=0.025,na.rm=TRUE)
  bootstrap.ucl<-apply(bootstrap.results[,-1],2,quantile,probs=0.975,na.rm=TRUE)
  n.overall<-round(length(data$gestage)/reps)
  n.term<-round(length(data$gestage[data$gestage>=37])/reps)
  n.preterm<-round(length(data$gestage[data$gestage<37])/reps)
  n.sga10<-round(length(data$sga10[data$sga10==1])/reps)
  n.sga3<-round(length(data$sga3[data$sga3==1])/reps)
  n.lbw<-round(length(data$lbw[data$lbw=="<2500 g"])/reps)
  n.multiple<-round(length(data$multiple[data$multiple=="multiple"])/reps)
  n.overall<-paste0("n= ",n.overall)
  n.term<-paste0("n.term= ",n.term)
  n.preterm<-paste0("n.preterm= ",n.preterm)
  n.sga10<-paste0("n.sga10= ",n.sga10)
  n.sga3<-paste0("n.sga3= ",n.sga3)
  n.lbw<-paste0("n.lbw= ",n.lbw)
  n.multiple<-paste0("n.multiple= ",n.multiple)
  
  mae<-paste0("mae: ",round(bootstrap.median["mae"],digits)," (",round(bootstrap.lcl["mae"],digits),",",round(bootstrap.ucl["mae"],digits),")")
  rmse<-paste0("rmse: ",round(bootstrap.median["rmse"],digits)," (",round(bootstrap.lcl["rmse"],digits),",",round(bootstrap.ucl["rmse"],digits),")")
  one_wk<-paste0("one_wk: ",round(bootstrap.median["one_wk"],digits)," (",round(bootstrap.lcl["one_wk"],digits),",",round(bootstrap.ucl["one_wk"],digits),")")
  mae.term<-paste0("mae.term: ",round(bootstrap.median["mae.term"],digits)," (",round(bootstrap.lcl["mae.term"],digits),",",round(bootstrap.ucl["mae.term"],digits),")")
  rmse.term<-paste0("rmse.term: ",round(bootstrap.median["rmse.term"],digits)," (",round(bootstrap.lcl["rmse.term"],digits),",",round(bootstrap.ucl["rmse.term"],digits),")")
  one_wk.term<-paste0("one_wk.term: ",round(bootstrap.median["one_wk.term"],digits)," (",round(bootstrap.lcl["one_wk.term"],digits),",",round(bootstrap.ucl["one_wk.term"],digits),")")
  mae.preterm<-paste0("mae.preterm: ",round(bootstrap.median["mae.preterm"],digits)," (",round(bootstrap.lcl["mae.preterm"],digits),",",round(bootstrap.ucl["mae.preterm"],digits),")")
  rmse.preterm<-paste0("rmse.preterm: ",round(bootstrap.median["rmse.preterm"],digits)," (",round(bootstrap.lcl["rmse.preterm"],digits),",",round(bootstrap.ucl["rmse.preterm"],digits),")")
  one_wk.preterm<-paste0("one_wk.preterm: ",round(bootstrap.median["one_wk.preterm"],digits)," (",round(bootstrap.lcl["one_wk.preterm"],digits),",",round(bootstrap.ucl["one_wk.preterm"],digits),")")
  mae.sga10<-paste0("mae.sga10: ",round(bootstrap.median["mae.sga10"],digits)," (",round(bootstrap.lcl["mae.sga10"],digits),",",round(bootstrap.ucl["mae.sga10"],digits),")")
  rmse.sga10<-paste0("rmse.sga10: ",round(bootstrap.median["rmse.sga10"],digits)," (",round(bootstrap.lcl["rmse.sga10"],digits),",",round(bootstrap.ucl["rmse.sga10"],digits),")")
  one_wk.sga10<-paste0("one_wk.sga10: ",round(bootstrap.median["one_wk.sga10"],digits)," (",round(bootstrap.lcl["one_wk.sga10"],digits),",",round(bootstrap.ucl["one_wk.sga10"],digits),")")
  mae.sga3<-paste0("mae.sga3: ",round(bootstrap.median["mae.sga3"],digits)," (",round(bootstrap.lcl["mae.sga3"],digits),",",round(bootstrap.ucl["mae.sga3"],digits),")")
  rmse.sga3<-paste0("rmse.sga3: ",round(bootstrap.median["rmse.sga3"],digits)," (",round(bootstrap.lcl["rmse.sga3"],digits),",",round(bootstrap.ucl["rmse.sga3"],digits),")")
  one_wk.sga3<-paste0("one_wk.sga3: ",round(bootstrap.median["one_wk.sga3"],digits)," (",round(bootstrap.lcl["one_wk.sga3"],digits),",",round(bootstrap.ucl["one_wk.sga3"],digits),")")
  mae.lbw<-paste0("mae.lbw: ",round(bootstrap.median["mae.lbw"],digits)," (",round(bootstrap.lcl["mae.lbw"],digits),",",round(bootstrap.ucl["mae.lbw"],digits),")")
  rmse.lbw<-paste0("rmse.lbw: ",round(bootstrap.median["rmse.lbw"],digits)," (",round(bootstrap.lcl["rmse.lbw"],digits),",",round(bootstrap.ucl["rmse.lbw"],digits),")")
  one_wk.lbw<-paste0("one_wk.lbw: ",round(bootstrap.median["one_wk.lbw"],digits)," (",round(bootstrap.lcl["one_wk.lbw"],digits),",",round(bootstrap.ucl["one_wk.lbw"],digits),")")
  mae.multiple<-paste0("mae.multiple: ",round(bootstrap.median["mae.multiple"],digits)," (",round(bootstrap.lcl["mae.multiple"],digits),",",round(bootstrap.ucl["mae.multiple"],digits),")")
  rmse.multiple<-paste0("rmse.multiple: ",round(bootstrap.median["rmse.multiple"],digits)," (",round(bootstrap.lcl["rmse.multiple"],digits),",",round(bootstrap.ucl["rmse.multiple"],digits),")")
  one_wk.multiple<-paste0("one_wk.multiple: ",round(bootstrap.median["one_wk.multiple"],digits)," (",round(bootstrap.lcl["one_wk.multiple"],digits),",",round(bootstrap.ucl["one_wk.multiple"],digits),")")
  
  cat("___________________________",title,"___________________________",n.overall,mae,rmse,one_wk,n.term,mae.term,rmse.term,one_wk.term,n.preterm,mae.preterm,rmse.preterm,
      one_wk.preterm,n.sga10,mae.sga10,rmse.sga10,one_wk.sga10,n.sga3,mae.sga3,rmse.sga3,one_wk.sga3,n.lbw,mae.lbw,rmse.lbw,
      one_wk.lbw,n.multiple,mae.multiple,rmse.multiple,one_wk.multiple,"___________________________",sep="\n")

  }

mae_plots<-function(data,model,alpha,title){
  pred<-predict(object=model,newdata=data,alpha=alpha)
  resid<-pred-data$gestage
  avg_ga<-(data$gestage+pred)/2
  plotdata<-data.frame(data$gestage,data$bweight,data$bweight_raw,pred,resid,avg_ga)
  names(plotdata)<-c("gestage","bweight","bweight_raw","pred","resid","avg_ga")

  ba.plot1<-ggplot(data=plotdata,mapping=aes(y=resid,x=avg_ga),na.rm=TRUE) +
    geom_point(size=1,shape=1,alpha=.2) +
    geom_abline(intercept=-1,slope=0,linetype=2,linewidth=1) +
    geom_abline(intercept=0,slope=0,linetype=1,linewidth=1) + 
    geom_abline(intercept=1,slope=0,linetype=2,linewidth=1) +
    geom_smooth(se=FALSE) +
    geom_jitter(alpha=.2) +
    labs(title=paste(title))
  ba.plot2<-ggplot(data=plotdata,mapping=aes(y=resid,x=gestage),na.rm=TRUE) +
    geom_point(size=1,shape=1,alpha=.2) +
    geom_abline(intercept=-1,slope=0,linetype=2,linewidth=1) +
    geom_abline(intercept=0,slope=0,linetype=1,linewidth=1) + 
    geom_abline(intercept=1,slope=0,linetype=2,linewidth=1) +
    geom_smooth(se=FALSE) +
    geom_jitter(alpha=.2)+
    labs(title=paste(title))
  ba.plot3<-ggplot(data=plotdata,mapping=aes(y=resid,x=pred),na.rm=TRUE) +
    geom_point(size=1,shape=1,alpha=.2) +
    geom_abline(intercept=-1,slope=0,linetype=2,linewidth=1) +
    geom_abline(intercept=0,slope=0,linetype=1,linewidth=1) + 
    geom_abline(intercept=1,slope=0,linetype=2,linewidth=1) +
    geom_smooth(se=FALSE) +
    geom_jitter(alpha=.2)+
    labs(title=paste(title))
  ba.plot4<-ggplot(data=plotdata,mapping=aes(y=resid,x=bweight_raw),na.rm=TRUE) +
    geom_point(size=1,shape=1,alpha=.2) +
    geom_abline(intercept=-1,slope=0,linetype=2,linewidth=1) +
    geom_abline(intercept=0,slope=0,linetype=1,linewidth=1) + 
    geom_abline(intercept=1,slope=0,linetype=2,linewidth=1) +
    geom_smooth(se=FALSE) +
    geom_jitter(alpha=.2)+
    labs(title=paste(title))
  ba.plot.combined<-gridExtra::grid.arrange(ba.plot1,ba.plot2,ba.plot3,ba.plot4,nrow=2)
#return(ba.plot.combined)
  
}

# Bootstrap Sampling plus single imputation at each step
bootsamples<-function(impnum,data,seed) {
  set.seed(seed)  
  for (imp in 1:impnum)  {
    print(paste("imputation # ",imp))
    boot<-data[sample(nrow(data),nrow(data),replace=TRUE),]
    boot.im = aregImpute(~ bweight + af_ratio + N_17OHP + TYR + C5 + C4DC + 
                           C5DC + TREC + Ala + Sex + af_ratio + 
                           bweight + N_17OHP + TREC + C4DC + TYR + C5 + C5DC + 
                           Ala + Arg + bio + C0 + C2 + C3 + C4 + C6 + C8 + 
                           C10 + C12 + C14 + C16 + C18 + C10_1 + C12_1 + 
                           C14OH + C14_1 + C14_2 + C16OH + C16_1OH + C18OH + 
                           C18_1 + C18_2 + C18_1OH + C3DC + C4OH + C5OH + 
                           C5_1 + C6DC + C8_1 + Cit + GALT + Gly + IRT + 
                           LEU + MET + Orn + PHE + TSH + Val,
                         data = boot,nk=3,n.impute=1)
    
    boot.imputed<-cbind.data.frame(impute.transcan(x = boot.im, 
                                                   imputation = 1, 
                                                   data = boot, 
                                                   list.out = TRUE, 
                                                   pr = FALSE))
    boot.combined <- data.frame(select(boot,study_id,gestage,gest_wk,sga3,sga10,lga90,
                                       premie,preterm,lbw,multiple,bweight_raw),boot.imputed,imp)
    boot.combined <-select(boot.combined,study_id,gestage,N_17OHP, 
                           TREC,af_ratio, C4DC, TYR, C5, C5DC,Ala, bweight,bweight_raw,
                           Arg, bio, C0, C2, C3, C4,C5, C6, C8, C10, C12, C14, C16, 
                           C18, C10_1, C12_1, C14OH, C14_1, C14_2, C16OH, C16_1OH, C18OH, 
                           C18_1, C18_2, C18_1OH, C3DC, C4OH, C5OH, C5_1, C6DC, C8_1, Cit, 
                           GALT, Gly, IRT, LEU, MET, Orn, PHE, TSH, Val,gest_wk, premie,
                           preterm,lbw,Sex, multiple,sga3,sga10,lga90,imp)
    boot.combined[,2:54]<-apply(boot.combined[,2:54],2,num)                   
    if (imp==1) { boot.collected<- boot.combined }
    if (imp!=1) {boot.collected<-rbind(boot.collected,boot.combined) }
    rm("boot","boot.imputed","boot.combined")
  }
  return(boot.collected)  
}



