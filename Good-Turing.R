
# Appendix S78 of Chao et al. (2016) "Deciphering the Enigma of Undetected Species, Phylogenetic, and Functional
# Diversity Based on Good-Turing Theory". Manuscript submitted to Ecology.

#### R scripts for estimating true species richness, phylogenetic diversity (PD) and FAD (functional attribute diversity) 
# in one assemblage, as well as for estimating shared species richness, shared PD and shared FAD between two 
# assemblages based on the Good-Turing frequency formula and its generalizations.

#----------------------------------------------
# Species richness estimator in one assemblage 
#----------------------------------------------
# Good_Turing_Richness is a function of obtaining estimators of species richness in one assemblage based on 
# abundance data.
# @para. abunDATA is a vector of species observed frequencies in a sample.
# @return species richness estimator and the associated 95% confidence interval (C.I.) 

Good_Turing_Richness=function(abunDATA){
  da=abunDATA;
  n=sum(da);obs=sum(da>0);
  
  f1=sum(da==1); f2=sum(da==2);
  us=round(ifelse(f1>0,f1^2/(2*f2),f1*(f1-1)/2 ));
  est=us+obs;
  ###95% C.I. of S
  fi=sapply(1:max(da),function(k) sum(da==k));
  V=fi[2]*(0.25*(fi[1]/fi[2])^4+(fi[1]/fi[2])^3+0.5*(fi[1]/fi[2])^2)
  R=exp(1.96*(log(1+V/(us)^2))^(1/2));
  est.CI=c(obs+us/R,obs+us*R);
  
  return(list( Richness=est, CI=est.CI))
}


#-----------------------------------------------------------
# Shared species richness estimator between two assemblages 
#-----------------------------------------------------------
# Good_Turing_shared.Richness is a function of obtaining estimators of shared species richness between two 
# assemblages based on abundance data.
# @para. abunDATA is an observed species-by-assemblage frequency matrix.
# @return shared species richness estimator and the associated 95% confidence interval. 

Good_Turing_shared.Richness=function(abunDATA){
  da1=abunDATA[,1];da2=abunDATA[,2];
  ##shared richness estimaiton
  I=which(da1*da2>0);sda1=da1[I];sda2=da2[I];obs12=length(I)
  f11=sum(sda1==1 & sda2==1);f22=sum(sda1==2 & sda2==2);f22=max(1,f22);
  f1_=sum(sda1==1); f2_=sum(sda1==2);f2_=max(1,f2_);
  f_1=sum(sda2==1); f_2=sum(sda2==2);f_2=max(1,f_2);
  us12=f_1^2/(2*f_2)+f1_^2/(2*f2_)+f11^2/(4*f22);
  s12=obs12+us12;
  ###95%CI of shared richness
  fii=c(obs12,f11,f22,f1_,f2_,f_1,f_2);
  df=c(1,f11/(2*f22),-(f11/f22)^2/4,f_1/f_2,-(f_1/f_2)^2/2,f1_/f2_,-(f_1/f_2)^2/2);
  COV=matrix(0,7,7);COV[7,7]=fii[7]*(1-fii[7]/s12);
  for(i in 1:6){
    COV[i,i]=fii[i]*(1-fii[i]/s12);
    for(j in (i+1):7){
      COV[i,j]=-fii[i]*fii[j]/s12;COV[i,j]=COV[j,i];}
  }
  V=t(df)%*%COV%*%df;
  R=exp(1.96*(log(1+V/(s12)^2))^(1/2));
  s12.CI=c(obs12+us12/R,obs12+us12*R) 
  return(list(shared.Richness=s12,CI=s12.CI))
}

#--------------------------------------------
# Estimator of Faith's PD in one assemblage 
#--------------------------------------------
# Good_Turing_PD is a function of obtaining estimators of PD based on abundance data.
# @para. tree is in the Newick format of a phylogenetic tree.
# @para. abunDATA is a vector of species observed frequencies in a sample.
# **NOTE: species names in abunDATA should be exactly the same as those in the Newick tree file.
# @return PD estimator and the associated 95% confidence interval. 

library(ade4);
Good_Turing_PD=function(tree,abunDATA){
  phyloData <- newick2phylog(tree);  
  abun <- abunDATA[names(phyloData$leaves)];
  nodenames=c(names(phyloData$leaves),names(phyloData$nodes));
  
  M=matrix(0,nrow=length(phyloData$leaves),ncol=length(nodenames),dimnames=list(names(phyloData$leaves),nodenames))
  for(i in 1:length(phyloData$leaves)){
    M[i,][unlist(phyloData$paths[i])]=rep(1,length(unlist(phyloData$paths[i])))}
  pda=t(abun)%*%M;
  BL=c(phyloData$leaves,phyloData$nodes)
  ######
  obspd=sum(BL[pda>0]);
  f1=sum(pda==1);f2=sum(pda==2);
  g1=sum(BL[pda==1]); g2=sum(BL[pda==2]);
  upd=ifelse( (2*g2*f1 <= g1*f2), g1*(f1-1)/(2*(f2+1)),(g1^2)/(2*g2))
  pd=obspd+upd;
  V=g2*((g1/g2)^4/4+(g1/g2)^3+(g1/g2)^2/2)
  R=exp(1.96*(log(1+V/(upd)^2))^(1/2));
  pd.CI=c(obspd+upd/R,obspd+upd*R) 
  return(list( PD=pd, CI=pd.CI))
}

#------------------------------------------------
# Shared PD estimator between two assemblages 
#------------------------------------------------
# Good_Turing_shared.PD is a function of obtaining estimators of shared PD between two assemblages
# based on abundance data.
# @para. tree is in the Newick format of a phylogenetic tree.
# @para. abunDATA is a species-by-assemblage frequency matrix.
# **NOTE: Species names in abunDATA should be exactly the same as those in the Newick tree file.
# @return shared PD estimator and the associated 95% confidence interval. 

library(ade4);
Good_Turing_shared.PD=function(tree,abunDATA){
  phyloData <- newick2phylog(tree)  
  abun <- as.matrix(abunDATA[names(phyloData$leaves), ])
  nodenames=c(names(phyloData$leaves),names(phyloData$nodes));
  
  M=matrix(0,nrow=length(phyloData$leaves),ncol=length(nodenames),dimnames=list(names(phyloData$leaves),nodenames))
  for(i in 1:length(phyloData$leaves)){
    M[i,][unlist(phyloData$paths[i])]=rep(1,length(unlist(phyloData$paths[i])))
  }
  phyloAbun=matrix(0,ncol=ncol(abunDATA),nrow=length(nodenames),dimnames=list(nodenames,colnames(abunDATA)))
  for(i in 1:ncol(abunDATA)){
    phyloAbun[,i]=abun[,i]%*%M;}
  BL=c(phyloData$leaves,phyloData$nodes)
  ######
  pda1=phyloAbun[,1];obspd1=sum(BL[pda1>0]);
  pda2=phyloAbun[,2];obspd2=sum(BL[pda2>0]);
  ###### CI of pd12
  I=which(pda1*pda2>0);sda1=pda1[I];sda2=pda2[I];sBL=BL[I];obspd12=sum(sBL);
  g11=sum(sBL[sda1==1 & sda2==1]);g22=sum(sBL[sda1==2 & sda2==2]);
  g1_=sum(sBL[sda1==1]);g2_=sum(sBL[sda1==2]);
  g_1=sum(sBL[sda2==1]);g_2=sum(sBL[sda2==2]);
  upd12=g_1^2/(2*g_2)+g1_^2/(2*g2_)+g11^2/(4*g22);spd12=obspd12+upd12;
  gii=c(obspd12,g11,g22,g1_,g2_,g_1,g_2);
  dg=c(1,g11/(2*g22),-(g11/g22)^2/4,g_1/g_2,-(g_1/g_2)^2/2,g1_/g2_,-(g_1/g_2)^2/2);
  COV=matrix(0,7,7);COV[7,7]=gii[7]*(1-gii[7]/spd12);
  for(i in 1:6){
    COV[i,i]=gii[i]*(1-gii[i]/spd12);
    for(j in (i+1):7){
      COV[i,j]=-gii[i]*gii[j]/spd12;COV[i,j]=COV[j,i];}
  }
  V=t(dg)%*%COV%*%dg;
  R=exp(1.96*(log(1+V/(spd12)^2))^(1/2));
  spd12.CI=c(obspd12+upd12/R,obspd12+upd12*R) 
  
  return(list(shared.PD=spd12,CI=spd12.CI))
}

#-------------------------------
# FAD estimator in an assemblage 
#-------------------------------
# Good_Turing_FD is a function of obtaining estimators of FAD in an assemblage based on abundance data.
# @para. Dis is a functional distance matrix of observed species.
# @para. abunDATA is a vector of species observed frequencies in a sample.

# **NOTE: species in functional distance matrix should be in the same ordering as those in species frequency matrix.
# @return FAD estimator and the associated 95% confidence interval. 


Good_Turing_FD=function(Dis,abunDATA){
  da=c(abunDATA);d=mean(Dis);
  I=which(da>0);D=Dis[I,I];ab=da[I];obsfd=sum(D);
  F1_=sum(D[ab==1,]);F2_=sum(D[ab==2,]);F2_=max(d,F2_);
  F_1=sum(D[,ab==1]);F_2=sum(D[,ab==2]);F_2=max(d,F_2);
  F11=sum(D[ab==1,ab==1]);F22=sum(D[ab==2,ab==2]);F22=max(d,F22)
  
  ufd=F_1^2/(2*F_2)+F1_^2/(2*F2_)+F11^2/(4*F22);
  fd=obsfd+ufd;
  
  Fii=c(obsfd,F11,F22,F1_,F2_,F_1,F_2);
  dF=c(1,F11/(2*F22),-(F11/F22)^2/4,F_1/F_2,-(F_1/F_2)^2/2,F1_/F2_,-(F_1/F_2)^2/2);
  COV=matrix(0,7,7);COV[7,7]=Fii[7]*(1-Fii[7]/fd);
  for(i in 1:6){
    COV[i,i]=Fii[i]*(1-Fii[i]/fd);
    for(j in (i+1):7){
      COV[i,j]=-Fii[i]*Fii[j]/fd;COV[i,j]=COV[j,i];}
  }
  V=t(dF)%*%COV%*%dF;
  R=exp(1.96*(log(1+V/(fd)^2))^(1/2));
  fd.CI=c(obsfd+ufd/R,obsfd+ufd*R) 
  
  return(list( FAD=fd, CI=fd.CI))
}

#-------------------------------
# Estimated bootstrap assemblage
#-------------------------------
# @para. Dis is a functional distance matrix of observed species.
# @para. abunDATA is an observed species-by-assemblage frequency matrix.
# **NOTE: species in functional distance matrix should be in the same ordering as those in species frequency matrix.
EstBootCommTwoFun = function(Dis,abunDATA){
  D=Dis;data=abunDATA;
  x1 = data[, 1]; x2 = data[, 2]
  n1 = sum(x1); n2 = sum(x2)
  index = which(x1>0 | x2>0)
  x1 = x1[index]; x2 = x2[index]
  D = D[index, index]
  S12.obs = sum(x1>0 & x2>0)
  S1.obs = sum(x1>0); S2.obs = sum(x2>0)
  #Basic information of relative abundance-------------------
  f.0 = sum(x1>0 & x2==0)
  f.1 = sum(x1>0 & x2==1)
  f.2 = sum(x1>0 & x2==2)
  f0. = sum(x1==0 & x2>0)
  f1. = sum(x1==1 & x2>0)
  f2. = sum(x1==2 & x2>0)
  f11 = sum(x1==1 & x2==1)
  f22 = sum(x1==2 & x2==2)
  f.0_hat = ceiling(ifelse(f.2>0, (n2-1)/n2*f.1^2/2/f.2, (n2-1)/n2*f.1*(f.1-1)/2))
  f0._hat = ceiling(ifelse(f2.>0, (n1-1)/n1*f1.^2/2/f2., (n1-1)/n1*f1.*(f1.-1)/2))
  f00_hat = ceiling(ifelse(f22>0, (n1-1)/n1*(n2-1)/n2*f11^2/4/f22, (n1-1)/n1*(n2-1)/n2*f11*(f11-1)/4))
  f1_1 = sum(x1 == 1)
  f2_1 = sum(x1 == 2)
  f0_1 = ceiling(ifelse(f2_1>0, (n1-1)/n1*f1_1^2/2/f2_1, (n1-1)/n1*f1_1*(f1_1-1)/2))
  f1_2 = sum(x2 == 1)
  f2_2 = sum(x2 == 2)
  f0_2 = ceiling(ifelse(f2_2>0, (n2-1)/n2*f1_2^2/2/f2_2, (n2-1)/n2*f1_2*(f1_2-1)/2))
  f0_1_p = max(f.0, f.0_hat)-f.0+f0._hat+f00_hat
  f0_2_p = max(f0., f0._hat)-f0.+f.0_hat+f00_hat
  f0_1_star = max(f0_1, f0_1_p)
  f0_2_star = max(f0_2, f0_2_p)
  n00 = f00_hat + (max(f.0, f.0_hat) - f.0) + (max(f0., f0._hat) - f0.)
  n.0 = min(S1.obs-S12.obs, f.0_hat)
  n0. = min(S2.obs-S12.obs, f0._hat)
  
  if(f2_1>0){
    A_1 = 2*f2_1/((n1-1)*f1_1+2*f2_1)
  }else if(f1_1>1){
    A_1 = 2/((n1-1)*(f1_1-1)+2)
  }else{
    A_1 = 0
  }
  if(f2_2>0){
    A_2 = 2*f2_2/((n2-1)*f1_2+2*f2_2)
  }else if(f1_2>1){
    A_2 = 2/((n2-1)*(f1_2-1)+2)
  }else{
    A_2 = 0
  }
  C1 = 1 - f1_1/n1*(1-A_1)
  C2 = 1 - f1_2/n2*(1-A_2)
  
  p0_1_total = (1-A_1)*f1./n1 + (1-A_1)*(1-A_2)*f11/n1/n2/A_2
  p0_2_total = (1-A_2)*f.1/n2 + (1-A_1)*(1-A_2)*f11/n1/n2/A_1
  
  if(f0_1_star > 0){
    if(p0_1_total<1-C1 & f0_1_star>f0_1_p){
      p0_1 = p0_1_total/f0_1_p
      pe_1 = (1-C1-p0_1_total)/(f0_1_star-f0_1_p)
    }else if(f0_1_star>f0_1_p){
      p0_1 = (1-C1)/f0_1_p
      pe_1 = 0
    }else{
      p0_1 = 0
      pe_1 = 0
    }
  }else{
    p0_1 = 0
    pe_1 = 0
  }
  
  if(f0_2_star > 0){
    if(p0_2_total<1-C2 & f0_2_star>f0_2_p){
      p0_2 = p0_2_total/f0_2_p
      pe_2 = (1-C2-p0_2_total)/(f0_2_star-f0_2_p)
    }else if(f0_2_star>f0_2_p){
      p0_2 = (1-C2)/f0_2_p
      pe_2 = 0
    }else{
      p0_2 = 0
      pe_2 = 0
    }
  }else{
    p0_2 = 0
    pe_2 = 0
  }
  
  # Construct bootstrap assemblage------------------------
  
  I.. = which(x1>0 & x2>0)
  I.0 = which(x1>0 & x2==0)
  I0. = which(x1==0 & x2>0)
  
  W1 = (1-C1)/sum(x1/n1*(1-x1/n1)^n1)
  W2 = (1-C2)/sum(x2/n2*(1-x2/n2)^n2)
  
  p1 = x1/n1*(1 - W1*(1-x1/n1)^n1)
  p2 = x2/n2*(1 - W2*(1-x2/n2)^n2)
  
  p1_add = numeric(); p2_add = numeric()
  if(f.0 >= f.0_hat){
    p2[sample(I.0, f.0_hat)] = rep(p0_2, f.0_hat)
  }else{
    p2[I.0] = rep(p0_2, f.0)
    p1_add = rep(p0_1, f.0_hat-f.0)
    p2_add = rep(p0_2, f.0_hat-f.0)
  }
  
  if(f0. >= f0._hat){
    p1[sample(I0., f0._hat)] = rep(p0_1, f0._hat)
  }else{
    p1[I0.] = rep(p0_1, f0.)
    p1_add = c(p1_add, rep(p0_1, f0._hat-f0.))
    p2_add = c(p2_add, rep(p0_2, f0._hat-f0.))
  }
  
  p1_add = c(p1_add, rep(p0_1, f00_hat))
  p2_add = c(p2_add, rep(p0_2, f00_hat))
  
  p1_add = c(p1_add, rep(pe_1, f0_1_star-f0_1_p), rep(0, f0_2_star-f0_2_p))
  p2_add = c(p2_add, rep(0, f0_1_star-f0_1_p), rep(pe_2, f0_2_star-f0_2_p))
  
  #Basic information of distance matrix------------------------
  Fn..1. = sum(as.matrix(D)[x1>0 & x2==1, x1>0 & x2>0])
  Fn..2. = sum(as.matrix(D)[x1>0 & x2==2, x1>0 & x2>0])
  Fn..0. = ifelse(Fn..2.*sqrt(f.1)>sqrt(f.2)*Fn..1./2, (n2-1)/n2*Fn..1.^2/2/Fn..2., (n2-1)/n2*Fn..1.*(sqrt(f.1*S12.obs)-1)/2/(sqrt(f.2*S12.obs)+1))
  
  Fn.1.. = sum(as.matrix(D)[x1>0 & x2>0, x1==1 & x2>0])
  Fn.2.. = sum(as.matrix(D)[x1>0 & x2>0, x1==2 & x2>0])
  Fn.0.. = ifelse(Fn.2..*sqrt(f1.)>sqrt(f2.)*Fn.1../2, (n1-1)/n1*Fn.1..^2/2/Fn.2.., (n1-1)/n1*Fn.1..*(sqrt(f1.*S12.obs)-1)/2/(sqrt(f2.*S12.obs)+1))
  
  Fn..11 = sum(as.matrix(D)[x1>0 & x2==1, x1>0 & x2==1])
  Fn..22 = sum(as.matrix(D)[x1>0 & x2==2, x1>0 & x2==2])
  Fn..00 = ifelse(Fn..22*f.1>f.2*Fn..11/2, (n2-2)*(n2-3)/n2/(n2-1)*Fn..11^2/4/Fn..22, (n2-2)*(n2-3)/n2/(n2-1)*Fn..11*(f.1-1)/4/(f.2+1))
  
  Fn11.. = sum(as.matrix(D)[x1==1 & x2>0, x1==1 & x2>0])
  Fn22.. = sum(as.matrix(D)[x1==2 & x2>0, x1==2 & x2>0])
  Fn00.. = ifelse(Fn22..*f1.>f2.*Fn11../2, (n1-2)*(n1-3)/n1/(n1-1)*Fn11..^2/4/Fn22.., (n1-2)*(n1-3)/n1/(n1-1)*Fn11..*(f1.-1)/4/(f2.+1))
  
  Fn.1.1 = sum(as.matrix(D)[x1>0 & x2>0, x1==1 & x2==1])
  Fn.2.2 = sum(as.matrix(D)[x1>0 & x2>0, x1==2 & x2==2])
  Fn.0.0 = ifelse(Fn.2.2*sqrt(f11)>sqrt(f22)*Fn.1.1/2, (n1-1)/n1*(n2-1)/n2*Fn.1.1^2/4/Fn.2.2, (n1-1)/n1*(n2-1)/n2*Fn.1.1*(sqrt(f11*S12.obs)-1)/4/(sqrt(f22*S12.obs)+1))
  
  Fn.11. = sum(as.matrix(D)[x1>0 & x2==1, x1==1 & x2>0])
  Fn.22. = sum(as.matrix(D)[x1>0 & x2==2, x1==2 & x2>0])
  Fn.00. = ifelse(Fn.22.*sqrt(f.1*f1.)>sqrt(f.2*f2.)*Fn.11./2, (n1-1)/n1*(n2-1)/n2*Fn.11.^2/4/Fn.22., (n1-1)/n1*(n2-1)/n2*Fn.11.*(sqrt(f.1*f1.)-1)/4/(sqrt(f.2*f2.)+1))
  
  Fn.111 = sum(as.matrix(D)[x1>0 & x2==1, x1==1 & x2==1])
  Fn.222 = sum(as.matrix(D)[x1>0 & x2==2, x1==2 & x2==2])
  Fn.000 = ifelse(Fn.222*sqrt(f.1*f11)>sqrt(f.2*f22)*Fn.111/2, (n1-1)/n1*(n2-2)*(n2-3)/n2/(n2-1)*Fn.111^2/8/Fn.222, (n1-1)/n1*(n2-2)*(n2-3)/n2/(n2-1)*Fn.111*(sqrt(f.1*f11)-1)/8/(sqrt(f.2*f22)+1))
  
  Fn11.1 = sum(as.matrix(D)[x1==1 & x2>0, x1==1 & x2==1])
  Fn22.2 = sum(as.matrix(D)[x1==2 & x2>0, x1==2 & x2==2])
  Fn00.0 = ifelse(Fn22.2*sqrt(f1.*f11)>sqrt(f2.*f22)*Fn11.1/2, (n1-2)*(n1-3)/n1/(n1-1)*(n2-1)/n2*Fn11.1^2/8/Fn22.2, (n1-2)*(n1-3)/n1/(n1-1)*(n2-1)/n2*Fn11.1*(sqrt(f1.*f11)-1)/8/(sqrt(f2.*f22)+1))
  
  Fn1111 = sum(as.matrix(D)[x1==1 & x2==1, x1==1 & x2==1])
  Fn2222 = sum(as.matrix(D)[x1==2 & x2==2, x1==2 & x2==2])
  Fn0000 = ifelse(Fn2222*f11>f22*Fn1111/2, (n1-2)*(n1-3)/n1/(n1-1)*(n2-2)*(n2-3)/n2/(n2-1)*Fn1111^2/16/Fn2222, (n1-2)*(n1-3)/n1/(n1-1)*(n2-2)*(n2-3)/n2/(n2-1)*Fn1111*(f11-1)/16/(f22+1))
  
  fn..00 = sum(as.matrix(D)[x1>0 & x2==0, x1>0 & x2==0])
  fn00.. = sum(as.matrix(D)[x1==0 & x2>0, x1==0 & x2>0])
  fn..0. = sum(as.matrix(D)[x1>0 & x2==0, x1>0 & x2>0])
  fn.0.. = sum(as.matrix(D)[x1>0 & x2>0, x1==0 & x2>0])
  fn.00. = sum(as.matrix(D)[x1>0 & x2==0, x1==0 & x2>0])
  
  if(n00>1){ave.Fn0000 = (Fn0000 + max(0, Fn..00-fn..00) + max(0, Fn00..-fn00..) + 2*(max(0, Fn..0.-fn..0.)+max(0, Fn.0..-fn.0..)+max(0, Fn.00.-fn.00.)))/n00/(n00-1)}
  if(n00>0){
    ave.Fn.000 = ifelse(f.0>0, Fn.000/n00/n.0, 0)
    ave.Fn.0.0 = ifelse(S12.obs>0, Fn.0.0/n00/S12.obs, 0)
    ave.Fn00.0 = ifelse(f0.>0, Fn00.0/n00/n0., 0)
  }
  #------------------------------------------------------
  
  D_boot = matrix(0, length(p1)+length(p1_add), length(p1)+length(p1_add))
  D_boot[1:nrow(D), 1:ncol(D)] = D
  if(n00 > 1){
    D_boot[which(x1>0 & x2==0), (1:n00)+ncol(D)] = ave.Fn.000; D_boot[(1:n00)+nrow(D), which(x1>0 & x2==0)] = ave.Fn.000
    D_boot[which(x1>0 & x2>0), (1:n00)+ncol(D)] = ave.Fn.0.0; D_boot[(1:n00)+nrow(D), which(x1>0 & x2>0)] = ave.Fn.0.0
    D_boot[which(x1==0 & x2>0), (1:n00)+ncol(D)] = ave.Fn00.0; D_boot[(1:n00)+nrow(D), which(x1==0 & x2>0)] = ave.Fn00.0
    D_boot[(1:n00)+nrow(D), (1:n00)+ncol(D)] = ave.Fn0000
  }else if(n00 > 0){
    D_boot[which(x1>0 & x2==0), (1:n00)+ncol(D)] = ave.Fn.000; D_boot[(1:n00)+nrow(D), which(x1>0 & x2==0)] = ave.Fn.000
    D_boot[which(x1>0 & x2>0), (1:n00)+ncol(D)] = ave.Fn.0.0; D_boot[(1:n00)+nrow(D), which(x1>0 & x2>0)] = ave.Fn.0.0
    D_boot[which(x1==0 & x2>0), (1:n00)+ncol(D)] = ave.Fn00.0; D_boot[(1:n00)+nrow(D), which(x1==0 & x2>0)] = ave.Fn00.0
  }
  diag(D_boot) = 0
  #out =
  return(list("p_boot" = cbind(c(p1, p1_add), c(p2, p2_add)), 
              "D_boot" = D_boot));
}

#---------------------------------------------
# Shared FAD estimator between two assemblages
#---------------------------------------------
# @para. Dis is functional distance matrix of observed species.
# @para. abunDATA is an observed species-by-assemblage frequency matrix.
# **NOTE: species in functional distance matrix should be in the same ordering as those in species frequency matrix.
# @return the shared FAD 

shared.FAD=function(Dis,abunDATA){
  D=Dis;data=abunDATA;
  X1 = data[, 1]; X2 = data[, 2]
  n1 = sum(X1); n2 = sum(X2)
  index = which(X1>0 | X2>0)
  X1 = X1[index]; X2 = X2[index]
  D = D[index, index]
  D12 = D[which(X1*X2>0),which(X1*X2>0)]#D[X1>0 & X2>0, X1>0 & X2>0]
  FAD12.obs = sum(D12)
  S12.obs = sum(X1>0 & X2>0)
  f.1 = sum(X1>0 & X2==1)
  f.2 = sum(X1>0 & X2==2)
  f1. = sum(X1==1 & X2>0)
  f2. = sum(X1==2 & X2>0)
  f11 = sum(X1==1 & X2==1)
  f22 = sum(X1==2 & X2==2)
  f.0_hat = ifelse(f.2>0, (n2-1)/n2*f.1^2/2/f.2, (n2-1)/n2*f.1*(f.1-1)/2)
  f0._hat = ifelse(f2.>0, (n1-1)/n1*f1.^2/2/f2., (n1-1)/n1*f1.*(f1.-1)/2)
  f00_hat = ifelse(f22>0, (n1-1)/n1*(n2-1)/n2*f11^2/4/f22, (n1-1)/n1*(n2-1)/n2*f11*(f11-1)/4)
  
  Fn..1. = sum(as.matrix(D)[X1>0 & X2==1, X1>0 & X2>0])
  Fn..2. = sum(as.matrix(D)[X1>0 & X2==2, X1>0 & X2>0])
  Fn..0. = ifelse(Fn..2.*sqrt(f.1)>sqrt(f.2)*Fn..1./2, (n2-1)/n2*Fn..1.^2/2/Fn..2., (n2-1)/n2*Fn..1.*(sqrt(f.1*S12.obs)-1)/2/(sqrt(f.2*S12.obs)+1))
  
  Fn.1.. = sum(as.matrix(D)[X1>0 & X2>0, X1==1 & X2>0])
  Fn.2.. = sum(as.matrix(D)[X1>0 & X2>0, X1==2 & X2>0])
  Fn.0.. = ifelse(Fn.2..*sqrt(f1.)>sqrt(f2.)*Fn.1../2, (n1-1)/n1*Fn.1..^2/2/Fn.2.., (n1-1)/n1*Fn.1..*(sqrt(f1.*S12.obs)-1)/2/(sqrt(f2.*S12.obs)+1))
  
  Fn..11 = sum(as.matrix(D)[X1>0 & X2==1, X1>0 & X2==1])
  Fn..22 = sum(as.matrix(D)[X1>0 & X2==2, X1>0 & X2==2])
  Fn..00 = ifelse(Fn..22*f.1>f.2*Fn..11/2, (n2-2)*(n2-3)/n2/(n2-1)*Fn..11^2/4/Fn..22, (n2-2)*(n2-3)/n2/(n2-1)*Fn..11*(f.1-1)/4/(f.2+1))
  
  Fn11.. = sum(as.matrix(D)[X1==1 & X2>0, X1==1 & X2>0])
  Fn22.. = sum(as.matrix(D)[X1==2 & X2>0, X1==2 & X2>0])
  Fn00.. = ifelse(Fn22..*f1.>f2.*Fn11../2, (n1-2)*(n1-3)/n1/(n1-1)*Fn11..^2/4/Fn22.., (n1-2)*(n1-3)/n1/(n1-1)*Fn11..*(f1.-1)/4/(f2.+1))
  
  Fn.1.1 = sum(as.matrix(D)[X1>0 & X2>0, X1==1 & X2==1])
  Fn.2.2 = sum(as.matrix(D)[X1>0 & X2>0, X1==2 & X2==2])
  Fn.0.0 = ifelse(Fn.2.2*sqrt(f11)>sqrt(f22)*Fn.1.1/2, (n1-1)/n1*(n2-1)/n2*Fn.1.1^2/4/Fn.2.2, (n1-1)/n1*(n2-1)/n2*Fn.1.1*(sqrt(f11*S12.obs)-1)/4/(sqrt(f22*S12.obs)+1))
  
  Fn.11. = sum(as.matrix(D)[X1>0 & X2==1, X1==1 & X2>0])
  Fn.22. = sum(as.matrix(D)[X1>0 & X2==2, X1==2 & X2>0])
  Fn.00. = ifelse(Fn.22.*sqrt(f.1*f1.)>sqrt(f.2*f2.)*Fn.11./2, (n1-1)/n1*(n2-1)/n2*Fn.11.^2/4/Fn.22., (n1-1)/n1*(n2-1)/n2*Fn.11.*(sqrt(f.1*f1.)-1)/4/(sqrt(f.2*f2.)+1))
  
  Fn.111 = sum(as.matrix(D)[X1>0 & X2==1, X1==1 & X2==1])
  Fn.222 = sum(as.matrix(D)[X1>0 & X2==2, X1==2 & X2==2])
  Fn.000 = ifelse(Fn.222*sqrt(f.1*f11)>sqrt(f.2*f22)*Fn.111/2, (n1-1)/n1*(n2-2)*(n2-3)/n2/(n2-1)*Fn.111^2/8/Fn.222, (n1-1)/n1*(n2-2)*(n2-3)/n2/(n2-1)*Fn.111*(sqrt(f.1*f11)-1)/8/(sqrt(f.2*f22)+1))
  
  Fn11.1 = sum(as.matrix(D)[X1==1 & X2>0, X1==1 & X2==1])
  Fn22.2 = sum(as.matrix(D)[X1==2 & X2>0, X1==2 & X2==2])
  Fn00.0 = ifelse(Fn22.2*sqrt(f1.*f11)>sqrt(f2.*f22)*Fn11.1/2, (n1-2)*(n1-3)/n1/(n1-1)*(n2-1)/n2*Fn11.1^2/8/Fn22.2, (n1-2)*(n1-3)/n1/(n1-1)*(n2-1)/n2*Fn11.1*(sqrt(f1.*f11)-1)/8/(sqrt(f2.*f22)+1))
  
  Fn1111 = sum(as.matrix(D)[X1==1 & X2==1, X1==1 & X2==1])
  Fn2222 = sum(as.matrix(D)[X1==2 & X2==2, X1==2 & X2==2])
  Fn0000 = ifelse(Fn2222*f11>f22*Fn1111/2, (n1-2)*(n1-3)/n1/(n1-1)*(n2-2)*(n2-3)/n2/(n2-1)*Fn1111^2/16/Fn2222, (n1-2)*(n1-3)/n1/(n1-1)*(n2-2)*(n2-3)/n2/(n2-1)*Fn1111*(f11-1)/16/(f22+1))
  
  FAD12 = FAD12.obs + Fn..00 + Fn00.. + Fn0000 + 2*(Fn..0. + Fn.0.. + Fn.0.0 + Fn.00. + Fn.000 + Fn00.0)
  ###
  return(FAD12)
}


#----------------------------------------------
# Shared FAD estimator between two assemblages
#----------------------------------------------
# Two functions "shared.FAD" and "EstBootCommTwoFun" should be run before running this function.
# @para. Dis is a functional distance matrix of observed species.
# @para. abunDATA is an observed species-by-assemblage frequency matrix.
# **NOTE: species in functional distance matrix should be in the same ordering as those in species frequency matrix.
# @return shared FAD estimator and the associated 95% confidence interval.

Good_Turing_shared.FD=function(Dis,abunDATA){
  da1=abunDATA[,1];da2=abunDATA[,2];
  n1=sum(da1);n2=sum(da2);
  I=which(da1*da2>0);
  obsfd12=sum(Dis[I,I]);
  sfd12=shared.FAD(Dis,abunDATA);
  ufd12=max(0,sfd12-obsfd12);
  o=EstBootCommTwoFun(Dis,abunDATA)
  boots=50;est=numeric(boots)
  for(i in 1:boots){
    data1=rmultinom(1,size=n1,prob=o$p_boot[,1])
    data2=rmultinom(1,size=n1,prob=o$p_boot[,2])
    D.boot=o$D_boot;
    est[i]=shared.FAD(D.boot,cbind(data1,data2))
  }
  V=var(est);
  R=exp(1.96*(log(1+V/(sfd12)^2))^(1/2));
  sfd12.CI=c(obsfd12+ufd12/R,obsfd12+ufd12*R) 
  return(list(shared_FD=sfd12,CI=sfd12.CI))
}
########
