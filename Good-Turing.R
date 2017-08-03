#*******************************************************************************
#********************************************************************************
  ## R scipts "Good-Turing" for Chao et al. (2017) Ecology paper. 
  ## Please cite the following Chao et al. (2017) paper to use these R scripts:
  ## Chao, A, Chiu, C.-H., Colwell, R.K., Magnago, L.F.S., Chazdon,R.L. and Gotelli, N. J. (2017) 
  ## Deciphering the Enigma of Undetected Species, Phylogenetic, and Functional Diversity Based on Good-Turing Theory. Ecology.
  #
  # The following R scripts include 6 parts:
  # (1). Script for estimating species richness for each individual assemblage (Table 2(a) of Chao et al. 2017 paper)
  # (2). Script for estimating shared species richness between two assemblages (Table 2(b) of Chao et al. 2017 paper
  # (3). Script for estimating Faith PD for each individual assemblage (Table 3(a) of Chao et al. 2017 paper)
  # (4). Script for estimating shared PD between two assemblages (Table 3(b) of Chao et al. 2017 paper)
  # (5). Script for estimating FAD for each individual assemblage (Table 4(a) of Chao et al. 2017 paper)
  # (6). Script for estimating shared FAD between two assemblages (Table 4(b) of Chao et al. 2017 paper)
# 
# NOTE: The packages "ade4", "phytools", "ape" and "knitr" must be installed and loaded before running the scripts. 
# These four packages can be downloaded from CRAN.
#
#*******************************************************************************
#*******************************************************************************



####################################################################################
#
# (1).  R function for estimating species richness in one assemblage: Richness(data)
#
####################################################################################

#' Estimating species richness in one assemblage
#' Richness(data) is a function of obtaining estimator of species richness in one assemblage based on abundance data.
#' @param data is an observed species-by-assemblage frequency matrix. The number of assemblages is allowed to be any positive integer.
#' @return the Chao1 species richness estimate and its 95% confidence interval.

Richness <- function(data){
  Good_Turing_Richness=function(data){
    da=data;
    n=sum(da);obs=sum(da>0);
    f1=sum(da==1); f2=sum(da==2);
    us=round(ifelse(f1>0,f1^2/(2*f2),f1*(f1-1)/2 ));
    est=us+obs;
    ###95% C.I. of S
    fi=sapply(1:max(da),function(k) sum(da==k));
    k = (n-1)/n
    V = ifelse(fi[2] != 0, fi[2]*(0.25*k^2*(fi[1]/fi[2])^4+k^2*(fi[1]/fi[2])^3+0.5*k*(fi[1]/fi[2])^2),
               k*fi[1]*(fi[1]-1)/2+k^2*fi[1]*(2*fi[1]-1)/4-k^2*fi[1]^4/(4*est))
    R=exp(1.96*(log(1+V/(us)^2))^(1/2));
    est.CI=round(c(obs+us/R,obs+us*R));
    c("Sample size" = n, "f1" = f1, "f2" = f2, "Observed richness" = obs,
      "Undetected richness" = us, "Chao1 richness" = est,
      "95% conf. interval" = paste0("(", est.CI[1], ", ", est.CI[2], ")"))
  }
  ans = t(apply(data, 2, function(K) Good_Turing_Richness(K)))
  return(kable(ans, align = 'c', digits = 0))
}



##########################################################################################################
#
# (2).  R function for estimating shared species richness between two assemblages: Shared_richness(data)
#
##########################################################################################################

#' Estimating shared species richness between two assemblages
#' Shared_richness(data) is a function of obtaining estimator of shared species richness between two assemblages based on abundance data.
#' @param data is an observed species-by-assemblage frequency matrix. The number of assemblages must be two.
#' @return the Chao1-shared richness estimate and its 95% confidence interval.

Shared_richness <- function(data){
  da1 = data[,1];da2=data[,2];n1 = sum(da1);n2 = sum(da2)
  ##shared richness estimaiton
  I = which(da1*da2>0);sda1 = da1[I];sda2 = da2[I];obs12 = length(I)
  f11 = sum(sda1==1 & sda2==1);f22 = sum(sda1==2 & sda2==2);f22 = max(1,f22)
  f1_ = sum(sda1==1); f2_ = sum(sda1==2);f2_ = max(1,f2_)
  f_1 = sum(sda2==1); f_2 = sum(sda2==2);f_2 = max(1,f_2)
  k1 = (n1-1)/n1; k2 = (n2-1)/n2
  f0_ = k1*f1_^2/(2*f2_);f_0 = k2*f_1^2/(2*f_2)
  f00 = k1*k2*f11^2/(4*f22)
  us12 = f00+f0_+f_0
  s12=obs12+us12
  ###95%CI of shared richness
  fii=c(obs12,f11,f22,f1_,f2_,f_1,f_2)
  df=c(1,k1*k2*f11/(2*f22),-k1*k2*(f11/f22)^2/4,k1*f1_/f2_,-k1*(f1_/f2_)^2/2,k2*f_1/f_2,-k2*(f_1/f_2)^2/2)
  COV=matrix(0,7,7);COV[7,7]=fii[7]*(1-fii[7]/s12)
  for(i in 1:6){
    COV[i,i]=fii[i]*(1-fii[i]/s12)
    for(j in (i+1):7){
      COV[i,j]=-fii[i]*fii[j]/s12;COV[i,j]=COV[j,i]}
  }
  V=t(df)%*%COV%*%df
  R=exp(1.96*(log(1+V/(us12)^2))^(1/2))
  s12.CI=round(c(obs12+us12/R,obs12+us12*R))
  ans = c("Observed" = obs12, "f+1" = f_1, "f+2" = f_2, "f1+" = f1_, "f2+" = f2_,
          f11 = f11, f22 = f22, "f+0" = round(f_0), "f0+" = round(f0_), f00 = round(f00), 
          "Undetected" = round(us12), "Chao1 Shared" = round(s12),
          "95% conf. interval" = paste0("(", s12.CI[1], ", ", s12.CI[2], ")"))
  return(kable(t(as.matrix(ans)), align = 'c', digits = 0))
}



#################################################################################
#
# (3).  R function for estimating Faith PD in one assemblage: PD(data, tree)
#
#################################################################################

#' Estimating Faith PD in one assemblage
#' PD(data, tree) is a function of obtaining estimator of Faith PD based on abundance data.
#' @param data is an observed species-by-assemblage frequency matrix. The number of assemblages is allowed to be any positive integer.
#' @param tree is in the Newick format of the phylogenetic tree spanned by the observed species in data.
#' @return the Chao1-PD estimate and its 95% confidence interval.

PD <- function(data, tree) {
  Good_Turing_PD=function(tree,data){
    data1 <- data[data>0]
    tip <- tree$tip.label[!(tree$tip.label %in% names(data1))]
    subtree <- drop.tip(tree,tip)
    phyloData <- newick2phylog(write.tree(subtree))
    
    #phyloData <- newick2phylog(write.tree(tree)); 
    abun <- data1[names(phyloData$leaves)];
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
    n=sum(data)
    k = (n-1)/n
    upd=ifelse( (2*g2*f1 <= g1*f2), k*g1*(f1-1)/(2*(f2+1)),k*(g1^2)/(2*g2))
    pd=obspd+upd;
    V = ifelse(g2 != 0, g2*(0.25*k^2*(g1/g2)^4+k^2*(g1/g2)^3+0.5*k*(g1/g2)^2),
               k*g1*(g1-1)/2+k^2*g1*(2*g1-1)/4-k^2*g1^4/(4*est))
    R=exp(1.96*(log(1+V/(upd)^2))^(1/2));
    pd.CI=round(c(obspd+upd/R,obspd+upd*R))
    #return(matrix(c(sum(abunDATA), g1, g2, obspd, upd, pd, pd.CI), nrow = 1))
    c("Sample size" = sum(data), g1 = g1, g2 = g2, "Observed" = obspd, 
      "Undetected" = round(upd), "Chao1-PD" = round(pd), 
      "95% conf. interval" = paste0("(", pd.CI[1], ", ", pd.CI[2], ")"))
  }
  N <- ncol(data)
  ans <- t(apply(data, 2, function(K) Good_Turing_PD(tree, K)))
  kable(ans, align = 'c')
}




############################################################################################
#
# (4).  R function for estimating shared PD between two assemblages: Shared_PD(data, tree)
#
############################################################################################

#' Shared PD between two assemblages
#' Shared_PD(data, tree) is a function of obtaining estimator of shared PD between two assemblages based on abundance data.
#' @param data is an observed species-by-assemblage frequency matrix. The number of assemblages must be two.
#' @param tree is in the Newick format of the phylogenetic tree spanned by the observed species in data.
#' @return the Chao1-PD-shared estimate and its 95% confidence interval.

Shared_PD=function(data, tree){
  tip <- tree$tip.label[!(tree$tip.label %in% rownames(data))]
  subtree <- drop.tip(tree,tip)
  phyloData <- newick2phylog(write.tree(subtree))
  abun <- as.matrix(data[names(phyloData$leaves), ])
  nodenames=c(names(phyloData$leaves),names(phyloData$nodes));
  
  M=matrix(0,nrow=length(phyloData$leaves),ncol=length(nodenames),dimnames=list(names(phyloData$leaves),nodenames))
  for(i in 1:length(phyloData$leaves)){
    M[i,][unlist(phyloData$paths[i])]=rep(1,length(unlist(phyloData$paths[i])))
  }
  phyloAbun=matrix(0,ncol=ncol(data),nrow=length(nodenames),dimnames=list(nodenames,colnames(data)))
  for(i in 1:ncol(data)){
    phyloAbun[,i]=abun[,i]%*%M;}
  BL=c(phyloData$leaves,phyloData$nodes)
  ######
  pda1=phyloAbun[,1];obspd1=sum(BL[pda1>0])
  pda2=phyloAbun[,2];obspd2=sum(BL[pda2>0])
  ###### CI of pd12
  I=which(pda1*pda2>0);sda1=pda1[I];sda2=pda2[I];sBL=BL[I];obspd12=sum(sBL)
  g11=sum(sBL[sda1==1 & sda2==1]);g22=sum(sBL[sda1==2 & sda2==2])
  g1_=sum(sBL[sda1==1]);g2_=sum(sBL[sda1==2])
  g_1=sum(sBL[sda2==1]);g_2=sum(sBL[sda2==2])
  n1 = sum(data[, 1]);n2 = sum(data[, 2])
  k1 = (n1-1)/n1; k2 = (n2-1)/n2
  gp0h = k2*g_1^2/(2*g_2)
  g0ph = k1*g1_^2/(2*g2_)
  g00h = k1*k2*g11^2/(4*g22)
  upd12=gp0h+g0ph+g00h;spd12=obspd12+upd12;
  gii=c(obspd12,g11,g22,g1_,g2_,g_1,g_2);
  dg=c(1,k1*k2*g11/(2*g22),-k1*k2*(g11/g22)^2/4,k1*g1_/g2_,-k1*(g1_/g2_)^2/2,k2*g_1/g_2,-k2*(g_1/g_2)^2/2);
  COV=matrix(0,7,7);COV[7,7]=gii[7]*(1-gii[7]/spd12);
  for(i in 1:6){
    COV[i,i]=gii[i]*(1-gii[i]/spd12);
    for(j in (i+1):7){
      COV[i,j]=-gii[i]*gii[j]/spd12;COV[i,j]=COV[j,i];}
  }
  V=t(dg)%*%COV%*%dg;
  R=exp(1.96*(log(1+V/(upd12)^2))^(1/2));
  spd12.CI=c(obspd12+upd12/R,obspd12+upd12*R) 
  final <- c("Observed" = obspd12, "g+1" = g_1, "g+2" = g_2, "g1+" = g1_, "g2+" = g2_, g11 = g11, 
             g22 = g22, "g+0" = round(gp0h), "g0+" = round(g0ph), g00 = round(g00h), "Undetected" = round(upd12), 
             "Chao1-PD-shared" = round(spd12), "95% conf. interval" = paste0("(", round(spd12.CI[1]), ", ", round(spd12.CI[2]), ")"))
  
  kable(t(final))
}




#################################################################################
#
# (5).  R function for estimating FAD in one assemblage: FAD(data, dis_matrix)
#
#################################################################################

#' Estimating FAD in one assemblage
#' FAD(data, dis_matrix) is a function of obtaining estimator of FAD based on abundance data.
#' @param data is an observed species-by-assemblage frequency matrix. The number of assemblages is allowed to be any positive integer.
#' @param dis_matrix is functional distance matrix of pairs of observed species.
#' @return the Chao1-FAD estimate and its 95% confidence interval.

FAD <- function(data, dis_matrix) {
  Good_Turing_FD=function(dis_matrix,data){
    da=c(data);d=mean(dis_matrix);
    I=which(da>0);D=dis_matrix[I,I];ab=da[I];obsfd=sum(D);
    F1_=sum(D[ab==1,]);F2_=sum(D[ab==2,]);F2_=max(d,F2_);
    F_1=sum(D[,ab==1]);F_2=sum(D[,ab==2]);F_2=max(d,F_2);
    F11=sum(D[ab==1,ab==1]);F22=sum(D[ab==2,ab==2]);F22=max(d,F22)
    n=sum(data)
    k=(n-1)/n;k_star = (n-2)*(n-3)/(n*(n-1))
    Fp0h = k*F_1^2/(2*F_2)
    F0ph = k*F1_^2/(2*F2_)
    F00h = k_star*F11^2/(4*F22)
    ufd=Fp0h+F0ph+F00h;
    fd=obsfd+ufd;
    
    Fii=c(obsfd,F11,F22,F1_,F2_,F_1,F_2);
    dF=c(1,k_star*F11/(2*F22),-k_star*(F11/F22)^2/4,k*F1_/F2_,-k*(F1_/F2_)^2/2,k*F_1/F_2,-k*(F_1/F_2)^2/2);
    COV=matrix(0,7,7);COV[7,7]=Fii[7]*(1-Fii[7]/fd);
    for(i in 1:6){
      COV[i,i]=Fii[i]*(1-Fii[i]/fd);
      for(j in (i+1):7){
        COV[i,j]=-Fii[i]*Fii[j]/fd;COV[i,j]=COV[j,i];}
    }
    V=t(dF)%*%COV%*%dF;
    R=exp(1.96*(log(1+V/(ufd)^2))^(1/2));
    fd.CI=c(obsfd+ufd/R,obsfd+ufd*R) 
    
    final <- c("Observed" = round(obsfd), "F+1" = round(F_1), "F+2" = round(F_2), 
               "F11" = round(F11), "F22" = round(F22),"F+0" = round(Fp0h), "F00" = round(F00h),
               "Undetected" = round(ufd), "Chao1-FAD" = round(fd), 
               "95% conf. interval" = paste0("(", round(fd.CI[1]), ", ", round(fd.CI[2]), ")"))
    final
  }
  ans <- apply(data,2, function(K) {
    #names(da) <- rownames(data)
    names(K) <- rownames(data)
    Good_Turing_FD(as.matrix(dis), K)
  })
  
  kable(t(as.data.frame(ans)))
}



##################################################################################################
#
# (6).  R function for estimating shared FAD between two assemblages: Shared_FAD(data, dis_matrix)
#
##################################################################################################

#' Estimating shared FAD between two assemblages
#' Shared_FAD(data, dis_matrix) is a function of obtaining estimator of shared FAD between two assemblages based on abundance data.
#' @param data is an observed species-by-assemblage frequency matrix. The number of assemblages must be two.
#' @param dis_matrix is functional distance matrix of all pairs of observed species.
#' @return the Chao1-FAD-shared estimate and its 95% confidence interval.

Shared_FAD=function(data,dis_matrix){
  
  EstBootCommTwoFun = function(dis_matrix,data){
    D = dis_matrix
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
  shared.FAD=function(data, dis_matrix){
    D = dis_matrix
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
    final = data.frame(F..00 = Fn..00, F00.. = Fn00.., F0000 = Fn0000, F..0. = Fn..0., F.0.. = Fn.0.., 
                       F.0.0 = Fn.0.0, F.00. = Fn.00., F.000 = Fn.000, F00.0 = Fn00.0, FAD12 = FAD12)
    return(final)
  }
  da1=data[,1];da2=data[,2];
  n1=sum(da1);n2=sum(da2);
  I=which(da1*da2>0);
  obsfd12=sum(dis[I,I]);
  share = shared.FAD(data, dis_matrix)
  sfd12=share$FAD12;
  ufd12=max(0,sfd12-obsfd12);
  o=EstBootCommTwoFun(dis_matrix, data)
  boots=200;est=numeric(boots)
  for(i in 1:boots){
    data1=rmultinom(1,size=n1,prob=o$p_boot[,1])
    data2=rmultinom(1,size=n1,prob=o$p_boot[,2])
    D.boot=o$D_boot;
    est[i]=shared.FAD(cbind(data1,data2), D.boot)$FAD12
  }
  V=var(est);
  R=exp(1.96*(log(1+V/(ufd12)^2))^(1/2));
  sfd12.CI=c(obsfd12+ufd12/R,obsfd12+ufd12*R) 
  
  final = data.frame(obs = round(obsfd12), Fpp00 = round(share$F..00), Fpp0p = round(share$F..0.), 
                     F00pp = round(share$F00..), Fp0pp = round(share$F.0..), Fp00p = round(share$F.00.), 
                     Fp000 = round(share$F.000), Fp0p0 = round(share$F.0.0), F00p0 = round(share$F00.0), 
                     F0000 = round(share$F0000), und = round(ufd12), fad = round(sfd12), 
                     "95% conf. interval" = paste0("(", round(sfd12.CI[1]), ", ", round(sfd12.CI[2]), ")"))
  colnames(final) = c("Observed", "F(++)(00)", "F(++)(0+)", "F(00)(++)", 
                      "F(+0)(++)", "F(+0)(0+)", "F(+0)(00)", 
                      "F(+0)(+0)", "F(00)(+0)", "F(00)(00)", 
                      "Undetected", "Chao1-FAD-shared", paste("95%", "CI"))
  kable(final)
}
