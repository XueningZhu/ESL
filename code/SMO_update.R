Update.alphas<-function(alpha1_ind,alpha2_ind,alpha,y,E,C,eps,ker_mat)
{
  yy=y[c(alpha1_ind,alpha2_ind)]
  ee=E[c(alpha1_ind,alpha2_ind)]
  if (yy[1]!=yy[2])
  {
    L=max(0,alpha[alpha2_ind]-alpha[alpha1_ind])
    H=min(C,C+alpha[alpha2_ind]-alpha[alpha1_ind])
  }
  else
  {
    L=max(0,alpha[alpha2_ind]+alpha[alpha1_ind]-C)
    H=min(C,alpha[alpha2_ind]+alpha[alpha1_ind])
  }
  kk=ker_mat[alpha1_ind,alpha1_ind]+ker_mat[alpha2_ind,alpha2_ind]-
    2*ker_mat[alpha1_ind,alpha2_ind]
  
  a2=alpha[alpha2_ind]+yy[2]*(ee[1]-ee[2])/kk
  alpha2=a2
  if (a2>H)
  {
    alpha2=H
  }
  if (a2<L)
  {
    alpha2=L
  }
  alpha1=alpha[alpha1_ind]+yy[1]*yy[2]*(alpha[alpha2_ind]-alpha2)
  
  alpha_old=alpha
  alpha_new=c(alpha1,alpha2)
  if (abs(ee[1]-ee[2])<eps|abs(L-H)<eps|sqrt(mean(sum((alpha_old[c(alpha1_ind,alpha2_ind)]-alpha_new)^2)))<eps)
  {
    return(list(changed=F,alpha_new=alpha_new))
  }
  return(list(changed=T,alpha_new=alpha_new))
}

Update<-function(y,alpha,C,E,ker_mat,eps,b)
{
  alpha_changed=1
  examineAll=T
  iter=1
  count=0
  alpha_count=0
  
  while((examineAll|alpha_changed>0)&iter<=1000)
  {
    show(iter)
    iter=iter+1
    if (examineAll)
    {
      alpha1_inds=1:length(alpha)
    }
    else
    {
      alpha1_inds=is.support(alpha,C,eps)
    }
    alpha1_inds=sample(alpha1_inds,length(alpha1_inds))
    alpha_changed=0
    alpha_oo=alpha
    ######################################
    for (i in 1:length(alpha1_inds))
    {
      alpha1_ind=alpha1_inds[i]
      alpha2_ind=second.alpha(E,alpha1_ind)
      alphas_up=Update.alphas(alpha1_ind,alpha2_ind,alpha,y,E,C,eps,ker_mat)
      if (alphas_up$changed)
      {
        alpha_changed=alpha_changed+1
      }
      else
      {
        for (j in sample(setdiff(1:length(alpha),alpha1_ind),length(alpha)-1))
        {
          alpha2_ind=j
          alphas_up=Update.alphas(alpha1_ind,alpha2_ind,alpha,y,E,C,eps,ker_mat)
          if (alphas_up$changed)
          {
            alpha_changed=alpha_changed+1
            break
          }
        }
        if (!alphas_up$changed)
          next
      }
      
      b=update.b(E,y,b,ker_mat,alpha1_ind,alpha2_ind,
                 alpha_old=alpha[c(alpha1_ind,alpha2_ind)],alphas_up$alpha_new,eps,C)
      
      alpha[c(alpha1_ind,alpha2_ind)]=alphas_up$alpha_new
      gg=sapply(1:length(alpha),g,alpha=alpha,y=y,
                ker_mat=ker_mat,b=b)
      E=update.E(E,y,alpha,ker_mat,b)

      #     cat("delta alpha:",sqrt(mean(sum((alpha_old[c(alpha1_ind,alpha2_ind)]-alpha_new)^2))),"\n")
      
      #     cat("alpha changed:",alpha_changed,"\n")
      #     cat("\n")
    }
    #############################################
    if (examineAll)
    {
      cat("examineAll=T","\n")
      cat("alpha changed:",alpha_changed,"\n")
      cat("delta alpha:",sqrt(mean(sum((alpha_oo-alpha)^2))),"\n")
    }
    ############################################
    if (examineAll)
    {
      examineAll = F
      if (alpha_changed<=10&count>=5)
      {
        count=count+1
        break
      }
    }
    else 
    {
      cat("examineAll=F","\n")
      cat("alpha changed:",alpha_changed,"\n")
      cat("delta alpha:",sqrt(mean(sum((alpha_oo-alpha)^2))),"\n")
      
      if (alpha_changed <= 10)
      {
        alpha_count=alpha_count+1
        if (alpha_count>=5|alpha_changed==0)
        {
          alpha_count=0
          examineAll = T
        }
        cat("alpha isn't changed!","\n")
      }
    }
  }
  cat("iter:",iter,"count:",count,"\n")
  return(list(Alpha=alpha,B=b,E=E,Fit=gg))
}