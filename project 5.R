load("ncaaf_2010.rdata")

library(tidyverse)


load("games2010.Rdata")

gendata = function(q){
  awaygames = summarise(group_by(q,awayteam),num=n(),numwins=sum(!homewin))
  homegames = summarise(group_by(q,hometeam),num=n(),numwins=sum(homewin))
  totalgames = full_join(awaygames,homegames,by=c(awayteam="hometeam"))
  totalgames$num.x[is.na(totalgames$num.x)]=0
  totalgames$num.y[is.na(totalgames$num.y)]=0
  totalgames$numwins.x[is.na(totalgames$numwins.x)]=0
  totalgames$numwins.y[is.na(totalgames$numwins.y)]=0
  totalgames$wins = totalgames$numwins.x+totalgames$numwins.y
  totalgames$losses=totalgames$num.x+totalgames$num.y-totalgames$numwins.x-totalgames$numwins.y
  
  
  homespread = summarise(group_by(q,hometeam),o1=first(awayteam),o2=nth(awayteam,2),o3=nth(awayteam,3),o4=nth(awayteam,4),o5=nth(awayteam,5),o6=nth(awayteam,6),o7=nth(awayteam,7),o8=nth(awayteam,8),o9=nth(awayteam,9)) 
  awayspread = summarise(group_by(q,awayteam),o1=first(hometeam),o2=nth(hometeam,2),o3=nth(hometeam,3),o4=nth(hometeam,4),o5=nth(hometeam,5),o6=nth(hometeam,6),o7=nth(hometeam,7),o8=nth(hometeam,8),o9=nth(hometeam,9)) 
  totalspread = full_join(awayspread,homespread,by=c(awayteam="hometeam"))
  
  opps=list()
  home=list()
  for(i in 1:length(totalgames$awayteam)){
    tmp=as.character(totalspread[i,])
    tmp=tmp[-1]
    h=c(F,F,F,F,F,F,F,F,F,T,T,T,T,T,T,T,T,T)
    h=h[!is.na(tmp)]
    tmp=tmp[!is.na(tmp)]
    
    opps[[i]]=tmp
    home[[i]]=h
  }
  colnames(totalgames)=c("awayteam","numaway","numawaywins","numhome","numhomewins","wins","losses")
  
  totalgames$opponents=opps
  totalgames$home=home
  return(totalgames)
}


ll1 = function(param,q,totalgames){
  lambda=param[1]
  
  C <- matrix(integer(length(totalgames$awayteam)^2),nrow=length(totalgames$awayteam))
  rownames(C)=totalgames$awayteam
  colnames(C)=totalgames$awayteam
  b <- numeric(length(totalgames$awayteam))
  for(i in 1:length(totalgames$awayteam)){
    C[i,i]=2+totalgames$wins[i]+totalgames$losses[i]
    for(j in totalgames$opponents[[i]]){
      C[i,j]=C[i,j]-1
    }
    b[i]=1+(totalgames$wins[i]-totalgames$losses[i])/2
  }
  r <- solve(C,b)
  
  l=0
  for(i in 1:nrow(q)){
    ra=r[q$awayteam[i]]
    rh=r[q$hometeam[i]]
    l = l+lambda*as.numeric(q$homewin[i])*rh+lambda*(1-as.numeric(q$homewin[i]))*ra-log(exp(lambda*rh)+exp(lambda*ra))
  }
  return(l)
  
}

ll1min<-function(b,q,totalgames){
  -ll1(b,q,totalgames)
}

totalgames=gendata(games)
(result <- nlm(ll1min,1,hessian=TRUE,q=games,totalgames=totalgames))

mse = function(param,q,totalgames){
  lambda=param[1]
  
  C <- matrix(integer(length(totalgames$awayteam)^2),nrow=length(totalgames$awayteam))
  rownames(C)=totalgames$awayteam
  colnames(C)=totalgames$awayteam
  b <- numeric(length(totalgames$awayteam))
  for(i in 1:length(totalgames$awayteam)){
    C[i,i]=2+totalgames$wins[i]+totalgames$losses[i]
    for(j in totalgames$opponents[[i]]){
      C[i,j]=C[i,j]-1
    }
    b[i]=1+(totalgames$wins[i]-totalgames$losses[i])/2
  }
  r <- solve(C,b)
  
  mse=0
  for(i in 1:nrow(q)){
    ra=r[q$awayteam[i]]
    rh=r[q$hometeam[i]]
    p=exp(rh*lambda)/(exp(lambda*rh)+exp(lambda*ra))
    mse=mse+(as.numeric(q$homewin[i])-p)^2
  }
  names(mse) = NULL
  return(mse/nrow(q))
  
}
set.seed(20)
folds=sample(1:5,nrow(games),replace=T)
testmse=0
for(i in 1:5){
  f=which(folds==i)
  tgames=gendata(games[-f,])
  result=nlm(ll1min,1,hessian=TRUE,q=games[-f,],totalgames=tgames)
  testmse=testmse+mse(result$estimate,games[f,],tgames)
  
}
testmse=testmse/5




ll2 = function(param,q,totalgames){
  lambda=param[1]
  gamma=param[2]
  
  C <- matrix(integer(length(totalgames$awayteam)^2),nrow=length(totalgames$awayteam))
  rownames(C)=totalgames$awayteam
  colnames(C)=totalgames$awayteam
  b <- numeric(length(totalgames$awayteam))
  for(i in 1:length(totalgames$awayteam)){
    C[i,i]=2+totalgames$wins[i]+totalgames$losses[i]
    for(j in totalgames$opponents[[i]]){
      C[i,j]=C[i,j]-(1-gamma)*2
      
    }
    
    b[i]=1+(gamma*totalgames$wins[i]-(1-gamma)*totalgames$losses[i])
    
  }
  r <- solve(C,b)
  
  l=0
  for(i in 1:nrow(q)){
    ra=r[q$awayteam[i]]
    rh=r[q$hometeam[i]]
    l = l+lambda*as.numeric(q$homewin[i])*rh+lambda*(1-as.numeric(q$homewin[i]))*ra-log(exp(lambda*rh)+exp(lambda*ra))
  }
  return(l)
  
}

llmin2<-function(b,q,totalgames){
  -ll2(b,q,totalgames)
}

(result <- nlm(llmin2,c(1,0.5),hessian=TRUE,q=games,totalgames=totalgames))

mse2 = function(param,q,totalgames){
  lambda=param[1]
  gamma=param[2]
  
  C <- matrix(integer(length(totalgames$awayteam)^2),nrow=length(totalgames$awayteam))
  rownames(C)=totalgames$awayteam
  colnames(C)=totalgames$awayteam
  b <- numeric(length(totalgames$awayteam))
  for(i in 1:length(totalgames$awayteam)){
    C[i,i]=2+totalgames$wins[i]+totalgames$losses[i]
    for(j in totalgames$opponents[[i]]){
      C[i,j]=C[i,j]-(1-gamma)*2
    }
    
    b[i]=1+(gamma*totalgames$wins[i]-(1-gamma)*totalgames$losses[i])
    #b[i]=1+(gamma*totalgames$wins[i]-gamma*totalgames$losses[i])
  }
  r <- solve(C,b)
  
  mse=0
  for(i in 1:nrow(q)){
    ra=r[q$awayteam[i]]
    rh=r[q$hometeam[i]]
    p=exp(rh*lambda)/(exp(lambda*rh)+exp(lambda*ra))
    mse=mse+(as.numeric(q$homewin[i])-p)^2
  }
  names(mse) = NULL
  return(mse/nrow(q))
}
#set.seed(20)
#folds=sample(1:5,nrow(q),replace=T)
testmse2=0
for(i in 1:5){
  f=which(folds==i)
  tgames=gendata(games[-f,])
  result=nlm(llmin2,c(1,0.5),hessian=TRUE,q=games[-f,],totalgames=tgames)
  testmse2=testmse2+mse2(result$estimate,games[f,],tgames)
  
}
testmse2=testmse2/5





rpi = function(param,q,totalgames){
  lambda=param[1]
  
  wp = totalgames$wins/(totalgames$wins+totalgames$losses)
  r = numeric(length(wp))
  names(r)=totalgames$awayteam
  for(i in 1:length(totalgames$awayteam)){
    
    w=0
    l=0
    w2=0
    l2=0
    for(j in totalgames$opponents[[i]]){
      jind=which(j==totalgames$awayteam)
      w = w + totalgames$wins[jind]
      l = l + totalgames$losses[jind]
      
      ind=which(totalgames$awayteam[i]==q$awayteam&j==q$hometeam)
      w=w-sum(q$homewin[ind])
      l=l-(length(q$homewin[ind])-sum(q$homewin[ind]))
      
      ind=which(totalgames$awayteam[i]==q$hometeam&j==q$awayteam)
      l=l-sum(q$homewin[ind])
      w=w-(length(q$homewin[ind])-sum(q$homewin[ind]))
      
      for (k in totalgames$opponents[[jind]]){
        kind=which(k==totalgames$awayteam)
        w2 = w2 + totalgames$wins[kind]
        l2 = l2 + totalgames$losses[kind]
        
        ind=which(totalgames$awayteam[i]==q$awayteam&k==q$hometeam)
        if(length(ind)>0){
          w2=w2-sum(q$homewin[ind])
          l2=l2-(length(q$homewin[ind])-sum(q$homewin[ind]))
        }
        
        ind=which(totalgames$awayteam[i]==q$hometeam&k==q$awayteam)
        if(length(ind)>0){
          l2=l2-sum(q$homewin[ind])
          w2=w2-(length(q$homewin[ind])-sum(q$homewin[ind]))
        }
      }
    }
    
    
    r[i] = 0.25*wp[i] + 0.5*w/(w+l) + 0.25*w2/(w2+l2)
  }
  
  l=0
  for(i in 1:nrow(q)){
    ra=r[q$awayteam[i]]
    rh=r[q$hometeam[i]]
    l = l+lambda*as.numeric(q$homewin[i])*rh+lambda*(1-as.numeric(q$homewin[i]))*ra-log(exp(lambda*rh)+exp(lambda*ra))
  }
  return(l)
  
}

rpimin<-function(b,q,totalgames){
  -rpi(b,q,totalgames)
}

(result <- nlm(rpimin,1,hessian=TRUE,q=games,totalgames=totalgames))



rpimse = function(param,q,totalgames,qtest){
  lambda=param[1]
  
  wp = totalgames$wins/(totalgames$wins+totalgames$losses)
  r = numeric(length(wp))
  names(r)=totalgames$awayteam
  for(i in 1:length(totalgames$awayteam)){
    
    w=0
    l=0
    w2=0
    l2=0
    for(j in totalgames$opponents[[i]]){
      jind=which(j==totalgames$awayteam)
      w = w + totalgames$wins[jind]
      l = l + totalgames$losses[jind]
      
      ind=which(totalgames$awayteam[i]==q$awayteam&j==q$hometeam)
      w=w-sum(q$homewin[ind])
      l=l-(length(q$homewin[ind])-sum(q$homewin[ind]))
      
      ind=which(totalgames$awayteam[i]==q$hometeam&j==q$awayteam)
      l=l-sum(q$homewin[ind])
      w=w-(length(q$homewin[ind])-sum(q$homewin[ind]))
      
      for (k in totalgames$opponents[[jind]]){
        kind=which(k==totalgames$awayteam)
        w2 = w2 + totalgames$wins[kind]
        l2 = l2 + totalgames$losses[kind]
        
        ind=which(totalgames$awayteam[i]==q$awayteam&k==q$hometeam)
        if(length(ind)>0){
          w2=w2-sum(q$homewin[ind])
          l2=l2-(length(q$homewin[ind])-sum(q$homewin[ind]))
        }
        
        ind=which(totalgames$awayteam[i]==q$hometeam&k==q$awayteam)
        if(length(ind)>0){
          l2=l2-sum(q$homewin[ind])
          w2=w2-(length(q$homewin[ind])-sum(q$homewin[ind]))
        }
      }
    }
    
    
    r[i] = 0.25*wp[i] + 0.5*w/(w+l) + 0.25*w2/(w2+l2)
  }
  
  mse=0
  for(i in 1:nrow(qtest)){
    ra=r[qtest$awayteam[i]]
    rh=r[qtest$hometeam[i]]
    p=exp(rh*lambda)/(exp(lambda*rh)+exp(lambda*ra))
    mse=mse+(as.numeric(qtest$homewin[i])-p)^2
  }
  names(mse) = NULL
  return(mse/nrow(qtest))
  
}
#set.seed(20)
#folds=sample(1:5,nrow(q),replace=T)
testmse3=0
for(i in 1:5){
  f=which(folds==i)
  tgames=gendata(games[-f,])
  result=nlm(rpimin,1,hessian=TRUE,q=games[-f,],totalgames=tgames)
  testmse3=testmse3+rpimse(result$estimate,games[-f,],tgames,games[f,])
  
}
testmse3=testmse3/5

rpi2 = function(param,q,totalgames){
  lambda=param[1]
  theta=param[2]
  gamma=param[3]
  
  wp = totalgames$wins/(totalgames$wins+totalgames$losses)
  r = numeric(length(wp))
  names(r)=totalgames$awayteam
  for(i in 1:length(totalgames$awayteam)){
    
    w=0
    l=0
    w2=0
    l2=0
    for(j in totalgames$opponents[[i]]){
      jind=which(j==totalgames$awayteam)
      w = w + totalgames$wins[jind]
      l = l + totalgames$losses[jind]
      
      ind=which(totalgames$awayteam[i]==q$awayteam&j==q$hometeam)
      w=w-sum(q$homewin[ind])
      l=l-(length(q$homewin[ind])-sum(q$homewin[ind]))
      
      ind=which(totalgames$awayteam[i]==q$hometeam&j==q$awayteam)
      l=l-sum(q$homewin[ind])
      w=w-(length(q$homewin[ind])-sum(q$homewin[ind]))
      
      for (k in totalgames$opponents[[jind]]){
        kind=which(k==totalgames$awayteam)
        w2 = w2 + totalgames$wins[kind]
        l2 = l2 + totalgames$losses[kind]
        
        ind=which(totalgames$awayteam[i]==q$awayteam&k==q$hometeam)
        if(length(ind)>0){
          w2=w2-sum(q$homewin[ind])
          l2=l2-(length(q$homewin[ind])-sum(q$homewin[ind]))
        }
        
        ind=which(totalgames$awayteam[i]==q$hometeam&k==q$awayteam)
        if(length(ind)>0){
          l2=l2-sum(q$homewin[ind])
          w2=w2-(length(q$homewin[ind])-sum(q$homewin[ind]))
        }
      }
    }
    
    
    r[i] = theta*wp[i] + (1-theta)*(gamma*w/(w+l) + (1-gamma)*w2/(w2+l2))
  }
  
  l=0
  for(i in 1:nrow(q)){
    ra=r[q$awayteam[i]]
    rh=r[q$hometeam[i]]
    l = l+lambda*as.numeric(q$homewin[i])*rh+lambda*(1-as.numeric(q$homewin[i]))*ra-log(exp(lambda*rh)+exp(lambda*ra))
  }
  return(l)
  
}

rpi2min<-function(b,q,totalgames){
  -rpi2(b,q,totalgames)
}


(result <- nlminb(c(1,0.5,0.5),rpi2min,hessian=TRUE,upper=c(Inf,1,1),lower=c(0,0,0),q=games,totalgames=totalgames))

rpimse2 = function(param,q,totalgames,qtest){
  lambda=param[1]
  theta=param[2]
  gamma=param[3]
  
  wp = totalgames$wins/(totalgames$wins+totalgames$losses)
  r = numeric(length(wp))
  names(r)=totalgames$awayteam
  for(i in 1:length(totalgames$awayteam)){
    
    w=0
    l=0
    w2=0
    l2=0
    for(j in totalgames$opponents[[i]]){
      jind=which(j==totalgames$awayteam)
      w = w + totalgames$wins[jind]
      l = l + totalgames$losses[jind]
      
      ind=which(totalgames$awayteam[i]==q$awayteam&j==q$hometeam)
      w=w-sum(q$homewin[ind])
      l=l-(length(q$homewin[ind])-sum(q$homewin[ind]))
      
      ind=which(totalgames$awayteam[i]==q$hometeam&j==q$awayteam)
      l=l-sum(q$homewin[ind])
      w=w-(length(q$homewin[ind])-sum(q$homewin[ind]))
      
      for (k in totalgames$opponents[[jind]]){
        kind=which(k==totalgames$awayteam)
        w2 = w2 + totalgames$wins[kind]
        l2 = l2 + totalgames$losses[kind]
        
        ind=which(totalgames$awayteam[i]==q$awayteam&k==q$hometeam)
        if(length(ind)>0){
          w2=w2-sum(q$homewin[ind])
          l2=l2-(length(q$homewin[ind])-sum(q$homewin[ind]))
        }
        
        ind=which(totalgames$awayteam[i]==q$hometeam&k==q$awayteam)
        if(length(ind)>0){
          l2=l2-sum(q$homewin[ind])
          w2=w2-(length(q$homewin[ind])-sum(q$homewin[ind]))
        }
      }
    }
    
    
    r[i] = theta*wp[i] + (1-theta)*(gamma*w/(w+l) + (1-gamma)*w2/(w2+l2))
  }
  
  mse=0
  for(i in 1:nrow(qtest)){
    ra=r[qtest$awayteam[i]]
    rh=r[qtest$hometeam[i]]
    p=exp(rh*lambda)/(exp(lambda*rh)+exp(lambda*ra))
    mse=mse+(as.numeric(qtest$homewin[i])-p)^2
  }
  names(mse) = NULL
  return(mse/nrow(qtest))
  
}
#set.seed(20)
#folds=sample(1:5,nrow(q),replace=T)
testmse4=0
for(i in 1:5){
  f=which(folds==i)
  tgames=gendata(games[-f,])
  result=nlminb(c(1,0.5,0.5),rpi2min,hessian=TRUE,upper=c(Inf,1,1),lower=c(0,0,0),q=games[-f,],totalgames=tgames)
  testmse4=testmse4+rpimse2(result$par,games[-f,],tgames,games[f,])
  
}
testmse4=testmse4/5
