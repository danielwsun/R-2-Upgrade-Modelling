#In this script, the model is with two absorbing states#
########################################################

## Global parameters
R<-20# the minimum number of servers requested, below this the service will crash or face a serious problem
MaxFailure<-6 #Maximum number of failures
NumVMs<-30 # number of all instances
Granu<-2#Granularity of rolling upgrade
Pops<-0.9#probability or ratio of sucessful operations

###########  system failures in the above #################
#Function to create the system failure probability.
#Non-zero failure probabilities generated from exponential distribution
#Pf is the function to know the probability of appearing a certain number of failures

FailureExpNum<-c(1:MaxFailure)
FailureExpDis<-c(0.3, 0.1, 0.07, 0.06, 0.04, 0.03)
Pf<-function(k){
  if(k<0) return(NaN)
  if(k == 0 || k > MaxFailure) 
    return(1-sum(FailureExpDis,FALSE))
  else  
    return(FailureExpDis[k])
}
###########  system failures in the above #################

############  operation failures    #######################
#Define operation failures
PopsDis<-function(X,x){
  return(dbinom(x,size=X,prob=Pops))
}
############  operation failures    #######################

#############System failure configuration distribution#####
#System failure configuration
Pconf<-function(w,u,z,k,i,g){
  numerator<-choose(g,w)*choose(i,u)*choose(NumVMs-g-i,z)
  denominator<-choose(NumVMs, k)
  temp<-numerator/denominator
  return(temp)
}
#############System failure configuration distribution#####

################Transition Probabilities ##################
#transition probability
#Forward
PForward<-function(i,j){
  temp<-0
  G<-Granu
  for(k in 0:MaxFailure){
      if(min(G-j,k)<0) next
      for(w in 0:min(G-j,k)){
        if(min(G-j-w,k-w,i)<0) next
        for(u in 0:min(G-j-w,k-w,i)){
          temp<-temp+PopsDis(G-w,u+j)*Pf(k)*Pconf(w,u,k-w-u,k,i,G)
        }
    }
    return(temp)
  }
}

#backward
PBackward<-function(i,j){
  temp<-0
  G<-Granu
  for(k in j:MaxFailure){
      if(min(G,k-j)<0) next
      for(w in 0:min(G,k-j)){
        if(min(G+j-w,k-w,i)<0) next
        for(u in j:min(G+j-w,k-w,i)){
          temp<-temp+PopsDis(G-w,u-j)*Pf(k)*Pconf(w,u,k-w-u,k,i,G)
        }
    }
    return(temp)
  } 
}

#staying
PStay<-function(i){
  temp<-0
  G<-Granu
  for(k in 0:MaxFailure){
      for(w in 0:min(G,k)){
        for(u in 0:min(G-w,k-w,i)){
          temp<-temp+PopsDis(G-w,u)*Pf(k)*Pconf(w,u,k-w-u,k,i,G)
        }
    }
    return(temp)
  }
}
################Transition Probabilities############################

####################Generate P matrix for no death chain#######################
##Constructing matrices
#Remember in R, index to a matrix starts from 1 but not 0
#Construct P matrix
PMatrix<-function(){
  TempMatrix <- matrix(nrow=NumVMs+1, ncol=NumVMs+1)
  for(i in 1:(NumVMs+1)){
    for(j in 1:(NumVMs+1)){
      if(i==NumVMs+1){
        TempMatrix[i, j]<-0
        next
      }
      if(i-MaxFailure > j){
        TempMatrix[i,j]<-0
        next
      }else{
        if(j > i + Granu){
          TempMatrix[i,j]<-0
          next
        }else{
          J <- abs(i-j) #In the model, j is the relative offset, but here j is absolute position in matrices
          I <- i-1
          if(j<i)TempMatrix[i,j]<-PBackward(I,J)
          if(j==i)TempMatrix[i,j]<-PStay(I)
          if(j>i)TempMatrix[i,j]<-PForward(I,J)
        }
      }
    }
  }
  TempMatrix[NumVMs+1, NumVMs+1]<-1
  return(TempMatrix)  
}

P<-PMatrix()
####################Generate P matrix################################



####################Generate PP and QQ matrices##############################
PP<-P[(R+1):(NumVMs+1),(R+1):(NumVMs+1)]
ToDeath<-rep(0,(NumVMs-R+1))
for(i in (R+1):(R+MaxFailure)){
  for(j in (R-MaxFailure+1):R){
    ToDeath[i-R]<-ToDeath[i-R]+P[i,j]
  }
}
PP<-cbind(PP,ToDeath)
InDeath<-rep(0,(NumVMs-R+2))
InDeath[NumVMs-R+2]<-1
PP<-rbind(PP,InDeath)
QQ<-PP[1:(NumVMs-R),1:(NumVMs-R)]
RR<-PP[1:(NumVMs-R),(NumVMs-R+1):(NumVMs-R+2)]
####################Generate PP and QQ matrices##############################



####################Calculate fundamental matrix#############################
IoGen<-diag(NumVMs-R)#the indendity matrix for general model
FMatrix<-solve(IoGen-QQ)#inverse a matrix
####################Calculate fundamental matrix#############################

###################to be dead or alive##############
#Expected steps to the final from any state
EStep<-rowSums(FMatrix)
#The variance
VStep<-(2*FMatrix-IoGen)%*%EStep - EStep^2
#eg. EStep[1] is the expected number of steps to the final from the state 0
B<-FMatrix%*%RR #two columns in B, the left column is the probabilities to finish successfully and the right one is the probabilities to the death
###################to be dead or alivePP##############

