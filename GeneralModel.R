## Global parameters
MaxFailure<-6 #Maximum number of faults
NumVMs<-30 # number of all instances
Granu<-2#Granularity of rolling upgrade
Pops<-0.9#probability or ratio of sucessful operations

###########  system failures in the above #################
#Function to create the system failure probability.
#Non-zero failure probabilities generated from exponential distribution
#Pf is the function to know the probability of appearing a certain number of failures

FailureExpNum<-c(1:MaxFailure)
FailureExpDis<-c(0.3, 0.1, 0.07, 0.06, 0.04, 0.03)#Set probabiliteis for different numbers of faults
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

####################Generate P and Q matrices#######################
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

QMatrix<-function(){
  TempMatrix <- matrix(nrow=NumVMs, ncol=NumVMs)
  for(i in 1:NumVMs){
    for(j in 1:NumVMs){      
      if(i-MaxFailure > j){
        TempMatrix[i,j]<-0
        next
      }else{
        if(j > i + Granu){
          TempMatrix[i,j]<-0
          next
        }else{
          J <- abs(i-j) #In the model, j is the relative offset, but here j is absolute position in matrices
          I <- i -1
          if(j<i)TempMatrix[i,j]<-PBackward(I,J)
          if(j==i)TempMatrix[i,j]<-PStay(I)
          if(j>i)TempMatrix[i,j]<-PForward(I,J)
        }
      }
    }
  }
  return(TempMatrix)  
}
P<-PMatrix()
Q<-QMatrix()
####################Generate P and Q matrices################################


####################Calculate fundamental matrix#############################
IoGen<-diag(NumVMs)#the indendity matrix for general model
FMatrix<-solve(IoGen-Q)#inverse a matrix
####################Calculate fundamental matrix#############################

###################some simple numeric results from genral model##############
#Expected steps to the final from any state
EStep<-rowSums(FMatrix)
#The variance
VStep<-(2*FMatrix-IoGen)%*%EStep - EStep^2
#eg. EStep[1] is the expected number of steps to the final from the state 0
###################some simple numeric results from genral model##############

####################Data Collection##########################################
#do your own data collection code here
