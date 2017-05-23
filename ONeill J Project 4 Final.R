# ################################
# Jack ONeill
# ECON 524
# Project 4
# 12/1/16
# Last Edit: 12/12/2016
# Due Wednesday 12/14/2016
# ################################
rm(list=ls())
###################################

data = read.csv("pddata.csv")
df = data
sv = c(0,0) # start values for parameters

################################################################
# log likelihood function for entire data frame

ll_rl = function(parms, data){
  h = parms[1]
  d = parms[2]
  ppa = as.numeric(2000)
  ppb = as.numeric(2000)
  ll = as.numeric(2000)
  for(i in (0:19)*100){
    for(j in 1:100){
      if(j == 1){
        ppa[j+i] = 0
        ppb[j+i] = 0
        }else{
        ppa[j+i] = d * ppa[j+i-1] + ifelse(data[j+i-1,1]==0, data[j+i-1,4],0)
        ppb[j+i] = d * ppb[j+i-1] + ifelse(data[j+i-1,1]==1, data[j+i-1,4],0)
        }
      if(data[j+i, 1]==0){
        ll[j+i] = h*ppa[j+i] - log(exp(h*ppa[j+i]) + exp(h*ppb[j+i]))
        }else{
        ll[j+i] = h*ppb[j+i] - log(exp(h*ppa[j+i]) + exp(h*ppb[j+i]))
      }
    }
  }
nloglike = -sum(ll)
return(nloglike)
}

ll_r2 = function(parms, data){
  h = parms[1]
  d = parms[2]
  ppa = as.numeric(2000)
  ppb = as.numeric(2000)
  ll = as.numeric(2000)
  for(i in (0:19)*100){
    for(j in 1:100){
      if(j == 1){
        ppa[j+i] = 0
        ppb[j+i] = 0
      }else{
        ppa[j+i] = d * ppa[j+i-1] + ifelse(data[j+i-1,1]==0, data[j+i-1,9],0)
        ppb[j+i] = d * ppb[j+i-1] + ifelse(data[j+i-1,1]==1, data[j+i-1,9],0)
      }
      if(data[j+i,1]==0){
        ll[j+i] = h*ppa[j+i] - log(exp(h*ppa[j+i]) + exp(h*ppb[j+i]))
      }else{
        ll[j+i] = h*ppb[j+i] - log(exp(h*ppa[j+i]) + exp(h*ppb[j+i]))
      }
    }
  }
  nloglike = -sum(ll)
  return(nloglike)
}

ste = function(x){
  se = sd(x)/sqrt(length(x))
  return(se)
}

################################################################
# Estimation 1
################################################################

maxll = nlminb(sv, ll_rl, lower = c(0,0), upper = c(Inf,1), data = data)
result = maxll$par
result
llresult = maxll$objective
llresult

# Bootstrap 1

set.seed(100)
bmatrix = matrix(0, nrow=10, ncol=2)

for(i in 1:10){
  bdata = data.frame()
  ind = sample(1:20, replace = T)
  for(j in ind){
    bdata = rbind(bdata,data[data$id==j,])
  }
  est = nlminb(sv, ll_rl, lower=c(0,0), upper = c(Inf,1), data = bdata ) 
  bmatrix[i,] = est$par
}

serror = apply(bmatrix, 2, sd)

serror


###################################################################
# Estimation 2
###################################################################

maxll2 = nlminb(sv, ll_r2, lower=c(0,0), upper = c(Inf,1), data = data)
result2 = maxll2$par
result2
llresult2 = maxll2$objective
llresult2

set.seed(100)

bmatrix = matrix(0, nrow=10, ncol=2)

for(i in 1:10){
  bdata = data.frame()
  ind = sample(1:20, replace = T)
  for(j in ind){
    bdata = rbind(bdata,data[data$id==j,])
  }
  est = nlminb(sv, ll_r2, lower=c(0,0), upper = c(Inf,1), data=bdata ) 
  bmatrix[i,] = est$par
}



serror2 = apply(bmatrix, 2, sd)

serror2

###########################################################################################