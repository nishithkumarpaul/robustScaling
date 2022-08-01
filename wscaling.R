#######################################
#####                              ####
##### Function of Weighted Scaling ####
#####                              ####
#######################################

#####################################################
##### Programmer:                                ####
##### Nishith Kumar and Biplab Biswas            ####
##### Dept. of Statistics, BSMRSTU               ####
#####                                            ####
#####################################################

wScaling<-function(x){
  if (sum(is.na(x))>0){
  print("Data contains missing Values")
  }
  nRow<-dim(x)[1]
  nCol<-dim(x)[2]
  m<-matrix(rep(1,nRow*nCol),ncol=nCol)
  w<-matrix(rep(NA,nRow*nCol),ncol=nCol)
  rmed<-apply(x,1,median, na.rm=TRUE)
  rmad<- apply(x,1,mad, na.rm=TRUE)
    for (i in 1:nRow){
      for (j in 1:nCol)
        {
          w[i,j]<-min(m[i,j],((1.96^2)/((x[i,j]-rmed[i])/rmad[i])^2))
        }
     }

###################################################
####                                           ####
####    Weighted Mean Calculation              ####
####                                           ####
###################################################

wx<-w*x
sum_wx<-as.matrix(apply(wx,1,sum))
sum_w<-as.matrix(apply(w,1,sum))
x_bar<-sum_wx/sum_w

###################################################
####                                           ####
####    Weighted Standard Deviation            ####
####                                           ####
###################################################

sum_devM<-matrix(rep(NA,nRow*nCol),ncol=nCol)

for (i in 1:nRow){
  for (j in 1:nCol){
  sum_devM[i,j]<-((w[i,j]*x[i,j])-x_bar[i])^2
  }
}
sumsq<-as.matrix(apply(sum_devM,1,sum))
s_w1<-sumsq/sum_w
s_w<-sqrt(s_w1)

###################################################
####                                           ####
####    Weighted Scaling Matrix Calculation    ####
####                                           ####
###################################################

wsp<-matrix(rep(NA,nRow*nCol),ncol=nCol)
for (i in 1:nRow){
  for (j in 1:nCol){
  if (w[i,j]==1)wsp[i,j]<-(x[i,j]-x_bar[i])/s_w[i]
  else if(w[i,j]!=1)wsp[i,j]<-(((w[i,j]*x[i,j])-x_bar[i])/s_w[i]) 
  }
}

dataWS<-wsp
return(dataWS)
}

###################################################
####                                           ####
####    Artificial Data                        ####
####                                           ####
#### Row indicates: Different Metabolite       ####
#### Column indicates: Different Subjects/sample ##
###################################################

set.seed(20)
x11<-matrix(rnbinom(106*118,200,0.60),ncol=106)
x12<-matrix(rnbinom(91*118,75,0.55),ncol=91)
x22<-matrix(rnbinom(197*118,50,0.40),ncol=197)
x1c<-as.matrix(cbind(x11,x12))
dummyData<-rbind(x1c,x22)
nRow<-dim(dummyData)[1]
nCol<-dim(dummyData)[2]
dataM<-as.matrix(dummyData)

###################################################
####                                           ####
####    Applying Weighted Scaling in           ####
####    Artificial Data                        ####
####                                           ####
###################################################

wData<-wScaling(dataM)

###################################################
####  The End                                  ####
###################################################
