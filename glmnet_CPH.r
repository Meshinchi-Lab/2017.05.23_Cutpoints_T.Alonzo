#Jenny Smith


#March 13, 2017 


#purpose: Use cox proportional hazards generalized linear model to identify EFS due to miRNA expression


# install.packages("glmnet")

library("glmnet")
library("survival")

# setwd("H:/miRNAseq_Analysis/")
# setwd("/Volumes/jlsmith3/miRNAseq_Analysis/")
setwd(file.path(TARGET,"RNA/miRNAseq/2017.05.23_Cutpoints_T.Alonzo"))


#Use this for the gaussian vignettes
# data(QuickStartExample)
 
#use this for the Cox models. 
#https://web.stanford.edu/~hastie/glmnet/glmnet_alpha.html#cox
data(CoxExample)

#another dataset for the cox model examples
#https://cran.r-project.org/web/packages/glmnet/vignettes/Coxnet.pdf
# data(VignetteExample)
# 
# load("VignetteExample.rdata") #cant find file 



#IMPORTANT NOTE: input expression values should be transformed 
#Emilia Lim used log2 RPM (was this log2(x+1)?)
#others have used the VST method of transformation of the raw counts


#example of y - which is the survival data
y[1:5,]

#1000 obs, 2 variables (time and status)
dim(y)

#sort into patietns are in same order  
#example of x - this is the input matrix of expression values (log2 RPM for Emilia)
x[1:5,]


#1000 obs, and 30 samples
dim(x)


#fit the expression data with the survival data. defaults used
fit <- glmnet(x,y, family="cox")

#plot the coeff
par(pty="s")
plot(fit, label = TRUE)


#summary of the glment path at each step
# It shows from left to right the number of nonzero coefficients (Df), the percent (of null) deviance explained (%dev) 
# and the value of λλ (Lambda). 
print(fit)


#column names are s0 through s48
#the variables (30 of them) become obs(not correct terminology)
#this is the multivariate analysis
colnames(coef(fit))


#extract the coefficients. WHY ARE WE USING S=0.05????
#answer: obtain the actual coefficients at one or more s's within the range of the sequence:
#second question: uhhh what range of sequence?
#second answer:
# (why s and not lambda? In case later we want to allow one to specify the model size in other ways.) 
#
coef(fit, s=0.001)

fit$beta

#extract the coefficients for each variable 1-10
coef(fit, s=lambda.min)

#CV fit
cvfit = cv.glmnet(x, y, family = "cox")

plot(cvfit)


coef.min = coef(cvfit, s = "lambda.min")

active.min = which(coef.min != 0)
index.min = coef.min[active.min]















