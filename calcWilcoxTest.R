#Jenny Smith


#original Author: Emilia Lim 

#Puprpose: DE analysis of miRNA

setwd(file.path(TARGET,"RNA/miRNAseq/2017.05.23_Cutpoints_T.Alonzo"))

calcWilcoxTest = function(df,libsA,libsB,log=T,paired=F,aname="A",bname="B",dcores=F){
  #df = your df of expression values
  #libsA/libsB = vector of strings of colnames in either group A or group B
  #dcores = number of cores to use for analysis. "F" will use only one. 
  libsA = intersect(colnames(df),libsA)
  libsB = intersect(colnames(df),libsB)
  df = df[,c(libsA,libsB)]
  if(log){
    df[df<1] = 1
    df = log2(df)
  }
  libAidx = c(1:length(libsA))
  libBidx = c((length(libsA)+1):(dim(df)[2]))
  
  if(dcores){
    library(doMC)
    library(foreach)
    registerDoMC(dcores)
    wt_res = foreach(i = 1:nrow(df))%dopar%{
      values = as.numeric(df[i,])
      x = values[libAidx]
      y = values[libBidx]
      wt = wilcox.test(x,y,paired=paired)
      if(log){
        fc = mean(x)-mean(y)
        c(mean(values),mean(x),mean(y),fc,wt$p.value)
      }else{
        fc = log2(mean(x))-log2(mean(y))
        c(log2(mean(values)),log2(mean(x)),log2(mean(y)),fc,wt$p.value)
      }
    }
    wt_res = do.call(rbind,wt_res)
    rownames(wt_res) = rownames(df)
  }else{
    wt_res = apply(df,1,function(values){
      x = values[libAidx]
      y = values[libBidx]
      wt = wilcox.test(x,y,paired=paired)
      if(log){
        fc = mean(x)-mean(y)
        c(mean(values),mean(x),mean(y),fc,wt$p.value)
      }else{
        fc = log2(mean(x))-log2(mean(y))
        c(log2(mean(values)),log2(mean(x)),log2(mean(y)),fc,wt$p.value)
      }
    })
    wt_res = t(wt_res)
  }
  colnames(wt_res) = c("log2_base_mean",paste("log2_mean_",aname,"_n",length(libsA),sep=""), 
                       paste("log2_mean_",bname,"_n",length(libsB),sep=""), 
                       paste("log2_fold_change_",aname,"..",bname,sep=""),"p_val")
  wt_res = data.frame(wt_res)
  adj_p_val = p.adjust(wt_res$p_val,method="BH")
  wt_res = cbind(wt_res,adj_p_val)
}


#wilcox test in R . This used a un-paired, two-sample, wilcox test, with two sided alternative hypothesis.
#this is a wilcoxon rank-sum test 
#the null hypothesis is that the distributions of x and y differ by a location shift of mu
# and the alternative is that they differ by some other location shift
#the null hypothesis that data in x and y are samples from continuous distributions with equal medians
#This takes the measurements for each group (A,B) and ranks them when combined in a ordered numberline
#then it sums the ranks PER group. So if group A had all larger expression values, it would have
#larger ranks and would have sum of ranks > group B. 
#in two sided test, this would reject the null and A is greater than B.
