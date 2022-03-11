#Jenny Smith 

#Original Author:Emilia Lim

#Purpose: Use xtile method to seperate clinical groups *rather than more arbitray median/quantile cut offs*

setwd(file.path(TARGET,"RNA/miRNAseq/2017.05.23_Cutpoints_T.Alonzo"))

findExpGroupsDetail = function(names,exp_vals,sdfp,survival_type){
  require('dplyr')
  #names = patient names
  #exp_vals = expression levels of that gene for each patient listed in names
  #sdfp = survival data (Event and Time to Event) for each patient listed in names
  #survival_type = "OS" or "EFS"
  names(exp_vals) = as.character(names)
  exp_vals_sorted = sort(exp_vals)
  sdfp = sdfp[names(exp_vals_sorted),]
  best_p_val = 1
  best_i = 1
  p_values = c()
  for(i in floor(length(exp_vals_sorted)*0.1):(ceiling(length(exp_vals_sorted)*0.9))){
    c_clusters = c(rep("A",i),rep("B",(length(exp_vals_sorted)-i)))
    if(survival_type=="OS"){
      sd = survdiff(Surv(OS_time_years, OS_event_ID==1 ) ~ c_clusters, data=sdfp)
      coxvar = try(coxph(formula = Surv(OS_time_years, OS_event_ID==1 ) ~ c_clusters, data = sdfp))
    }
    if(survival_type=="EFS"){
      sd = survdiff(Surv(EFS_time_years, EFS_event_type_ID==1 ) ~ c_clusters, data=sdfp)
      coxvar = try(coxph(formula = Surv(EFS_time_years, EFS_event_type_ID==1 ) ~ c_clusters, data = sdfp))
    }
    p_values = c(p_values,summary(coxvar)$coefficient[5])
    p_val = 1 - pchisq(sd$chisq, 1)
    if(as.numeric(p_val)<as.numeric(best_p_val)){
      best_i = i
      best_p_val = p_val
    }
  }
  
  GROUPA = names(exp_vals_sorted[1:best_i])
  GROUPB = names(exp_vals_sorted[(best_i+1):length(exp_vals_sorted)])
  groupings = ifelse(names%in%GROUPA,"Low","High")
  cut_point = mean(exp_vals_sorted[best_i:(best_i+1)])
  r_data = list(groupings,cut_point)
  names(r_data) = c("Groupings","Cut_Point")
  return(r_data)
}






