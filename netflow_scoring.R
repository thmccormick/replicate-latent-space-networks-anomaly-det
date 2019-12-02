#Setup
######

# help function
vmv <- function(m, v) {
  return(m[1,1]*v[1]^2 + (m[2,1] + m[1,2])*v[1]*v[2] + m[2,2]*v[2]^2)
}

# auxillary data
ids <- read.csv("data/new_id_dictionary.csv")

n = max(ids$new_ids)

load("data/period_weights.RData")



# Score network activity
########################

results_logit <- list()
results_logit_base <- list()

start = Sys.time()

for(day in 2:90) {
  print(paste(day, difftime(Sys.time(),start,units="mins")))
  data <- read.csv(paste("data/netflow_re4h-d",formatC(day,width=2,flag="0"),".csv",sep=""))
  if(day == 2) {
    periods = 9:12-1
  } else {
    periods = 1:7 + day*6 - 7 -1
  }
  
  score_logit = rep(0,dim(data)[1])
  score_logit_base = rep(0,dim(data)[1])
  
  for(period in periods) {
    load(paste("data/results_logit_p",period,".RData",sep=""))
    results_logit[[period]] <- results
    load(paste("data/results_basic_logit_p",period,".RData",sep=""))
    results_logit_base[[period]] <- results
  }
  
  for(row in 1:dim(data)[1]) {
    i = data[row,1]
    j = data[row,2]
    period = data[row,3]-1
    
    # Account for case-control
    logit_adj = period_weights$mu_mod[period_weights$tod == ((period+1) %% 6) & period_weights$dow == (((period+1) %/% 6) %% 7)]
    
    if(period != 7) {
      #Score for logit model
      f = results_logit[[period+1]]$f
      u_mean = results_logit[[period]]$u_mean[i,]
      v_mean = results_logit[[period]]$v_mean[j,]
      u_var = results_logit[[period]]$u_var[,,i]
      v_var = results_logit[[period]]$v_var[,,j]
      mean_f = results_logit[[period]]$mu_mean + results_logit[[period]]$alpha_mean[i] + results_logit[[period]]$beta_mean[j] + sum(u_mean * v_mean) + logit_adj
      sigma_f = results_logit[[period]]$mu_var*f[1] + results_logit[[period]]$alpha_var[i]*f[2] + results_logit[[period]]$beta_var[j]*f[2] + 
        f[3]^2*sum((u_var) * (v_var)) + f[3]*vmv(v_var, u_mean) + f[3]*vmv(u_var, v_mean)
      score_logit[row] = (1+pi/8*sigma_f)^(-0.5)*mean_f
      
      #Score for base logit model
      f = results_logit_base[[period+1]]$f
      mean_f = results_logit_base[[period]]$mu_mean + results_logit_base[[period]]$alpha_mean[i] + results_logit_base[[period]]$beta_mean[j] + logit_adj
      sigma_f = results_logit_base[[period]]$mu_var*f[1] + results_logit_base[[period]]$alpha_var[i]*f[2] + results_logit_base[[period]]$beta_var[j]*f[2]
      score_logit_base[row] = (1+pi/8*sigma_f)^(-0.5)*mean_f
    }
  }
  
  for(period in periods) {
    results_linear[[period]] <- NULL
    results_linear_base[[period]] <- NULL
    results_logit[[period]] <- NULL
    results_logit_base[[period]] <- NULL
  }
  data$score_logit = score_logit
  data$score_logit_base = score_logit_base
  
  save(data, file = paste("data/scores_d",day,".RData",sep=""))
}
