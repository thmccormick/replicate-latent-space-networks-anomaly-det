# Setup
#######

# helper functions
reltol <- function(x,y) {
  sum(abs(x-y))/sum(abs(x))
}

expit <- function(x) {y = exp(x); y/(1+y)}

weighted_sum <- function(x1, x2, w) 1/(1+w)*x1 + w/(1+w)*x2



# Auxillary data
ids <- read.csv("data/new_id_dictionary.csv")

n = max(ids$new_ids)

load("data/period_weights.RData")





# Estimation
############

sigma_prior = c(0.5,0.25,0.25,0)

damping = 2
epsilon = 3e-04
ldim = 2
average_zeros = 500000 #on average, sample 500000 zeros in each period


start = Sys.time()

for(day in 2:90) {
  data <- read.csv(paste("data/netflow_re4h-d",formatC(day,width=2,flag="0"),".csv",sep=""))
  
  # First day is shortened
  if(day == 2) {
    periods = 8:11
  } else {
    periods = 1:6 + day*6 - 7 
  }
  
  for(period in periods) {
    period_weight = 1/period_weights$weight[period_weights$tod == (period %% 6) & period_weights$dow == ((period %/% 6) %% 7)]
    
    # Select sample of data for inference
    rows = as.matrix(data[data$TimeFrame == period,1:2])
    num_ones = dim(rows)[1]
    num_zeros = round(average_zeros*period_weight)
    temp.vals <- sample.int(n^2-1, num_zeros)
    sampled.vals = cbind(temp.vals %/% n + 1, temp.vals %% n + 1)
    rows = rbind(rows, sampled.vals)
    rows = rows[!duplicated(rows),]
    num_zeros = dim(rows)[1] - num_ones
    rows = cbind(rows, c(rep(1,num_ones), rep(0,num_zeros)))
    rows = rows[sample(dim(rows)[1]),]

    # initiate priors for period
    if(period == 8) {
      q_mu_mean <- -2
      q_mu_var <- sigma_prior[1]
      q_alpha_mean <- rep(0, n)
      q_alpha_var <- rep(sigma_prior[2], n)
      q_beta_mean <- rep(0, n)
      q_beta_var <- rep(sigma_prior[2], n)

      f = 0
    } else {
      forgetting <- as.matrix(expand.grid(0:4, 0:4))
      colnames(forgetting) = NULL
      forgetting[forgetting == 0] = 1.1
      forgetting[forgetting == 1] = 1
      forgetting[forgetting == 2] = 1.01
      forgetting[forgetting == 3] = 2
      forgetting[forgetting == 4] = 1.001
      
      # Estimate amount of forgetting for period
      ll = rep(0, dim(forgetting)[1])
      for(row in 1:n) {
        i = rows[row,1]
        j = rows[row,2]
        sign = 1-2*rows[row,3]
        
        mean_f = q_mu_mean + q_alpha_mean[i] + q_beta_mean[j]
        for(row2 in 1:dim(forgetting)[1]) {
          f = forgetting[row2,]
          sigma_f = q_mu_var*f[1] + q_alpha_var[i]*f[2] + q_beta_var[j]*f[2]
          ll[row2] = ll[row2] + log(expit(-sign*(1+pi/8*sigma_f)^(-0.5)*mean_f))
        }
      }
      
      f = forgetting[which.max(ll),]
      print(f)
      
      q_mu_var <- q_mu_var*f[1]
      q_alpha_var <- q_alpha_var*f[2]
      q_beta_var <- q_beta_var*f[2]
    }
    
    # Initiate messages
    msg_mu_mean <- rep(0, dim(rows)[1])
    msg_mu_var <- rep(Inf, dim(rows)[1])
    msg_alpha_mean <- rep(0, dim(rows)[1])
    msg_alpha_var <- rep(Inf, dim(rows)[1])
    msg_beta_mean <- rep(0, dim(rows)[1])
    msg_beta_var <- rep(Inf, dim(rows)[1])

    counter = 0
    
    while(TRUE) {
      q_mu_old_mean = q_mu_mean
      q_alpha_old_mean = q_alpha_mean
      q_beta_old_mean = q_beta_mean
      q_mu_old_var = q_mu_var
      q_alpha_old_var = q_alpha_var
      q_beta_old_var = q_beta_var
      
      for(row in 1:dim(rows)[1]) {
        i = rows[row,1]
        j = rows[row,2]
        sign = 1-2*rows[row,3]

        q_tmp_mu_mean = q_mu_mean
        q_tmp_mu_var = q_mu_var
        q_tmp_alpha_mean = q_alpha_mean[i]
        q_tmp_alpha_var = q_alpha_var[i]
        q_tmp_beta_mean = q_beta_mean[j]
        q_tmp_beta_var = q_beta_var[j]

        msg_tmp_mu_mean = msg_mu_mean[row]
        msg_tmp_mu_var = msg_mu_var[row]
        msg_tmp_alpha_mean = msg_alpha_mean[row]
        msg_tmp_alpha_var = msg_alpha_var[row]
        msg_tmp_beta_mean = msg_beta_mean[row]
        msg_tmp_beta_var = msg_beta_var[row]

        # Define "g functions" (see bottom of page 8)
        g_mu_var = (1/q_mu_var + 1/msg_mu_var[row])^(-1)
        g_mu_mean = g_mu_var*(q_mu_mean/q_mu_var + msg_mu_mean[row]/msg_mu_var[row])
        g_alpha_var = (1/q_alpha_var[i] + 1/msg_alpha_var[row])^(-1)
        g_alpha_mean = g_alpha_var*(q_alpha_mean[i]/q_alpha_var[i] + msg_alpha_mean[row]/msg_alpha_var[row])
        g_beta_var = (1/q_beta_var[j] + 1/msg_beta_var[row])^(-1)
        g_beta_mean = g_beta_var*(q_beta_mean[j]/q_beta_var[j] + msg_beta_mean[row]/msg_beta_var[row])

        tmp_weight = exp(sign*(g_mu_mean + g_alpha_mean + g_beta_mean) + g_mu_var/2 + g_alpha_var/2 + g_beta_var/2)
        
        # Equation 8 via equations 21-23
        tmp_mean = g_mu_mean
        tmp_mean2 = g_mu_mean+sign*g_mu_var
        q_prime_mu_mean = weighted_sum(tmp_mean,tmp_mean2,tmp_weight)
        q_prime_mu_var = weighted_sum(tmp_mean^2+g_mu_var,tmp_mean2^2+g_mu_var,tmp_weight) - q_prime_mu_mean^2
        
        tmp_mean = g_alpha_mean
        tmp_mean2 = tmp_mean+sign*g_alpha_var
        q_prime_alpha_mean = weighted_sum(tmp_mean,tmp_mean2,tmp_weight)
        q_prime_alpha_var = weighted_sum(tmp_mean^2+g_alpha_var,tmp_mean2^2+g_alpha_var,tmp_weight) - q_prime_alpha_mean^2
        
        tmp_mean = g_beta_mean
        tmp_mean2 = tmp_mean+sign*g_beta_var
        q_prime_beta_mean = weighted_sum(tmp_mean,tmp_mean2,tmp_weight)
        q_prime_beta_var = weighted_sum(tmp_mean^2+g_beta_var,tmp_mean2^2+g_beta_var,tmp_weight) - q_prime_beta_mean^2
        
        # Equation 9
        q_mu_var = (damping/q_tmp_mu_var + (1-damping)/q_prime_mu_var)^(-1)
        q_mu_mean = q_mu_var*(damping/q_tmp_mu_var*q_tmp_mu_mean + (1-damping)/q_prime_mu_var*q_prime_mu_mean)
        q_alpha_var[i] = (damping/q_tmp_alpha_var + (1-damping)/q_prime_alpha_var)^(-1)
        q_alpha_mean[i] = q_alpha_var[i]*(damping/q_tmp_alpha_var*q_tmp_alpha_mean + (1-damping)/q_prime_alpha_var*q_prime_alpha_mean)
        q_beta_var[j] = (damping/q_tmp_beta_var + (1-damping)/q_prime_beta_var)^(-1)
        q_beta_mean[j] = q_beta_var[j]*(damping/q_tmp_beta_var*q_tmp_beta_mean + (1-damping)/q_prime_beta_var*q_prime_beta_mean)
        
        # Equation 10
        msg_mu_var[row] = (1/q_mu_var + 1/msg_mu_var[row] - 1/q_tmp_mu_var)^(-1)
        if(msg_mu_var[row] != Inf) {
          msg_mu_mean[row] = msg_mu_var[row]*(q_mu_mean/q_mu_var + msg_mu_mean[row]/msg_tmp_mu_var - q_tmp_mu_mean/q_tmp_mu_var)  
        }
        
        msg_alpha_var[row] = (1/q_alpha_var[i] + 1/msg_alpha_var[row] - 1/q_tmp_alpha_var)^(-1)
        if(msg_alpha_var[row] != Inf) {
          msg_alpha_mean[row] = msg_alpha_var[row]*(q_alpha_mean[i]/q_alpha_var[i] + msg_alpha_mean[row]/msg_tmp_alpha_var - q_tmp_alpha_mean/q_tmp_alpha_var)  
        }
        msg_beta_var[row] = (1/q_beta_var[j] + 1/msg_beta_var[row] - 1/q_tmp_beta_var)^(-1)
        if(msg_beta_var[row] != Inf) {
          msg_beta_mean[row] = msg_beta_var[row]*(q_beta_mean[j]/q_beta_var[j] + msg_beta_mean[row]/msg_tmp_beta_var - q_tmp_beta_mean/q_tmp_beta_var)  
        }
      }
      
      criteria = reltol(c(q_mu_var,q_alpha_var,q_beta_var), c(q_mu_old_var,q_alpha_old_var,q_beta_old_var))
      print(paste("Period: ",period,", Iter: ",counter,", Time: ",format(round(difftime(Sys.time(),start,units="mins"),3),nsmall=4),", Rel Tol: ", format(criteria, nsmall=4),sep=""))
      
      if(criteria < epsilon) {
        break
      }
      counter = counter+1
    }
    
    results = list(mu_var = q_mu_var, 
                   alpha_var = q_alpha_var,
                   beta_var = q_beta_var,
                   mu_mean = q_mu_mean,
                   alpha_mean = q_alpha_mean,
                   beta_mean = q_beta_mean,
                   f = f)
    
    save(results, file = paste("data/results_basic_logit_p",period,".RData",sep=""))
  }
}
