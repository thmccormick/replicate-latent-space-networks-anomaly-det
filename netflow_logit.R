# Setup
#######

# helper functions
reltol <- function(x,y) {
  sum(abs(x-y))/sum(abs(x))
}

expit <- function(x) {y = exp(x); y/(1+y)}

weighted_sum <- function(x1, x2, w) 1/(1+w)*x1 + w/(1+w)*x2


# linear algebra helper (speed) functions for ldim=2 case
mdet <- function(mat) {
  return(mat[1,1]*mat[2,2] - mat[1,2]*mat[2,1])
}

msolve <- function(mat) {
  return({mat = c(mat[2,2], -mat[2,1], -mat[1,2], mat[1,1])/(mat[1,1]*mat[2,2] - mat[1,2]*mat[2,1]); dim(mat) = c(2,2); mat})
}

vmv <- function(m, v) {
  return(m[1,1]*v[1]^2 + (m[2,1] + m[1,2])*v[1]*v[2] + m[2,2]*v[2]^2)
}

vvt <- function(v) {
  return({mat = c(v[1]^2, v[1]*v[2], v[1]*v[2], v[2]^2); dim(mat) = c(2,2); mat})
}


# Auxillary data
ids <- read.csv("data/new_id_dictionary.csv")

n = max(ids$new_ids)

load("data/period_weights.RData")





# Estimation
############

sigma_prior = c(0.5,0.25,0.25,0)

damping = 2
epsilon = 3e-03
ldim = 2
average_zeros = 500000 #on average, sample 500000 zeros in each period


mat_prior = matrix(data = c(sigma_prior[3],sigma_prior[4],sigma_prior[4],sigma_prior[3]), nrow = 2)
sigma_prior_uv = array(0, c(ldim,ldim,n))
sigma_prior_uv[1,1,] = sigma_prior[3]
sigma_prior_uv[1,2,] = sigma_prior[4]
sigma_prior_uv[2,1,] = sigma_prior[4]
sigma_prior_uv[2,2,] = sigma_prior[3]


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
    num_zeros = round(500000*period_weight)
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
      q_u_mean <- matrix(data = rnorm(n*ldim,mean=0,sd=0.1), nrow = n, ncol = ldim)
      q_u_var = sigma_prior_uv
      q_v_mean <- matrix(data = rnorm(n*ldim,mean=0,sd=0.1), nrow = n, ncol = ldim)
      q_v_var <- q_u_var
      
      f = 0
    } else {
      forgetting <- as.matrix(expand.grid(0:3, 0:3, 0:3))
      colnames(forgetting) = NULL
      forgetting[forgetting == 0] = 1.1
      forgetting[forgetting == 1] = 1
      forgetting[forgetting == 2] = 1.01
      forgetting[forgetting == 3] = 2
      
      # Estimate amount of forgetting for period
      ll = rep(0, dim(forgetting)[1])
      for(row in 1:n) {
        i = rows[row,1]
        j = rows[row,2]
        sign = 1-2*rows[row,3]

        mean_f = q_mu_mean + q_alpha_mean[i] + q_beta_mean[j] + t(q_u_mean[i,]) %*% q_v_mean[j,]
        temp1 = sum((q_u_var[,,i]) * (q_v_var[,,j]))
        temp2 = vmv(q_v_var[,,j], q_u_mean[i,])
        temp3 = vmv(q_u_var[,,i], q_v_mean[j,])
        for(row2 in 1:dim(forgetting)[1]) {
          f = forgetting[row2,]
          sigma_f = q_mu_var*f[1] + q_alpha_var[i]*f[2] + q_beta_var[j]*f[2] + f[3]^2*temp1 + f[3]*temp2 + f[3]*temp3
          
          ll[row2] = ll[row2] + log(expit(-sign*(1+pi/8*sigma_f)^(-0.5)*mean_f))
        }
      }
      
      f = forgetting[which.max(ll),]
      print(f)
      
      q_mu_var <- q_mu_var*f[1]
      q_alpha_var <- q_alpha_var*f[2]
      q_beta_var <- q_beta_var*f[2]
      q_u_var <- q_u_var*f[3]
      q_v_var <- q_v_var*f[3]
    }
    
    # Initiate messages
    msg_mu_mean <- rep(0, dim(rows)[1])
    msg_mu_var <- rep(Inf, dim(rows)[1])
    msg_alpha_mean <- rep(0, dim(rows)[1])
    msg_alpha_var <- rep(Inf, dim(rows)[1])
    msg_beta_mean <- rep(0, dim(rows)[1])
    msg_beta_var <- rep(Inf, dim(rows)[1])
    msg_u_mean <- array(0, c(ldim,dim(rows)[1]))
    msg_u_prec <- array(0, c(ldim,ldim,dim(rows)[1]))
    msg_u_prec[1,1,] = 1e-10
    msg_u_prec[ldim,ldim,] = 1e-10
    msg_v_mean <- msg_u_mean
    msg_v_prec <- msg_u_prec
    
    counter = 0
    
    while(TRUE) {
      q_u_old_mean = q_u_mean
      q_v_old_mean = q_v_mean
      q_u_old_var = q_u_var
      q_v_old_var = q_v_var
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
        q_tmp_u_mean = q_u_mean[i,]
        tmp_inv_q_u = msolve(q_u_var[,,i])
        q_tmp_v_mean = q_v_mean[j,]
        tmp_inv_q_v = msolve(q_v_var[,,j])
        
        msg_tmp_mu_mean = msg_mu_mean[row]
        msg_tmp_mu_var = msg_mu_var[row]
        msg_tmp_alpha_mean = msg_alpha_mean[row]
        msg_tmp_alpha_var = msg_alpha_var[row]
        msg_tmp_beta_mean = msg_beta_mean[row]
        msg_tmp_beta_var = msg_beta_var[row]
        msg_tmp_u_mean = msg_u_mean[,row]
        msg_tmp_u_prec = msg_u_prec[,,row]
        msg_tmp_v_mean = msg_v_mean[,row]
        msg_tmp_v_prec = msg_v_prec[,,row]
        
        # Define "g functions" (see bottom of page 8)
        g_mu_var = (1/q_mu_var + 1/msg_mu_var[row])^(-1)
        g_mu_mean = g_mu_var*(q_mu_mean/q_mu_var + msg_mu_mean[row]/msg_mu_var[row])
        g_alpha_var = (1/q_alpha_var[i] + 1/msg_alpha_var[row])^(-1)
        g_alpha_mean = g_alpha_var*(q_alpha_mean[i]/q_alpha_var[i] + msg_alpha_mean[row]/msg_alpha_var[row])
        g_beta_var = (1/q_beta_var[j] + 1/msg_beta_var[row])^(-1)
        g_beta_mean = g_beta_var*(q_beta_mean[j]/q_beta_var[j] + msg_beta_mean[row]/msg_beta_var[row])
        g_u_prec = tmp_inv_q_u + msg_u_prec[,,row]
        g_u_var = msolve(g_u_prec)
        g_u_mean = g_u_var %*% (tmp_inv_q_u %*% q_u_mean[i,] + msg_u_prec[,,row] %*% msg_u_mean[,row])
        g_v_prec = tmp_inv_q_v + msg_v_prec[,,row]
        g_v_var = msolve(g_v_prec)
        g_v_mean = g_v_var %*% (tmp_inv_q_v %*% q_v_mean[j,] + msg_v_prec[,,row] %*% msg_v_mean[,row])
        
        tmp_mat = sign*g_u_mean + g_v_prec %*% g_v_mean
        tmp_mat2 = sign*g_v_mean + g_u_prec %*% g_u_mean
        tmp_var = msolve(g_v_prec - g_u_var)
        tmp_var2 = msolve(g_u_prec - g_v_var)
        tmp_weight = mdet(diag(ldim) - g_u_var %*% g_v_var)^(-0.5) * exp(sign*(g_mu_mean + g_alpha_mean + g_beta_mean) + g_mu_var/2 + g_alpha_var/2 + g_beta_var/2 + 0.5*vmv(tmp_var,tmp_mat) - 0.5*vmv(g_v_prec,g_v_mean))
        
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
        
        tmp_mean = g_u_mean
        tmp_mean2 = tmp_var2 %*% tmp_mat2
        q_prime_u_mean = weighted_sum(tmp_mean,tmp_mean2,tmp_weight)
        q_prime_u_var = weighted_sum(vvt(tmp_mean)+g_u_var,vvt(tmp_mean2)+tmp_var2,tmp_weight) - vvt(q_prime_u_mean)
        
        tmp_mean = g_v_mean
        tmp_mean2 = tmp_var %*% tmp_mat
        q_prime_v_mean = weighted_sum(tmp_mean,tmp_mean2,tmp_weight)
        q_prime_v_var = weighted_sum(vvt(tmp_mean)+g_v_var,vvt(tmp_mean2)+tmp_var,tmp_weight) - vvt(q_prime_v_mean)

        # Equation 9
        q_mu_var = (damping/q_tmp_mu_var + (1-damping)/q_prime_mu_var)^(-1)
        q_mu_mean = q_mu_var*(damping/q_tmp_mu_var*q_tmp_mu_mean + (1-damping)/q_prime_mu_var*q_prime_mu_mean)
        q_alpha_var[i] = (damping/q_tmp_alpha_var + (1-damping)/q_prime_alpha_var)^(-1)
        q_alpha_mean[i] = q_alpha_var[i]*(damping/q_tmp_alpha_var*q_tmp_alpha_mean + (1-damping)/q_prime_alpha_var*q_prime_alpha_mean)
        q_beta_var[j] = (damping/q_tmp_beta_var + (1-damping)/q_prime_beta_var)^(-1)
        q_beta_mean[j] = q_beta_var[j]*(damping/q_tmp_beta_var*q_tmp_beta_mean + (1-damping)/q_prime_beta_var*q_prime_beta_mean)
        tmp_inv = msolve(q_prime_u_var)
        tmp_inv_u = damping*tmp_inv_q_u + (1-damping)*tmp_inv
        q_u_var[,,i] = msolve(tmp_inv_u)
        q_u_mean[i,] = q_u_var[,,i] %*% (damping*tmp_inv_q_u%*%q_tmp_u_mean + (1-damping)*tmp_inv%*%q_prime_u_mean)
        tmp_inv = msolve(q_prime_v_var)
        tmp_inv_v = damping*tmp_inv_q_v + (1-damping)*tmp_inv
        q_v_var[,,j] = msolve(tmp_inv_v)
        q_v_mean[j,] = q_v_var[,,j] %*% (damping*tmp_inv_q_v%*%q_tmp_v_mean + (1-damping)*tmp_inv%*%q_prime_v_mean)
        
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
        msg_u_prec[,,row] = tmp_inv_u + msg_tmp_u_prec - tmp_inv_q_u
        msg_u_mean[,row] = msolve(msg_u_prec[,,row]) %*% (tmp_inv_u %*% q_u_mean[i,] + msg_tmp_u_prec %*% msg_u_mean[,row] - tmp_inv_q_u %*% q_tmp_u_mean)
        msg_v_prec[,,row] = tmp_inv_v + msg_tmp_v_prec - tmp_inv_q_v
        msg_v_mean[,row] = msolve(msg_v_prec[,,row]) %*% (tmp_inv_v %*% q_v_mean[j,] + msg_tmp_v_prec %*% msg_v_mean[,row] - tmp_inv_q_v %*% q_tmp_v_mean)
        
      }
      
      criteria = reltol(c(q_u_var,q_v_var), c(q_u_old_var,q_v_old_var))
      print(paste("Period: ",period,", Iter: ",counter,", Time: ",format(round(difftime(Sys.time(),start,units="mins"),3),nsmall=4),", Rel Tol: ", format(criteria, nsmall=4),sep=""))
      
      if(criteria < epsilon || heart_break == TRUE) {
        break
      }
      counter = counter+1
    }
    
    results = list(mu_var = q_mu_var, 
                   alpha_var = q_alpha_var,
                   beta_var = q_beta_var,
                   u_var = q_u_var,
                   v_var = q_v_var,
                   mu_mean = q_mu_mean,
                   alpha_mean = q_alpha_mean,
                   beta_mean = q_beta_mean,
                   u_mean = q_u_mean,
                   v_mean = q_v_mean,
                   f = f)
    
    save(results, file = paste("data/results_logit_p",period,".RData",sep=""))
  }
}
