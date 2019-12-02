# Measure weekly periodicity in netflow data
############################################

# Auxillary data
ids <- read.csv("data/new_id_dictionary.csv")
n = max(ids$new_ids)



# Count events in each period
obs_counts <- list()
counter = 0
for(period in 2:90) {
  data <- read.csv(paste("data/netflow_re4h-d",formatC(period,width=2,flag="0"),".csv",sep=""))
  tbl = table(data$TimeFrame)
  for(i in 1:length(tbl)) {
    obs_counts[[counter + i]] = tbl[i]
  }
  counter = counter + i
}

times = sapply(obs_counts, function(x) as.numeric(attr(x,"names")))
counts = as.numeric(sapply(obs_counts, function(x) x))

activity_4h <- data.frame(period = times, counts = counts)
activity_4h[nrow(activity_4h)+1,] = c(372,0) #period with no activity
activity_4h$tod = as.factor(activity_4h$period %% 6)
activity_4h$dow = as.factor((activity_4h$period %/% 6) %% 7)



# Estimate mean activity in each tod x dow period
period_weights <- expand.grid(tod = 0:5, dow = 0:6)
period_weights$mean_activity = 0
for(i in 1:nrow(period_weights)) {
  period_weights$mean_activity[i] = mean(activity_4h$counts[activity_4h$tod == period_weights$tod[i] & activity_4h$dow == period_weights$dow[i]])/n/n
}

period_weights$weight = mean(period_weights$mean_activity)/period_weights$mean_activity



# Case-control mean shift
odds = period_weights$mean_activity/(1-period_weights$mean_activity)
obs_odds = n*n*period_weights$mean_activity/round(500000/period_weights$weight)

period_weights$mu_mod <- log(odds/obs_odds)

save(period_weights, file = "data/period_weights.RData")
