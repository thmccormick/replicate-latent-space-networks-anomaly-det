#Setup
######

library(rjson)

ids <- read.csv("data/new_id_dictionary.csv")

n = max(ids$new_ids)


# Hyperparameter
init_threshold = log(0.5)





anomalies = list()
counter = 1

for(day in 2:90) {
  
  load(paste("data/scores_d",day,".RData",sep=""))
  data$score = data$score_logit
  #data$score = data$score_logit_base
  
  if(day == 2) {
    periods = 9:11
  } else {
    periods = 1:6 + day*6 - 7 
  }
  
  for(period in periods) {
    temp_data = data[data$TimeFrame == period & data$score < init_threshold,]
    
    src_counts <- table(temp_data$SrcDevice)
    
    # star3 (a-b;a-c;a-d)
    # record the most anomalous star for each source (a)
    sources = as.numeric(attr(which(src_counts >= 3), "names"))
    for(source in sources) {
      indices = order(temp_data$score[temp_data$SrcDevice == source])[1:3]
      anomalies[[counter]] = list(type = "star3",
                                  period = period,
                                  a = source,
                                  b = temp_data$DstDevice[temp_data$SrcDevice == source][indices[1]],
                                  c = temp_data$DstDevice[temp_data$SrcDevice == source][indices[2]],
                                  d = temp_data$DstDevice[temp_data$SrcDevice == source][indices[3]],
                                  score = sum(temp_data$score[temp_data$SrcDevice == source][indices]))
      counter = counter + 1
    }
    
    # fork (a-b;b-c;b-d)
    # record the most anomalous fork for each source (b)
    sources = as.numeric(attr(which(src_counts >= 2), "names"))
    sources = sources[sources %in% temp_data$DstDevice]
    for(source in sources) {
      indices1 = which.min(temp_data$score[temp_data$DstDevice == source])
      head = temp_data$SrcDevice[temp_data$DstDevice == source][indices1]
      indices2 = order(temp_data$score[temp_data$SrcDevice == source & temp_data$DstDevice != head])[1:2]
      
      if(any(is.na(indices2))) {
        next()
      }
      anomalies[[counter]] = list(type = "fork", 
                                  period = period,
                                  a = head,
                                  b = source,
                                  c = temp_data$DstDevice[temp_data$SrcDevice == source & temp_data$DstDevice != head][indices2[1]],
                                  d = temp_data$DstDevice[temp_data$SrcDevice == source & temp_data$DstDevice != head][indices2[2]],
                                  score = temp_data$score[temp_data$DstDevice == source][indices1] + sum(temp_data$score[temp_data$SrcDevice == source][indices2]))
      
      counter = counter + 1
    }
    
    # 3-path (a-b;b-c;c-d)
    # record the most anomalous fork for each pair b-c
    sources = unique(temp_data$SrcDevice)
    sources = sources[sources %in% temp_data$DstDevice]
    rows = which(temp_data$SrcDevice %in% sources & temp_data$DstDevice %in% sources)
    for (row in rows) {
      indices1 = which.min(temp_data$score[temp_data$DstDevice == temp_data$SrcDevice[row] & temp_data$SrcDevice != temp_data$DstDevice[row]])
      
      if(length(indices1) == 0) {
        next()
      }
      head = temp_data$SrcDevice[temp_data$DstDevice == temp_data$SrcDevice[row] & temp_data$SrcDevice != temp_data$DstDevice[row]][indices1]
      
      indices2 = which.min(temp_data$score[temp_data$SrcDevice == temp_data$DstDevice[row] & temp_data$DstDevice != head & temp_data$DstDevice != temp_data$SrcDevice[row]])
      
      if(length(indices2) == 0) {
        next()
      }
      tail = temp_data$DstDevice[temp_data$SrcDevice == temp_data$DstDevice[row] & temp_data$DstDevice != head & temp_data$DstDevice != temp_data$SrcDevice[row]][indices2]
      
      anomalies[[counter]] = list(type = "3path", 
                                  period = period,
                                  a = head,
                                  b = temp_data$SrcDevice[row],
                                  c = temp_data$DstDevice[row],
                                  d = tail,
                                  score = temp_data$score[row] + temp_data$score[temp_data$DstDevice == temp_data$SrcDevice[row] & temp_data$SrcDevice != temp_data$DstDevice[row]][indices1] + temp_data$score[temp_data$SrcDevice == temp_data$DstDevice[row] & temp_data$DstDevice != head & temp_data$DstDevice != temp_data$SrcDevice[row]][indices2])
      
      counter = counter + 1
    }
  }
}

save(anomalies, file = "data/anomalies_logit.RData")
#save(anomalies, file = "data/anomalies_logit_base.RData")