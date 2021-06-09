library(pROC)
setwd("C:/Users/Yanran/Desktop/Network anomaly detection/")

expit <- function(x) {y = exp(x); y/(1+y)}

load("dcsbm/simulation_dcsbm_data.RData")

periods = 100
n = 500

auc_dcsbm <- list()
auc_2latent <- list()
corr_dcsbm <- list()

for(t in 1:periods) {
  #epsilon 1.5
  load(paste("dcsbm/dcsbm_logit_p",t,".RData",sep=""))
  results_full = results
  mf = rep(results_full$mu_mean, times = n*n) + rep(results_full$alpha_mean, times = n) + rep(results_full$beta_mean, each = n) + as.numeric(results_full$u_mean %*% t(results_full$v_mean))
  never_y = which(y[[t]]==0)
  
  auc_roc = roc(as.factor(y[[t]]), mf)
  auc_dcsbm[[t]] = auc(auc_roc)
}

load("data/simulation_logit_data.RData")

for(t in 1:periods) {
  #epsilon 1.5
  load(paste("data/simulation_logit_p",t,".RData",sep=""))
  results_full = results
  mf = rep(results_full$mu_mean, times = n*n) + rep(results_full$alpha_mean, times = n) + rep(results_full$beta_mean, each = n) + as.numeric(results_full$u_mean %*% t(results_full$v_mean))
  never_y = which(y[[t]]==0)
  
  roc_2latent = roc(as.factor(y[[t]]), mf)
  auc_2latent[[t]] = auc(roc_2latent)
}
# Figure for AUC ROC comparison and Log-likelihood comparison
##########

pdf("dcsbm/simulation_dcsbm.pdf",8,3)
par(mfrow = c(1,1), mar = c(4.1,4.1,1.1,1.1))
plot(1:100,unlist(auc_dcsbm), type = "l", lty=1, lwd = 2, ylab = "AUC ROC", xlab = "Time Period", ylim = range(0.85,0.94))
lines(1:100,unlist(auc_2latent), type = "l", lty=2, lwd = 2, col = rgb(1,0,0,alpha=0.5))
legend("topright",legend = c("DCSBM", "Latent"), col = c("black","red"),lty=c(1,2),lwd = 2,cex = 0.8)
dev.off()
