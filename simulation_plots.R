library(pROC)

load("data/simulation_logit_data.RData")

periods = 100
ll = lapply(1:periods, function(p) sum(y[[p]]*log(expit(mean_func[[p]])) + (1-y[[p]])*log(1-expit(mean_func[[p]]))))
ll_full <- list()
ll_cc <- list()
auc_true <- list()
auc_full <- list()
auc_cc <- list()

corr_full <- list()
corr_cc <- list()
corr_never_full <- list()
corr_never_cc <- list()

for(t in 1:periods) {
  load(paste("data/simulation_logit_p",t,".RData",sep=""))
  results_full = results
  mf_full = rep(results_full$mu_mean, times = n*n) + rep(results_full$alpha_mean, times = n) + rep(results_full$beta_mean, each = n) + as.numeric(results_full$u_mean %*% t(results_full$v_mean))
  
  load(paste("data/simulation_logit_cc_p",t,".RData",sep=""))
  results_cc = results
  mf_cc = rep(results_cc$mu_mean+log(0.025), times = n*n) + rep(results_cc$alpha_mean, times = n) + rep(results_cc$beta_mean, each = n) + as.numeric(results_cc$u_mean %*% t(results_cc$v_mean))
  
  corr_cc[[t]] = cor(mean_func[[t]],mf_cc,method = c("pearson"))
  corr_full[[t]] = cor(mean_func[[t]],mf_full,method = c("pearson"))
  
  corr_never_cc[[t]] = cor(mean_func[[t]][never_y],mf_cc[never_y],method = c("pearson"))
  corr_never_full[[t]] = cor(mean_func[[t]][never_y],mf_full[never_y],method = c("pearson"))
}



# Figure 2
##########

pdf("simulation_comp.pdf",8,3)
par(mfrow = c(1,2), mar = c(4.1,4.1,1.1,1.1))
plot(1:100,unlist(auc_true), type = "l", lwd = 2, ylab = "AUC ROC", xlab = "Time Period", ylim = range(unlist(auc_full)))
lines(1:100,unlist(auc_full), type = "l", lwd = 2, col = rgb(1,0,0,alpha=0.5))
lines(1:100,unlist(auc_cc), type = "l", lwd = 2, col = rgb(0,1,1,alpha=0.5))
legend("topright",legend = c("Actual", "Full EP", "CC EP"),col = c("black","red","cyan"),lwd = 2,cex = 0.5)

plot(1:100,unlist(ll), type = "l", lwd = 2, ylab = "Log-likelihood", xlab = "Time Period")
lines(1:100,unlist(ll_full), type = "l", lwd = 2, col = rgb(1,0,0,alpha=0.5))
lines(1:100,unlist(ll_cc), type = "l", lwd = 2, col = rgb(0,1,1,alpha=0.5))
legend("bottomleft",legend = c("Actual", "Full EP", "CC EP"),col = c("black","red","cyan"),lwd = 2,cex = 0.5)
dev.off()



# Figure 3
##########

pdf("simulation_comp_corr.pdf",8,3,)
par(mfrow = c(1,2), mar = c(4.1,4.1,1.1,1.1))
plot(1:100,unlist(corr_full), type = "l", lwd = 2, col = rgb(1,0,0,alpha=0.5), ylab = "Correlation", xlab = "Time Period", ylim = c(0.6,1))
lines(1:100,unlist(corr_cc), type = "l", lwd = 2, col = rgb(0,1,1,alpha=0.5))
legend("bottomright",legend = c("Full EP", "CC EP"),col = c("red","cyan"),lwd = 2,cex = 0.5)

plot(1:100,unlist(corr_never_full), type = "l", lwd = 2, col = rgb(1,0,0,alpha=0.5), ylab = "Cor. (unobserved only)", xlab = "Time Period", , ylim = c(0.5,1))
lines(1:100,unlist(corr_never_cc), type = "l", lwd = 2, col = rgb(0,1,1,alpha=0.5))
legend("bottomright",legend = c("Full EP", "CC EP"),col = c("red","cyan"),lwd = 2,cex = 0.5)
dev.off()



# Figure 4
##########

pdf("scatter_corr.pdf",8,9)
par(mfrow = c(3,2), mar = c(4.1,4.1,1.1,1.1))

t = 10
load(paste("data/simulation_logit_p",t,".RData",sep=""))
results_full = results
mf_full = rep(results_full$mu_mean, times = n*n) + rep(results_full$alpha_mean, times = n) + rep(results_full$beta_mean, each = n) + as.numeric(results_full$u_mean %*% t(results_full$v_mean))

load(paste("data/simulation_logit_cc_p",t,".RData",sep=""))
results_cc = results
mf_cc = rep(results_cc$mu_mean+log(0.025), times = n*n) + rep(results_cc$alpha_mean, times = n) + rep(results_cc$beta_mean, each = n) + as.numeric(results_cc$u_mean %*% t(results_cc$v_mean))

plot(mean_func[[t]],mf_full, col = rgb(0.5,0.5,0.5,alpha=0.1), ylab = "Full EP", xlab = "Actual", main = "T = 10")
abline(a = 0, b = 1, lwd = 2)
plot(mean_func[[t]],mf_cc, col = rgb(0.5,0.5,0.5,alpha=0.1), ylab = "CC EP", xlab = "Actual", main = "T = 10")
abline(a = 0, b = 1, lwd = 2)

t = 25
load(paste("data/simulation_logit_p",t,".RData",sep=""))
results_full = results
mf_full = rep(results_full$mu_mean, times = n*n) + rep(results_full$alpha_mean, times = n) + rep(results_full$beta_mean, each = n) + as.numeric(results_full$u_mean %*% t(results_full$v_mean))

load(paste("data/simulation_logit_cc_p",t,".RData",sep=""))
results_cc = results
mf_cc = rep(results_cc$mu_mean+log(0.025), times = n*n) + rep(results_cc$alpha_mean, times = n) + rep(results_cc$beta_mean, each = n) + as.numeric(results_cc$u_mean %*% t(results_cc$v_mean))

plot(mean_func[[t]],mf_full, col = rgb(0.5,0.5,0.5,alpha=0.1), ylab = "Full EP", xlab = "Actual", main = "T = 25")
abline(a = 0, b = 1, lwd = 2)
plot(mean_func[[t]],mf_cc, col = rgb(0.5,0.5,0.5,alpha=0.1), ylab = "CC EP", xlab = "Actual", main = "T = 25")
abline(a = 0, b = 1, lwd = 2)

t = 100
load(paste("data/simulation_logit_p",t,".RData",sep=""))
results_full = results
mf_full = rep(results_full$mu_mean, times = n*n) + rep(results_full$alpha_mean, times = n) + rep(results_full$beta_mean, each = n) + as.numeric(results_full$u_mean %*% t(results_full$v_mean))

load(paste("data/simulation_logit_cc_p",t,".RData",sep=""))
results_cc = results
mf_cc = rep(results_cc$mu_mean+log(0.025), times = n*n) + rep(results_cc$alpha_mean, times = n) + rep(results_cc$beta_mean, each = n) + as.numeric(results_cc$u_mean %*% t(results_cc$v_mean))

plot(mean_func[[t]],mf_full, col = rgb(0.5,0.5,0.5,alpha=0.1), ylab = "Full EP", xlab = "Actual", main = "T = 100")
abline(a = 0, b = 1, lwd = 2)
plot(mean_func[[t]],mf_cc, col = rgb(0.5,0.5,0.5,alpha=0.1), ylab = "CC EP", xlab = "Actual", main = "T = 100")
abline(a = 0, b = 1, lwd = 2)
dev.off()

