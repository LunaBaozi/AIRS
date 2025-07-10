# SCRIPT FOCUSED ON AURKB NO BS

print(getwd())
setwd("/home/luna/Documents/Coding/GraphBP/GraphBP/analyses/results")
print(getwd())

aurka <- read.csv("tanimoto_simil_AURKA.csv")
aurkb <- read.csv("tanimoto_simil_AURKB.csv")
aurka_nobs <- read.csv("tanimoto_simil_AURKA_NOBS.csv")
aurkb_nobs <- read.csv("tanimoto_simil_AURKB_NOBS.csv")

# Basic file statistics
print(is.data.frame(aurkb_nobs))
print(ncol(aurkb_nobs))
print(nrow(aurkb_nobs))

# Basic analyses
summary(aurkb_nobs$simil)
summary(aurkb_nobs$cinfony.simil)
max_sim_morgan <- max(aurkb_nobs$simil)
max_sim_cinfony <- max(aurkb_nobs$cinfony.simil)
print(max_sim_morgan)
print(max_sim_cinfony)

# Sort the files to obtain highest scoring small molecules
order_aurkb_nobs_simil <- aurkb_nobs[order(aurkb_nobs$simil, decreasing = TRUE, na.last = TRUE),]
order_aurkb_nobs_cinfony <- aurkb_nobs[order(aurkb_nobs$cinfony.simil, decreasing = TRUE, na.last = TRUE),]
print(order_aurkb_nobs_simil[1:5,])
print(order_aurkb_nobs_cinfony[1:5,])
top5_aurkb_nobs_simil <- order_aurkb_nobs_simil[1:5,]
top5_aurkb_nobs_cinfony <- order_aurkb_nobs_cinfony[1:5,]

write.csv(top5_aurkb_nobs_simil, "top5_aurkb_nobs_simil", row.names = TRUE)
write.csv(top5_aurkb_nobs_cinfony, "top5_aurkb_nobs_cinfony", row.names = TRUE)


# Get the couple details having max similarity
retval1 <- subset(aurkb_nobs, simil == max_sim_morgan)
print(retval1)
retval2 <- subset(aurkb_nobs, cinfony.simil == max_sim_cinfony)
print(retval2)

# Start plotting something
x1 <- aurkb_nobs$simil
x2 <- aurkb_nobs$cinfony.simil
c1 <- rgb(204, 255, 122, max = 255, alpha = 80, names = "lt.green")
c2 <- rgb(255, 189, 122, max = 255, alpha = 80, names = "lt.orange")

histA <- hist(aurkb_nobs$simil, 30, plot=FALSE)
histB <- hist(aurkb_nobs$cinfony.simil, 50, plot=FALSE)

plot(histA, col = c1, xlim = c(-0.01, 0.250), ylim = c(0, 850), main = "Distribution of similarity scores (AurK A - NO BS)", xlab = "Tanimoto similarity")
plot(histB, col = c2, xlim = c(-0.01, 0.250), ylim = c(0, 850), add = TRUE) 
legend("topright", c("Morgan", "Cinfony"), lwd=2, col=c("green", "orange"))

xfit1 <- seq(min(x1), max(x1), length = length(x1))
yfit1 <- dnorm(xfit1, mean = mean(x1), sd = sd(x1))
yfit1 <- yfit1 * diff(histA$mids[1:2]) * length(x1)
xfit2 <- seq(min(x2), max(x2), length = length(x2))
yfit2 <- dnorm(xfit2, mean = mean(x2), sd = sd(x2))
yfit2 <- yfit2 * diff(histB$mids[1:2]) * length(x2)
lines(xfit1, yfit1, col = "green", lwd = 2) 
lines(xfit2, yfit2, col = "orange", lwd = 2) 
abline(v = c(mean(x1), mean(x2)), col=c("green", "orange"), lty=c(2,2), lwd=c(1, 1))

