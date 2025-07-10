# STARTING SCRIPT + FOCUS ON AURKA
print(getwd())
setwd("/home/luna/Documents/Coding/GraphBP/GraphBP/analyses/results")
print(getwd())

aurka <- read.csv("tanimoto_simil_AURKA_4byi.csv")
aurkb <- read.csv("tanimoto_simil_AURKB.csv")
aurka_nobs <- read.csv("tanimoto_simil_AURKA_NOBS.csv")
aurkb_nobs <- read.csv("tanimoto_simil_AURKB_NOBS.csv")

# Basic file statistics
print(is.data.frame(aurka))
print(ncol(aurka))
print(nrow(aurka))

# Basic analyses
summary(aurka$simil)
summary(aurka$cinfony.simil)
max_sim_morgan <- max(aurka$simil)
max_sim_cinfony <- max(aurka$cinfony.simil)
print(max_sim_morgan)
print(max_sim_cinfony)

# Sort the files to obtain highest scoring small molecules
order_aurka_simil <- aurka[order(aurka$simil, decreasing = TRUE, na.last = TRUE),]
order_aurka_cinfony <- aurka[order(aurka$cinfony.simil, decreasing = TRUE, na.last = TRUE),]
print(order_aurka_simil[1:5,])
print(order_aurka_cinfony[1:5,])
top5_aurka_simil <- order_aurka_simil[1:5,]
top5_aurka_cinfony <- order_aurka_cinfony[1:5,]

write.csv(top5_aurka_simil, "top5_aurka_simil.csv", row.names = TRUE)
write.csv(top5_aurka_cinfony, "top5_aurka_cinfony.csv", row.names = TRUE)


# Get the couple details having max similarity
retval1 <- subset(aurka, simil == max_sim_morgan)
print(retval1)
retval2 <- subset(aurka, cinfony.simil == max_sim_cinfony)
print(retval2)

# Start plotting something
x1 <- aurka$simil
x2 <- aurka$cinfony.simil
c1 <- rgb(173,216,230,max = 255, alpha = 80, names = "lt.blue")
c2 <- rgb(255,192,203, max = 255, alpha = 80, names = "lt.pink")

histA <- hist(aurka$simil, 30, plot=FALSE)
histB <- hist(aurka$cinfony.simil, 50, plot=FALSE)

plot(histA, col = c1, xlim = c(-0.01, 0.250), ylim = c(0, 850), main = "Distribution of similarity scores (AurK A)", xlab = "Tanimoto similarity")
plot(histB, col = c2, xlim = c(-0.01, 0.250), ylim = c(0, 850), add = TRUE) 
legend("topright", c("Morgan", "Cinfony"), lwd=2, col=c("blue", "red"))

xfit1 <- seq(min(x1), max(x1), length = length(x1))
yfit1 <- dnorm(xfit1, mean = mean(x1), sd = sd(x1))
yfit1 <- yfit1 * diff(histA$mids[1:2]) * length(x1)
xfit2 <- seq(min(x2), max(x2), length = length(x2))
yfit2 <- dnorm(xfit2, mean = mean(x2), sd = sd(x2))
yfit2 <- yfit2 * diff(histB$mids[1:2]) * length(x2)
lines(xfit1, yfit1, col = "blue", lwd = 2) 
lines(xfit2, yfit2, col = "red", lwd = 2) 
abline(v = c(mean(x1), mean(x2)), col=c("blue", "red"), lty=c(2,2), lwd=c(1, 1))

