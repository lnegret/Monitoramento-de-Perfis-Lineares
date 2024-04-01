##############################################################################
a0 <- 2
b0 <- 1
s20 <- 0.5^2
alpha <- 0.005
j <- 1:20
n <- 5
m <- 20
datos <- c(2.209, 2.308, 1.914, 2.500, 3.147, 2.395, 1.796, 2.418, 1.987, 2.872,
           2.751, 2.778, 2.481, 3.294, 3.371, 2.245, 1.923, 2.502, 3.263, 3.241,
           1.902, 1.307, 2.263, 1.740, 2.367, 2.013, 2.056, 2.164, 2.749, 2.873,
           1.273, 2.361, 3.084, 2.892, 2.310, 1.482, 2.581, 1.720, 2.638, 2.674,
           2.743, 2.019, 2.186, 3.217, 2.516, 2.186, 2.516, 2.449, 2.461, 3.328,
           2.900, 3.033, 3.984, 4.469, 3.753, 3.493, 2.749, 3.566, 3.077, 3.645,
           2.481, 2.972, 2.885, 3.570, 4.033, 3.708, 3.568, 3.059, 3.681, 2.877,
           3.290, 2.485, 2.776, 3.291, 2.755, 3.686, 2.460, 3.085, 2.874, 4.261,
           3.396, 3.089, 3.656, 3.758, 3.720, 3.179, 3.530, 4.410, 2.808, 3.463,
           2.892, 3.104, 3.335, 3.978, 3.797, 2.390, 2.397, 3.746, 3.474, 4.114)
y <- matrix(datos, m, n, byrow = T)
y<-Ejemplo
x <- c(0.2, 0.4, 0.6, 0.8, 1)
###############################################################################
a1 <- 6
b1 <- 9
y1 <- t(a1+b1*x+t(datos))
y[11:20,]<-y1[11:20,]
a0 <- 5
b0 <- 10
s20 <- 0.5^2
alpha <- 0.005
j <- 1:20
n <- 5
m <- 20
x <- c(0.2, 0.4, 0.6, 0.8, 1)
datos <- matrix(rnorm(20*5,0,s20), m, n)
y <- t(a0+b0*x+t(datos))
###############################################################################
ymean <- rowMeans(y)
xmean <- mean(x)
best <- colSums((x-xmean)*t(y-ymean))/sum((x-xmean)^2)
aest <- ymean-best*xmean
s2est <- rowSums((y-(aest+t(t(replicate(n,best))*x)))^2)/(n-2)
sxx <- sum((x-xmean)^2)
d0 <- c(a0, b0)
Sigma <- s20*matrix(c(1/n+xmean^2/sxx, -xmean/sxx, -xmean/sxx, 1/sxx),2,2)
T2 <- c()
for (pos in 1:m) {
  dest <- c(aest[pos], best[pos])
  T2 <- c(T2, t(dest-d0)%*%solve(Sigma)%*%(dest-d0))
}
LCS <- qchisq(1-alpha, 2)
plot(j, aest, type = 'p', pch = 20, cex = 1.5)
plot(j, best, type = 'p', pch = 20, cex = 1.5)
plot(j, s2est, type = 'p', pch = 20, cex = 1.5)
plot(j, T2, type = 'o', ylim = c(0,max(T2,LCS)), pch = 20, cex = 1.5)
abline(h = LCS, lty = 3)
##########################################################################
aast <- aest+best*xmean
s0 <- sqrt(s20)
lambdai <- 0.2
rhoi <- 3.0156
Ei <- a0+b0*xmean
for (pos in 1:m) {
  Ei <- c(Ei, lambdai*aast[pos]+(1-lambdai)*Ei[pos])
}
Ei <- Ei[2:(m+1)]
U <- a0+b0*xmean+rhoi*s0*sqrt(lambdai/((2-lambdai)*n))
C <- a0+b0*xmean
L <- a0+b0*xmean-rhoi*s0*sqrt(lambdai/((2-lambdai)*n))
plot(j, Ei, type = 'o', ylim = c(min(min(Ei),L),max(max(Ei),U)), pch = 20,
     cex = 1.5)
abline(h = U, lty = 3)
abline(h = C)
abline(h = L, lty = 3)

lambdas <- 0.2
rhos <- 3.0109
Es <- b0
for (pos in 1:m) {
  Es <- c(Es, lambdas*best[pos]+(1-lambdas)*Es[pos])
}
Es <- Es[2:(m+1)]
U <- b0+rhos*s0*sqrt(lambdas/((2-lambdas)*sxx))
C <- b0
L <- b0-rhos*s0*sqrt(lambdas/((2-lambdas)*sxx))
plot(j, Es, type = 'o', ylim = c(min(min(Es),L),max(max(Es),U)), pch = 20, cex = 1.5)
abline(h = U, lty = 3)
abline(h = C)
abline(h = L, lty = 3)

lambdav <- 0.2
rhov <- 1.3723
Ev <- 0
for (pos in 1:m) {
  Ev <- c(Ev, max(0, lambdav*log(s2est[pos]/s20)+(1-lambdav)*Ev[pos]))
}
Ev <- Ev[2:(m+1)]
sw <- sqrt(2/(n-2)+2/(n-2)^2+4/(3*(n-2)^3)-16/(15*(n-2)^5))
U <- rhov*sw*sqrt(lambdav/(2-lambdav))
plot(j, Ev, type = 'o', ylim = c(0,max(max(Ev),U)), pch = 20, cex = 1.5)
abline(h = U, lty = 3)
##########################################################################
a0est <- mean(aest)
a0astest <- mean(ymean)
s20est <- mean(s2est)
s0est <- sqrt(s20est)
U <- a0astest+qt(1-alpha/2,(n-2)*m)*s0est*sqrt((m-1)/(m*n))
C <- a0astest
L <- a0astest-qt(1-alpha/2,(n-2)*m)*s0est*sqrt((m-1)/(m*n))
plot(j, ymean, type = 'o', ylim = c(min(min(ymean),L),max(max(ymean),U)), pch = 20, cex = 1.5)
abline(h = U, lty = 3)
abline(h = C)
abline(h = L, lty = 3)

b0est <- mean(best)
U <- b0est+qt(1-alpha/2,(n-2)*m)*s0est*sqrt((m-1)/m/sxx)
C <- b0est
L <- b0est-qt(1-alpha/2,(n-2)*m)*s0est*sqrt((m-1)/m/sxx)
plot(j, best, type = 'o', ylim = c(min(min(best),L),max(max(best),U)), pch = 20, cex = 1.5)
abline(h = U, lty = 3)
abline(h = C)
abline(h = L, lty = 3)

Datos = data.frame(1:20,round(aest,3),round(best,3),round(s2est,3),round(T2,3))
write.table(Datos, file="Ejemplo2.txt", row.names=F)
T2

Datos = data.frame(1:20,round(y,3))
write.table(Datos, file="Ejemplo.txt", row.names=F)

