# install.packages("circular")
library(circular)

mu=0
kappa=c(1,10)
n=1000

x1=rvonmises(n, mu, kappa[1])
x2=rvonmises(n, mu, kappa[2])

par(mfrow=c(1,2))
hist(as.vector(x1), 
     breaks=50, main=paste("kappa =", kappa[1]), xlab="Angle (radians)", col="lightblue")
hist(as.vector(x2),
     breaks=50, main=paste("kappa =", kappa[2]), xlab="Angle (radians)", col="pink")

mu=pi/2
kappa=c(1,10)
n=1000

x1=rvonmises(n, mu, kappa[1])
x2=rvonmises(n, mu, kappa[2])

par(mfrow=c(1,2))
hist(as.vector(x1), 
     breaks=50, main=paste("kappa =", kappa[1]), xlab="Angle (radians)", col="lightblue")
hist(as.vector(x2),
     breaks=50, main=paste("kappa =", kappa[2]), xlab="Angle (radians)", col="pink")