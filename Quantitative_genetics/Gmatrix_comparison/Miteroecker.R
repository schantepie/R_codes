# *~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*
# DISTANCES
# *~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*

D12= sqrt(sum(log(eigen(solve(pG_Age1,pG_Age2))$values)^2))
D13= sqrt(sum(log(eigen(solve(pG_Age1,pG_Age3))$values)^2))
D14= sqrt(sum(log(eigen(solve(pG_Age1,pG_Age4))$values)^2))
D15= sqrt(sum(log(eigen(solve(pG_Age1,pG_Age5))$values)^2))
D16= sqrt(sum(log(eigen(solve(pG_Age1,pG_Age6))$values)^2))
D17= sqrt(sum(log(eigen(solve(pG_Age1,pG_Age7))$values)^2))
D18= sqrt(sum(log(eigen(solve(pG_Age1,pG_Age8))$values)^2))

D23= sqrt(sum(log(eigen(solve(pG_Age2,pG_Age3))$values)^2))
D24= sqrt(sum(log(eigen(solve(pG_Age2,pG_Age4))$values)^2))
D25= sqrt(sum(log(eigen(solve(pG_Age2,pG_Age5))$values)^2))
D26= sqrt(sum(log(eigen(solve(pG_Age2,pG_Age6))$values)^2))
D27= sqrt(sum(log(eigen(solve(pG_Age2,pG_Age7))$values)^2))
D28= sqrt(sum(log(eigen(solve(pG_Age2,pG_Age8))$values)^2))

D34= sqrt(sum(log(eigen(solve(pG_Age3,pG_Age4))$values)^2))
D35= sqrt(sum(log(eigen(solve(pG_Age3,pG_Age5))$values)^2))
D36= sqrt(sum(log(eigen(solve(pG_Age3,pG_Age6))$values)^2))
D37= sqrt(sum(log(eigen(solve(pG_Age3,pG_Age7))$values)^2))
D38= sqrt(sum(log(eigen(solve(pG_Age3,pG_Age8))$values)^2))

D45= sqrt(sum(log(eigen(solve(pG_Age4,pG_Age5))$values)^2))
D46= sqrt(sum(log(eigen(solve(pG_Age4,pG_Age6))$values)^2))
D47= sqrt(sum(log(eigen(solve(pG_Age4,pG_Age7))$values)^2))
D48= sqrt(sum(log(eigen(solve(pG_Age4,pG_Age8))$values)^2))

D56= sqrt(sum(log(eigen(solve(pG_Age5,pG_Age6))$values)^2))
D57= sqrt(sum(log(eigen(solve(pG_Age5,pG_Age7))$values)^2))
D58= sqrt(sum(log(eigen(solve(pG_Age5,pG_Age8))$values)^2))

D67= sqrt(sum(log(eigen(solve(pG_Age6,pG_Age7))$values)^2))
D68= sqrt(sum(log(eigen(solve(pG_Age6,pG_Age8))$values)^2))

D78= sqrt(sum(log(eigen(solve(pG_Age7,pG_Age8))$values)^2))

DistMat=matrix(c(0,D12, D13,D14,D15,D16,D17,D18,D12,0, D23,D24,D25,D26,D27,D28,D13,D23,0,D34,D35,D36,D37,D38,D14,D24,D34,0,D45,D46,D47,D38, D15,D25,D35,D45,0,D56,D57,D58,D16,D26,D36,D46,D56,0,D67,D68,D17,D27,D37,D47,D57,D67,0,D78,D18,D28,D38,D48,D58,D68,D78,0), ncol=8)


DistMat=as.dist (DistMat)
library(ape)
res <- pcoa(DistMat)
res$values
biplot(res)

lines(res$vectors[,1],res$vectors[,2])





