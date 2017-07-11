N <- 50
psi <- 0.90
foo <- function(W, col="red") lines(0:N, dbb(0:N, N, W*(1-psi), W*psi), col=col)
plot(0:N, dbinom(0:N, N, 1-psi), type='l')
foo(N+2)
foo(4*N, 'blue')


rm(N, psi ,foo)

G <- 10000
N1 <- 30
N2 <- 50
tailRankCutoff(G, N1, N2, 0.95, 0.50, model="bin")
tailRankCutoff(G, N1, N2, 0.95, 0.50, model="bb")

tailRankCutoff(G, N1, N2, 0.95, 0.950, model="bin")
tailRankCutoff(G, N1, N2, 0.95, 0.950, model="bb")

tailRankPower(G, N1, N2, 0.95, 0.25, 0.50, "bin")
tailRankPower(G, N1, N2, 0.95, 0.25, 0.90, "bin")
tailRankPower(G, N1, N2, 0.95, 0.25, 0.50, "bb")
tailRankPower(G, N1, 120, 0.95, 0.25, 0.50, "bb")
