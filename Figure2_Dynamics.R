# jonashaslbeck@gmail.com; October 2018

# Here we run the simulation shown in Figure 2, and plot Figure 2

#---------------------------------------------------------------------------------------- 
# ---------------- Aux Function: Ising Sampler ------------------------------------------
# ---------------------------------------------------------------------------------------

# Samples from a dynamic Ising model

IsSam <- function(G, 
                  th, 
                  n, 
                  domain = c(0,1)) {
  
  
  # Get aux variables
  p <- ncol(G)
  data <- matrix(NA, ncol = p, nrow = n)
  D <- domain
  
  # Set initial values
  data[1, ] <- sample(D, size = p, replace = T)
  
  # Loop over n cases  
  
  for(i in 2:n) {
    for(v in 1:p) {
      
      # Compute Potentials
      pot_A <- exp(sum(th[v]*D[1] +  G[v, -v] * D[1] * data[i-1, -v] ))
      pot_B <- exp(sum(th[v]*D[2] +  G[v, -v] * D[2] * data[i-1, -v] ))
      
      # Compute Probabilities
      Z <- pot_A + pot_B
      prob <- c(pot_A/Z, pot_B/Z)
      
      # Sample
      data[i, v] <- sample(D, size = 1, prob = prob)
      
    } # end for v
  } # end for i
  
  return(data)
  
} # eoF


#---------------------------------------------------------------------------------------- 
# ---------------- Sample Data ----------------------------------------------------------
# ---------------------------------------------------------------------------------------

# Setup variations
p <- 10
theta <- c(0, -.1, -.2)
th <- rep(0, p)
n <- 1000000

# Storage
d_all <- vector("list", length=2)
d_all <- lapply(d_all, function(x) vector("list", length=3))
l_G <- list()

set.seed(1)
for(w in 1:3) {
  l_G[[w]] <- G <- matrix(theta[w], p, p)
  for(d in 1:2) {
    ifelse(d==1, domain <- c(0,1), domain <- c(-1, 1))
    d_all[[d]][[w]] <- IsSam(G = G,
                             th = th, 
                             n = n, 
                             domain =domain)
  }
  print(w)
}

#---------------------------------------------------------------------------------------- 
# ---------------- Plot Figure 2 --------------------------------------------------------
# ---------------------------------------------------------------------------------------

sc <- 8

pdf("Dynamics_Illu_Figure_negative.pdf", width = sc, height = sc)

lmat <- cbind(1:3,4:6, 7:9)
lo <- layout(lmat)

# plot the graphs
library(qgraph)
maximum <- .4
for(w in 1:3 ) qgraph(l_G[[w]], maximum = maximum, 
                      layout = "circle", labels = FALSE, 
                      mar=c(5, 5, 5, 5), edge.color = "black")


# for {-1,1}
for(w in 1:3) {
  rs <- apply(d_all[[2]][[w]], 1, function(x) sum(x == 1))
  hist(rs, xlim=c(0, 10), main = "", ylim=c(0, .75), 
       probability  = T, xlab = "Number Variables in state {1}", 
       xaxt = "n", yaxt = "n", breaks = 20)
  axis(1, c(0, 5, 10))
  axis(2, c(0, .25, .5, .75), las=2)
  
  text(1, .75, paste0("Mean = ", round(mean(rs),2)), adj=0)
  text(1, .65, paste0("SD = ", round(sd(rs),2)), adj=0)
  
}

# for {0,1}
for(w in 1:3) {
  rs <- apply(d_all[[1]][[w]], 1, function(x) sum(x == 1))
  hist(rs, xlim=c(0, 10), main = "", ylim=c(0, .75), 
       probability  = T, xlab = "Number Variables in state {1}", 
       xaxt = "n", yaxt = "n", breaks = 20)
  axis(1, c(0, 5, 10))
  axis(2, c(0, .25, .5, .75), las=2)
  
  text(1, .75, paste0("Mean = ", round(mean(rs),2)), adj=0)
  text(1, .65, paste0("SD = ", round(sd(rs),2)), adj=0)
}


dev.off()





