# jonashaslbeck@gmail.com; November 2018

# Here we generate the data and estimate the models shown in Figure 1

# --------------------------------------------------------------------------------------- 
# ---------------- Generate two dependent binary variables ------------------------------
# ---------------------------------------------------------------------------------------

# ----- Data Generation -----

# Define joint distribution
probs <- matrix(c(.14, .18, .18, .5), 2, 2, byrow = T)

N <- 1000
data <- rbind(matrix(rep(c(0,0), probs[1,1]*N), ncol=2, byrow = T), 
              matrix(rep(c(1,0), probs[2,1]*N), ncol=2, byrow = T), 
              matrix(rep(c(0,1), probs[1,2]*N), ncol=2, byrow = T),
              matrix(rep(c(1,1), probs[2,2]*N), ncol=2, byrow = T))

tb <- table(data[, 1], data[, 2])
tb_n <- round(tb / sum(tb) , 3)

# To start out with:
data[data == 1] <- 2
data[data == 0] <- 1


# --------------------------------------------------------------------------------------- 
# ---------------- Estimate Different Ising Versions ------------------------------------
# ---------------------------------------------------------------------------------------

library(IsingSampler)

# ---------- Version I: A = 0, B = 1 ---------

data_i <- data
data_i[data_i == 1] <- 0
data_i[data_i == 2] <- 1

## Different Estimation Methods

# Loglinear model
fit_i <- EstimateIsing(data = data_i, responses = c(0,1), method = "ll")
fit_i$graph
fit_i$thresholds

# Look at empirical conditional probabilties
tb <- table(data_i[,1], data_i[,2])
tb/sum(tb)

# Compute potentials
pot_AA <- exp(fit_i$thresholds[1] * 0 + fit_i$thresholds[2] * 0 + fit_i$graph[2, 1] * 0 * 0)
pot_AB <- exp(fit_i$thresholds[1] * 0 + fit_i$thresholds[2] * 1 + fit_i$graph[2, 1] * 0 * 1)
pot_BA <- exp(fit_i$thresholds[1] * 1 + fit_i$thresholds[2] * 0 + fit_i$graph[2, 1] * 1 * 0)
pot_BB <- exp(fit_i$thresholds[1] * 1 + fit_i$thresholds[2] * 1 + fit_i$graph[2, 1] * 1 * 1)
Z <- pot_AA + pot_AB + pot_BA + pot_BB

# Compute probabilities
pot_AA / Z
pot_AB / Z
pot_BA / Z
pot_BB / Z


# 

# ---------- Version II: A = -1, B = 1 ---------

data_ii <- data
data_ii[data_ii==1] <- -1
data_ii[data_ii==2] <- 1

# Loglinear model
fit_ii <- EstimateIsing(data = data_ii, responses = c(-1,1), method = "ll")
fit_ii$graph
fit_ii$thresholds

## Recover conditional probabilities to get a feeling for the model:
tb <- table(data_ii[,1], data_ii[,2])
tb/sum(tb)

# Compute potentials
pot_AA <- exp(fit_ii$thresholds[1] * -1 + fit_ii$thresholds[2] * -1 + fit_ii$graph[2, 1] * -1 * -1)
pot_AB <- exp(fit_ii$thresholds[1] * -1 + fit_ii$thresholds[2] * 1 + fit_ii$graph[2, 1] * -1 * 1)
pot_BA <- exp(fit_ii$thresholds[1] * 1 + fit_ii$thresholds[2] * -1 + fit_ii$graph[2, 1] * 1 * -1)
pot_BB <- exp(fit_ii$thresholds[1] * 1 + fit_ii$thresholds[2] * 1 + fit_ii$graph[2, 1] * 1 * 1)
Z <- pot_AA + pot_AB + pot_BA + pot_BB

# Compute probabilities
pot_AA / Z
pot_AB / Z
pot_BA / Z
pot_BB / Z


