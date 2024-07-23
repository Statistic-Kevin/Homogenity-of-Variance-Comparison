library(car)   # For Levene's Test
library(stats) # For Bartlett Test
library(GAD)   # For Cochran's C Test

# Function to compute empirical errors for Normal Distribution
empirical_error_norm <- function(k1, k2, k3, alpha, nsim = 10000) {
  mu <- 10
  sigma <- 26
  npbart <- numeric(nsim)
  nplev <- numeric(nsim)
  npcoch <- numeric(nsim)
  
  for (i in 1:nsim) {
    g1 <- rnorm(k1, mu, sigma)
    g2 <- rnorm(k2, mu, sigma)
    g3 <- rnorm(k3, mu, sigma)
    x <- c(g1, g2, g3)
    group <- as.factor(c(rep(1, k1), rep(2, k2), rep(3, k3)))
    
    nplev[i] <- leveneTest(x, group)$`Pr(>F)`[1] < alpha
    npbart[i] <- bartlett.test(x, group)$p.value < alpha
    npcoch[i] <- C.test(lm(x ~ group))$p.value < alpha
  }
  
  list(
    bartlett_error = mean(npbart),
    levene_emp_error = mean(nplev),
    cochran_emp_error = mean(npcoch)
  )
}

# Function to compute power for Normal Distribution with specific variance ratios
power_norm <- function(k1, k2, k3, alpha, var_ratios, nsim = 10000) {
  mu <- 10
  sigma <- 26
  npbart <- numeric(nsim)
  nplev <- numeric(nsim)
  npcoch <- numeric(nsim)
  
  for (i in 1:nsim) {
    sigma2 <- var_ratios[2] * sigma
    sigma3 <- var_ratios[3] * sigma
    g1 <- rnorm(k1, mu, sigma)
    g2 <- rnorm(k2, mu, sigma2)
    g3 <- rnorm(k3, mu, sigma3)
    x <- c(g1, g2, g3)
    group <- as.factor(c(rep(1, k1), rep(2, k2), rep(3, k3)))
    
    nplev[i] <- leveneTest(x, group)$`Pr(>F)`[1] < alpha
    npbart[i] <- bartlett.test(x, group)$p.value < alpha
    npcoch[i] <- C.test(lm(x ~ group))$p.value < alpha
  }
  
  list(
    bartlett_power = mean(npbart),
    levene_power = mean(nplev),
    cochran_power = mean(npcoch)
  )
}

# Function to compute empirical errors for Chi-Squared Distribution
empirical_error_chisq <- function(k1, k2, k3, alpha, nsim = 10000) {
  npbart <- numeric(nsim)
  nplev <- numeric(nsim)
  npcoch <- numeric(nsim)
  
  for (i in 1:nsim) {
    g1 <- rchisq(k1, df = 2)
    g2 <- rchisq(k2, df = 2)
    g3 <- rchisq(k3, df = 2)
    x <- c(g1, g2, g3)
    group <- as.factor(c(rep(1, k1), rep(2, k2), rep(3, k3)))
    
    nplev[i] <- leveneTest(x, group)$`Pr(>F)`[1] < alpha
    npbart[i] <- bartlett.test(x, group)$p.value < alpha
    npcoch[i] <- C.test(lm(x ~ group))$p.value < alpha
  }
  
  list(
    bartlett_error = mean(npbart),
    levene_emp_error = mean(nplev),
    cochran_emp_error = mean(npcoch)
  )
}

# Function to compute power for Chi-Squared Distribution with specific degrees of freedom ratios
power_chisq <- function(k1, k2, k3, alpha, df_ratios, nsim = 10000) {
  npbart <- numeric(nsim)
  nplev <- numeric(nsim)
  npcoch <- numeric(nsim)
  
  for (i in 1:nsim) {
    df1 <- df_ratios[1]
    df2 <- df_ratios[2]
    df3 <- df_ratios[3]
    g1 <- rchisq(k1, df = df1)
    g2 <- rchisq(k2, df = df2)
    g3 <- rchisq(k3, df = df3)
    x <- c(g1, g2, g3)
    group <- as.factor(c(rep(1, k1), rep(2, k2), rep(3, k3)))
    
    nplev[i] <- leveneTest(x, group)$`Pr(>F)`[1] < alpha
    npbart[i] <- bartlett.test(x, group)$p.value < alpha
    npcoch[i] <- C.test(lm(x ~ group))$p.value < alpha
  }
  
  list(
    bartlett_power = mean(npbart),
    levene_power = mean(nplev),
    cochran_power = mean(npcoch)
  )
}



# Empirical Errors for Normal Distribution
alpha_values <- c(0.01, 0.05)
k_values <- list(c(15, 15, 15), c(30, 30, 30), c(45, 45, 45), c(15, 30, 45))
for (alpha in alpha_values) {
  for (k in k_values) {
    cat("Empirical Error - Normal, alpha =", alpha, ", k1 =", k[1], ", k2 =", k[2], ", k3 =", k[3], "\n")
    print(empirical_error_norm(k1 = k[1], k2 = k[2], k3 = k[3], alpha = alpha))
  }
}

# Power Evaluations for Normal Distribution
var_ratios_list <- list(c(1, 1, 2), c(1, 2, 4))
for (alpha in alpha_values) {
  for (k in k_values) {
    for (ratios in var_ratios_list) {
      cat("Power - Normal, alpha =", alpha, ", k1 =", k[1], ", k2 =", k[2], ", k3 =", k[3], ", Var Ratios =", paste(ratios, collapse = ":") , "\n")
      print(power_norm(k1 = k[1], k2 = k[2], k3 = k[3], alpha = alpha, var_ratios = ratios))
    }
  }
}

# Empirical Errors for Chi-Squared Distribution
alpha_values_chisq <- c(0.01, 0.05)
k_values_chisq <- list(c(15, 15, 15), c(30, 30, 30), c(45, 45, 45), c(15, 30, 45))
for (alpha in alpha_values_chisq) {
  for (k in k_values_chisq) {
    cat("Empirical Error - Chi-Squared, alpha =", alpha, ", k1 =", k[1], ", k2 =", k[2], ", k3 =", k[3], "\n")
    print(empirical_error_chisq(k1 = k[1], k2 = k[2], k3 = k[3], alpha = alpha))
  }
}

# Power Evaluations for Chi-Squared Distribution
df_ratios_list <- list(c(2, 2, 2), c(2, 4, 16))
for (alpha in alpha_values_chisq) {
  for (k in k_values_chisq) {
    for (ratios in df_ratios_list) {
      cat("Power - Chi-Squared, alpha =", alpha, ", k1 =", k[1], ", k2 =", k[2], ", k3 =", k[3], ", DF Ratios =", paste(ratios, collapse = ":") , "\n")
      print(power_chisq(k1 = k[1], k2 = k[2], k3 = k[3], alpha = alpha, df_ratios = ratios))
    }
  }
}
