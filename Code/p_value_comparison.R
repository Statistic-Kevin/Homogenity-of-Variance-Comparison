library(ggplot2)
library(car)
library(GAD)
library(stats)

# Set parameters
nsim <- 1000
alpha <- 0.05
k1 <- k2 <- k3 <- 20

# Function to compute p-values for different tests
compute_p_values <- function(k1, k2, k3, alpha, nsim) {
  p_values_bartlett <- numeric(nsim)
  p_values_levene <- numeric(nsim)
  p_values_cochran <- numeric(nsim)
  
  for (i in 1:nsim) {
    # Generate data
    g1 <- rnorm(k1, 0, 1)
    g2 <- rnorm(k2, 0, 1)
    g3 <- rnorm(k3, 0, 2)
    x <- c(g1, g2, g3)
    group <- as.factor(c(rep(1, k1), rep(2, k2), rep(3, k3)))
    
    # Compute p-values
    p_values_levene[i] <- leveneTest(x, group)$`Pr(>F)`[1]
    p_values_bartlett[i] <- bartlett.test(x, group)$p.value
    p_values_cochran[i] <- C.test(lm(x ~ group))$p.value
  }
  
  data.frame(
    Test = rep(c("Bartlett", "Levene", "Cochran"), each = nsim),
    p_value = c(p_values_bartlett, p_values_levene, p_values_cochran)
  )
}

# Simulate data and compute p-values
data <- compute_p_values(k1, k2, k3, alpha, nsim)

# Create the plot
ggplot(data, aes(x = p_value, fill = Test)) +
  geom_histogram(binwidth = 0.05, position = "dodge", alpha = 0.7) +
  labs(title = "Comparison of Test Statistics",
       x = "p-value",
       y = "Frequency") +
  theme_minimal() +
  facet_wrap(~ Test)
