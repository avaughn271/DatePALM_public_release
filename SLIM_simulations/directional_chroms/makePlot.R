MLE = c()
logp = c()

for (i in 1:30) {
      A =   read.csv(paste0("Results/output", toString(i), ".txt"), sep = "\t")
      s = A$sgrad_mle
      z = A$sgrad_z
      p = log10(2) +  pnorm(abs( z), lower.tail = F , log.p = T) / log(10)
      MLE = c(MLE, s)
      logp = c(logp, -p)
}

plot(MLE, logp, xlim = c(-0.005,0.015), ylab =  expression("-log"["10"]*"(p)" ) ,
     main = "Directional Selection",
     xlab  =  expression("Estimated Selection Gradient " * hat(omega)^"MLE" * ""))





library(ggplot2)

MLE <- c()
logp <- c()

for (i in 1:30) {
  A <- read.csv(paste0("Results/output", i, ".txt"), sep = "\t")
  s <- A$sgrad_mle
  z <- A$sgrad_z
  p <- log10(2) + pnorm(abs(z), lower.tail = FALSE, log.p = TRUE) / log(10)
  MLE <- c(MLE, s)
  logp <- c(logp, -p)
}

# Step 1: Put into a data frame
df <- data.frame(MLE = MLE, logp = logp)

# Step 2: Make ggplot
ggplot(df, aes(x = MLE, y = logp)) +
  geom_point(color = "royalblue") +  geom_hline(yintercept = -log10(0.05/30), color = "darkgreen", linetype = "dashed") + 
  xlim(-0.005, 0.015) + ylim(0, 680) + 
  labs(
    title = "Directional Selection",
    x = expression("Estimated Selection Gradient " * hat(omega)^"MLE"),
    y = expression("-log"["10"]*"(p)")
  ) +
  theme_minimal() + theme(plot.title = element_text(hjust = 0.5))