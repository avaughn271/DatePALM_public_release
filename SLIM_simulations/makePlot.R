
library(ggplot2)
library(cowplot)

MLE <- c()
logp <- c()

for (i in 1:30) {
  A <- read.csv(paste0("balancing_chroms/Results/output", i, ".txt"), sep = "\t")
  s <- A$sgrad_mle
  z <- A$sgrad_z
  p <- log10(2) + pnorm(abs(z), lower.tail = FALSE, log.p = TRUE) / log(10)
  MLE <- c(MLE, s)
  logp <- c(logp, -p)
}

# Step 1: Put into a data frame
df <- data.frame(MLE = MLE, logp = logp)

# Step 2: Make ggplot
p1 = ggplot(df, aes(x = MLE, y = logp)) +
  geom_point(color = "royalblue") +
  xlim(-0.015, 0.015) + ylim(0, 3) +
  geom_hline(yintercept = -log10(0.05/30), color = "darkgreen", linetype = "dashed") + 
  labs(
    title = "Stabilizing Selection", 
    x = expression("Estimated Selection Gradient " * hat(omega)^"MLE"),
    y = expression("-log"["10"]*"(p)")
  ) +
  theme_minimal() + theme(plot.title = element_text(hjust = 0.5)) 

 
MLE <- c()
logp <- c()

for (i in 1:30) {
  A <- read.csv(paste0("directional_chroms/Results/output", i, ".txt"), sep = "\t")
  s <- A$sgrad_mle
  z <- A$sgrad_z
  p <- log10(2) + pnorm(abs(z), lower.tail = FALSE, log.p = TRUE) / log(10)
  MLE <- c(MLE, s)
  logp <- c(logp, -p)
}

# Step 1: Put into a data frame
df <- data.frame(MLE = MLE, logp = logp)

# Step 2: Make ggplot
p2 = ggplot(df, aes(x = MLE, y = logp)) +
  geom_point(color = "royalblue") +  geom_hline(yintercept = -log10(0.05/30), color = "darkgreen", linetype = "dashed") + 
  xlim(-0.005, 0.015) + ylim(0, 700) + 
   geom_vline(xintercept = 0.012, color = "red", linetype = "dotted") + 
  labs(
    title = "Directional Selection",
    x = expression("Estimated Selection Gradient " * hat(omega)^"MLE"),
    y = expression("-log"["10"]*"(p)")
  ) +
  theme_minimal() + theme(plot.title = element_text(hjust = 0.5))











 



library(ggplot2)

MLE <- c()
logp <- c()

for (i in 1:30) {
  A <- read.csv(paste0("neutral_chroms/Results/output", i, ".txt"), sep = "\t")
  s <- A$sgrad_mle
  z <- A$sgrad_z
  p <- log10(2) + pnorm(abs(z), lower.tail = FALSE, log.p = TRUE) / log(10)
  MLE <- c(MLE, s)
  logp <- c(logp, -p)
}

# Step 1: Put into a data frame
df <- data.frame(MLE = MLE, logp = logp)

# Step 2: Make ggplot
p3 = ggplot(df, aes(x = MLE, y = logp)) +
  geom_point(color = "royalblue") +  geom_hline(yintercept = -log10(0.05/30), color = "darkgreen", linetype = "dashed") + 
  xlim(-0.005, 0.005) + ylim(0, 3) + 
  labs(
    title = "Neutrality",
    x = expression("Estimated Selection Gradient " * hat(omega)^"MLE"),
    y = expression("-log"["10"]*"(p)")
  ) +
  theme_minimal() + theme(plot.title = element_text(hjust = 0.5))






 
MLE <- c()
logp <- c()

for (i in 1:30) {
  A <- read.csv(paste0("shifting_chroms/Results/output", i, ".txt"), sep = "\t")
  s = A$delta_mle[1]
  z = A$delta_z[1]
  p <- log10(2) + pnorm(abs(z), lower.tail = FALSE, log.p = TRUE) / log(10)
  MLE <- c(MLE, s)
  logp <- c(logp, -p)
}

# Step 1: Put into a data frame
df <- data.frame(MLE = MLE, logp = logp)

# Step 2: Make ggplot
p4 = ggplot(df, aes(x = MLE, y = logp))  + 
  geom_point(color = "royalblue") +  geom_hline(yintercept = -log10(0.05/30), color = "darkgreen", linetype = "dashed") + 
  xlim(0, .2) + ylim(0, 200) + 
  labs(
    title = "Stabilizing with Changing Optimum",
    x = expression("Estimated Selection Gradient Change " * hat(Delta)^"MLE"),
    y = expression("-log"["10"]*"(p)")
  ) +
  theme_minimal() + theme(plot.title = element_text(hjust = 0.5))

plot_grid(p3, p2, p1, p4, labels = c("(a)", "(b)", "(c)", "(d)"), ncol = 2)
