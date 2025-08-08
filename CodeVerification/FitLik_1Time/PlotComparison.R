Approx = read.table("truevalues1d.txt")
width = 5
col2 = "#e69F00"
col1 = "#009E73"
FittedDistribution = scan("chr_pos.txt")

plot(Approx[[1]], exp(Approx[[2]] - max(Approx[[2]])),  xlab="s" , ylab="log-likelihood", 
     col=col1 , lwd=width, type = "l"  ,   xaxt = "n")
axis(1, at = c(seq(.012,0.04, by = 0.004)))
yvals = dnorm(Approx[[1]], mean =  FittedDistribution[1],  
              sd = sqrt(FittedDistribution[2]))
yvals = yvals / max(yvals)
lines(Approx[[1]] , yvals , col=col2 , lwd=width,lty = 2)

# Add a legend
legend("topleft", 
       legend =rev( c("True Likelihood Function", "Fitted Likelihood Function")), 
       col = rev( c(col1, col2)), 
       lty = rev( c(1,2)), lwd = 3, 
       cex = 1, inset = c(.01,.01))
