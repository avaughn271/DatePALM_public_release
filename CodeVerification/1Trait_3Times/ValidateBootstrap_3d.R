mu1 =  c(.125,.125,.125)
mu2 = c(0,0.16666666666,-0.2)
mu3 = c(0.5,0.25,-0.3)
sigma1 = c(0.1, 0.0, 0.0, 
          0.0, 1.5, 0.0,
          0.0, 0.0, 1.3)
sigma2 = c(0.7, 0.0, 0.4,
          0.0, 0.6, 0.3,
          0.4, 0.3, 0.8)
sigma3 = c(0.9, 0.3, 0.2,
          0.3, 0.2, 0.0,
          0.2, 0.0, 0.3)
sigma1 = matrix(sigma1, ncol = 3)
sigma2 = matrix(sigma2, ncol = 3)
sigma3 = matrix(sigma3, ncol = 3)

multivariate_normal_pdf <- function(x, mean_vector, covariance_matrix) {
    k <- length(mean_vector)  # Dimensionality of the random variable
    det_cov <- det(covariance_matrix)
    inv_cov <- solve(covariance_matrix)
    const_term <- 1 / ((2 * pi)^(k / 2) * sqrt(det_cov))
    
    mahalanobis_dist <- as.numeric(t(x - mean_vector) %*% inv_cov %*% (x - mean_vector))
    
    pdf_value <- const_term * exp(-0.5 * mahalanobis_dist)
    
    return(pdf_value)
}

TotalFunction = function(x) {
    A = log((multivariate_normal_pdf(x * -0.001, mu1, sigma1)))
    B = log((multivariate_normal_pdf(x *-0.01, mu2, sigma2)))
    C = log((multivariate_normal_pdf(x *0.02, mu3, sigma3)))
    return(-(A+B/2+C/2  ))
}
Bootstrap1 =  function(x) {
    A = log((multivariate_normal_pdf(x * -0.001, mu1, sigma1)))
     return(-A )
}
Bootstrap2 =  function(x) {
    B = log((multivariate_normal_pdf(x *-0.01, mu2, sigma2)))
    C = log((multivariate_normal_pdf(x *0.02, mu3, sigma3)))
    return(-(B/2+C/2  ))
}
#all of them are 0.0
optim(c(0,0,0),TotalFunction, control = list(abstol = 0, reltol = 0,maxit = 100000))$par
optim(c(0,0,0),Bootstrap1, control = list(abstol = 0, reltol = 0,maxit = 100000))$par
optim(c(0,0,0),Bootstrap2, control = list(abstol = 0, reltol = 0,maxit = 100000))$par
