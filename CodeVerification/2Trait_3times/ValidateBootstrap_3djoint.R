mu1 = c(.125,.125,.125)
mu2 = c(0,0.16666666666,-0.2)
mu3 = c(0.5,0.25,-0.3)
mu4 = c(-1,0.5,0.5)
mu5 = c(.125,.125,.125)

sigma1 = c(0.1, 0.0, 0.0, 
          0.0, 1.5, 0.0,
          0.0, 0.0, 1.3)
sigma2 = c(0.7, 0.0, 0.4,
          0.0, 0.6, 0.3,
          0.4, 0.3, 0.8)
sigma3 = c(0.9, 0.3, 0.2,
          0.3, 0.2, 0.0,
          0.2, 0.0, 0.3)
sigma4 = c(0.4, 0.0, 0.1, 
          0.0, 0.5, 0.2,
          0.1, 0.2, 0.9)
sigma5 = c(0.1, 0.0, 0.0, 
           0.0, 1.5, 0.0, 
           0.0, 0.0, 1.3)
    
sigma1 = matrix(sigma1, ncol = 3)
sigma2 = matrix(sigma2, ncol = 3)
sigma3 = matrix(sigma3, ncol = 3)
sigma4 = matrix(sigma4, ncol = 3)
sigma5 = matrix(sigma5, ncol = 3)

multivariate_normal_pdf <- function(x, mean_vector, covariance_matrix) {
    k <- length(mean_vector)
    det_cov <- det(covariance_matrix)
    inv_cov <- solve(covariance_matrix)
    const_term <- 1 / ((2 * pi)^(k / 2) * sqrt(det_cov))
    
    mahalanobis_dist <- as.numeric(t(x - mean_vector) %*% inv_cov %*% (x - mean_vector))
    
    pdf_value <- const_term * exp(-0.5 * mahalanobis_dist)
    
    return(pdf_value)
}

#################################################
dotprodvector  = function(x, entry1, entry2) {
    return(c( (x[1]*entry1 + x[2]*entry2), 
        (x[3]*entry1 + x[4]*entry2),
        (x[5]*entry1 + x[6]*entry2)))
}

TotalFunction = function(x) {
    A = log(multivariate_normal_pdf(dotprodvector(x, -0.1, 0.01), mu1, sigma1))
    B = log(multivariate_normal_pdf(dotprodvector(x, 0.01, -0.02), mu2, sigma2))
    C = log(multivariate_normal_pdf(dotprodvector(x, -0.15, -0.14), mu3, sigma3))
    D = log(multivariate_normal_pdf(dotprodvector(x, 0.07, -0.004), mu4, sigma4))
    E = log(multivariate_normal_pdf(dotprodvector(x, 0.01, 0.03), mu5, sigma5))
    
    return( (-A-B)/2  + (-C-D-E)/3)
}

Bootstrap1 =  function(x) {
    A = log(multivariate_normal_pdf(dotprodvector(x, -0.1, 0.01), mu1, sigma1))
    B = log(multivariate_normal_pdf(dotprodvector(x, 0.01, -0.02), mu2, sigma2))
     return(-A-B )
}
Bootstrap2 =  function(x) {
    C = log(multivariate_normal_pdf(dotprodvector(x, -0.15, -0.14), mu3, sigma3))
    D = log(multivariate_normal_pdf(dotprodvector(x, 0.07, -0.004), mu4, sigma4))
    E = log(multivariate_normal_pdf(dotprodvector(x, 0.01, 0.03), mu5, sigma5))
    return(-C-D-E)
}
optim(c(0,0,0,0,0,0),TotalFunction, control = list(abstol = 0, reltol = 0,maxit = 100000))$par
optim(c(0,0,0,0,0,0),Bootstrap1, control = list(abstol = 0, reltol = 0,maxit = 100000))$par
optim(c(0,0,0,0,0,0),Bootstrap2, control = list(abstol = 0, reltol = 0,maxit = 100000))$par

###################################
dotprodvector1  = function(x, entry1) {
    return(c( x[1]*entry1  ,  x[2]*entry1  , x[3]*entry1  ))
}

TotalFunction1 = function(x) {
    A = log(multivariate_normal_pdf(dotprodvector1(x, -0.1 ), mu1, sigma1))
    B = log(multivariate_normal_pdf(dotprodvector1(x, 0.01 ), mu2, sigma2))
    C = log(multivariate_normal_pdf(dotprodvector1(x, -0.15 ), mu3, sigma3))
    D = log(multivariate_normal_pdf(dotprodvector1(x, 0.07 ), mu4, sigma4))
    E = log(multivariate_normal_pdf(dotprodvector1(x, 0.01 ), mu5, sigma5))
    
    return( (-A-B)/2  + (-C-D-E)/3)
}

Bootstrap1_1 =  function(x) {
    A = log(multivariate_normal_pdf(dotprodvector1(x, -0.1 ), mu1, sigma1))
    B = log(multivariate_normal_pdf(dotprodvector1(x, 0.01 ), mu2, sigma2))
    return(-A-B )
}
Bootstrap2_1 =  function(x) {
    C = log(multivariate_normal_pdf(dotprodvector1(x, -0.15 ), mu3, sigma3))
    D = log(multivariate_normal_pdf(dotprodvector1(x, 0.07 ), mu4, sigma4))
    E = log(multivariate_normal_pdf(dotprodvector1(x, 0.01 ), mu5, sigma5))
    return(-C-D-E)
}
optim(c(0,0,0),TotalFunction1, control = list(abstol = 0, reltol = 0,maxit = 100000))$par
optim(c(0,0,0),Bootstrap1_1, control = list(abstol = 0, reltol = 0,maxit = 100000))$par
optim(c(0,0,0),Bootstrap2_1, control = list(abstol = 0, reltol = 0,maxit = 100000))$par
####################################
 
TotalFunction2 = function(x) {
    A = log(multivariate_normal_pdf(dotprodvector1(x, 0.01), mu1, sigma1))
    B = log(multivariate_normal_pdf(dotprodvector1(x,   -0.02), mu2, sigma2))
    C = log(multivariate_normal_pdf(dotprodvector1(x,  -0.14), mu3, sigma3))
    D = log(multivariate_normal_pdf(dotprodvector1(x,   -0.004), mu4, sigma4))
    E = log(multivariate_normal_pdf(dotprodvector1(x,  0.03), mu5, sigma5))
    
    return( (-A-B)/2  + (-C-D-E)/3)
}

Bootstrap12 =  function(x) {
    A = log(multivariate_normal_pdf(dotprodvector1(x, 0.01), mu1, sigma1))
    B = log(multivariate_normal_pdf(dotprodvector1(x,   -0.02), mu2, sigma2))
    return(-A-B )
}
Bootstrap22 =  function(x) {
    C = log(multivariate_normal_pdf(dotprodvector1(x,  -0.14), mu3, sigma3))
    D = log(multivariate_normal_pdf(dotprodvector1(x,   -0.004), mu4, sigma4))
    E = log(multivariate_normal_pdf(dotprodvector1(x,  0.03), mu5, sigma5))
    return(-C-D-E)
}
optim(c(0,0,0),TotalFunction2, control = list(abstol = 0, reltol = 0,maxit = 100000))$par
optim(c(0,0,0),Bootstrap12, control = list(abstol = 0, reltol = 0,maxit = 100000))$par
optim(c(0,0,0),Bootstrap22, control = list(abstol = 0, reltol = 0,maxit = 100000))$par
