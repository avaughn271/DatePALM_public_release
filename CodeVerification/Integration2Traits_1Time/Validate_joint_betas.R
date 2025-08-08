
#####Checking marginal 1
integrand1_beta1 =  function(omega, beta) {
    return(dnorm( omega * beta ,  1/8, (1/8)**0.5 , log =T)
            + dnorm(  beta,  -0.1, 0.25 , log =T))
}

integrand2_beta1 =  function(omega, beta) {
    return(dnorm(  omega * beta ,  0.1, (1/6)**0.5 , log =T)
            + dnorm(  beta, -0.15, 0.5 , log =T))
}

objective1 = function(omega) {
    print(omega)
    seq1 = seq(-5, 5, length.out = 20000)
    SUM = 0
    for (i in seq1) {
            SUM = SUM + exp(integrand1_beta1(omega, i) )
    }
    return(log(SUM*10/20000))
}

objective2 = function(omega) {
    print(omega)
    seq1 = seq(-5, 5, length.out = 20000)
    SUM = 0
    for (i in seq1) {
        SUM = SUM + exp(integrand2_beta1(omega, i) )
    }
    return(log(SUM*10/20000))
}

objective3 = function(omega) {-objective1(omega) - objective2(omega)}

optimize(objective3, c(-10,10))
#####Checking marginal 2
####################################
integrand1_beta1 =  function(omega, beta) {
    return(dnorm( omega * beta ,  1/8, (1/8)**0.5 , log =T)
           + dnorm(  beta,  0.01, 0.1 , log =T))
}

integrand2_beta1 =  function(omega, beta) {
    return(dnorm(  omega * beta ,  0.1, (1/6)**0.5 , log =T)
           + dnorm(  beta, -0.14, 0.3 , log =T))
}

objective1 = function(omega) {
    print(omega)
    seq1 = seq(-5, 5, length.out = 20000)
    SUM = 0
    for (i in seq1) {
        SUM = SUM + exp(integrand1_beta1(omega, i) )
    }
    return(log(SUM*10/20000))
}

objective2 = function(omega) {
    print(omega)
    seq1 = seq(-5, 5, length.out = 20000)
    SUM = 0
    for (i in seq1) {
        SUM = SUM + exp(integrand2_beta1(omega, i) )
    }
    return(log(SUM*10/20000))
}

objective3 = function(omega) {-objective1(omega) - objective2(omega)}

optimize(objective3, c(-10,10))
############################joint optimization
####################################
integrand1_beta1 =  function(omega, beta) {
    return(dnorm( sum(omega * beta) ,  1/8, (1/8)**0.5 , log =T)
           + dnorm(  beta[1],  -0.1, 0.25 , log =T)
           + dnorm(  beta[2], 0.01, 0.1 , log =T)
           )
}

integrand2_beta1 =  function(omega, beta) {
    return(dnorm(  sum(omega * beta ),  0.1, (1/6)**0.5 , log =T)
           + dnorm(  beta[2], -0.14, 0.3 , log =T)
           + dnorm(  beta[1], -0.15, 0.5 , log =T))
}

objective1 = function(omega) {
    print(omega)
    POINTS = 580
    BOUNDS = 4
    seq1 = seq(-BOUNDS, BOUNDS, length.out = POINTS)
    seq2 = seq(-BOUNDS, BOUNDS, length.out = POINTS)
    
    SUM = matrix(0, nrow = POINTS, ncol = POINTS)
    for (i in 1:length(seq1)) {
        for (j in 1:length(seq2)) {
        SUM[i,j] = exp(integrand1_beta1(omega, c(seq1[i],seq2[j]))) 
        }
    }
    SUM2 = matrix(0, nrow = POINTS, ncol = POINTS)
    
    for (i in 1:length(seq1)) {
        for (j in 1:length(seq2)) {
            SUM2[i,j] = exp(integrand2_beta1(omega, c(seq1[i],seq2[j])) )
        }
    }
    return(log(sum(SUM2)*(2*BOUNDS)/POINTS*(2*BOUNDS)/POINTS) + log(sum(SUM)*(2*BOUNDS)/POINTS*(2*BOUNDS)/POINTS) )
    
}

objective3 = function(omega) {-objective1(omega)}

optim(c(0,0), objective3)