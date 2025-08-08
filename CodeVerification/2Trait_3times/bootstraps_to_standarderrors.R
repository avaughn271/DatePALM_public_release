JOINT1 = c( -1.3157903  ,  -2.1929827 ,  -0.2631604  ,-0.6578985  ,-9.4298201 ,9.8684207)
JOINT2 = c(  -13.041135,  6.079011  ,  6.287562  , 9.964359  ,  -8.428696 , -4.653323)
JOINTBOTH = c(  -2.3349665  ,  2.2556910 ,  0.7916436   , 0.1866682, -3.8632317 ,  1.8060994)

for (i in 1:6) {
    vec = c(JOINT1[i], JOINT2[i], JOINTBOTH[i], JOINTBOTH[i])
   print(var(rep(vec, each = 10000))**0.5)
}

MARGBOTH = c(-2.1933301, -0.9128879 , 1.8622048 , -0.6200947, -0.9161500 , 3.0230119)
MARG1 = c( -1.2115093, -0.3423013, -2.0862898,   7.880026 , -6.579118 , 15.164073)
MARG2 = c( -5.722803, -2.173564 , 1.517481,   -0.9703177 ,  -0.9114315, 2.6964915)

for (i in 1:6) {
    vec = c(MARG1[i], MARG2[i], MARGBOTH[i], MARGBOTH[i])
    print(var(rep(vec, each = 10000))**0.5)
}

#######################################R value
R1 = JOINT1 - MARG1
R2 = JOINT2 - MARG2
Rboth = JOINTBOTH - MARGBOTH
for (i in 1:6) {
    vec = c(R1[i],  R2[i],Rboth[i], Rboth[i])
    print(var(rep(vec, each = 10000))**0.5)
}

##################################deltavalue
Delta1 = JOINT1[-1] - JOINT1[c(1,2,3,4,5)]
Delta2 = JOINT2[-1] - JOINT2[c(1,2,3,4,5)]
DeltaBoth = JOINTBOTH[-1] - JOINTBOTH[c(1,2,3,4,5)]
for (i in 1:5) {
    vec = c(Delta1[i],  Delta2[i],DeltaBoth[i], DeltaBoth[i])
    print(var(rep(vec, each = 10000))**0.5)
}
#################delta_r
RDelta1 = R1[-1] - R1[c(1,2,3,4,5)]
RDelta2 = R2[-1] - R2[c(1,2,3,4,5)]
RDeltaBoth = Rboth[-1] - Rboth[c(1,2,3,4,5)]
for (i in 1:5) {
    vec = c(RDelta1[i],  RDelta2[i],RDeltaBoth[i], RDeltaBoth[i])
    print(var(rep(vec, each = 10000))**0.5)
}