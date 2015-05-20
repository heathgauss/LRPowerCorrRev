LRPowerCorr <- function ( sampsize, nsims, p, 
                          a, b, c, d, 
                          A, B, C, D,
                          or4, or5, or6, 
                          or7, or8, or9, or10,
                          fullmodel, 
                          reducedmodel, 
                          alpha, 
                          dftest, pcx1, pcx2)
{
  #Define some things that don't need to be re-defined every iteration
  
  #The data.frame reject= '1' or '0' is stored in @ the end of each iteration
  base <<- data.frame(Reject=numeric())
  
  f = .00001
  #b0 = log(avep/(1-avep))
  or0 = c/a
  or1 = (a*d)/(b*c)
  or2 = (a*C)/(A*c)
  or3 = (A*D*b*c)/(B*C*a*d)
  
  #Cholesky's Decomposition
  corr <-chol(matrix(c(1, p, p, p, p, p, p, p, p, p,
                       p, 1, p, p, p, p, p, p, p, p,
                       p, p, 1, p, p, p, p, p, p, p,
                       p, p, p, 1, p, p, p, p, p, p,
                       p, p, p, p, 1, p, p, p, p, p,
                       p, p, p, p, p, 1, p, p, p, p,
                       p, p, p, p, p, p, 1, p, p, p,
                       p, p, p, p, p, p, p, 1, p, p,
                       p, p, p, p, p, p, p, p, 1, p,
                       p, p, p, p, p, p, p, p, p, 1),
                     nrow = 10, ncol = 10))
  
  #This function generates numbers each iteration
  Generate <- function(n){
    
    #Generate Binomial Variables
    x1=rbinom(n,1,.5)
    x2=rbinom(n,1,.5)
    
    #Normalize the binomial variables
    cx1= (x1 - pcx1)/sqrt((pcx1)*(1-pcx1))
    cx2= (x2 - pcx2)/sqrt((pcx2)*(1-pcx2))
    
    x3=x1*x2
    meanx3=mean(x3)
    stdx3=sd(x3)
    #cx3=(x3 - mean(x3)) / sd(x3)
    #cx3=x3
    cx3= (x3 - mean(x3)) / sd(x3)
    
    
    #Generate Uniform Variable
    #x3=runif(n,-3,3)
    x4=runif(n,-3,3)
    x5=runif(n,-3,3)
    x6=runif(n,-3,3)
    #cx3 = x3/sqrt(3)
    cx4 = x4/sqrt(3)
    cx5 = x5/sqrt(3)
    cx6 = x6/sqrt(3)
    
    #Generate Normal
    cx7=rnorm(n,0,1)
    cx8=rnorm(n,0,1)
    cx9=rnorm(n,0,1)
    cx10=rnorm(n,0,1)
    
    xy = data.frame(runif(n))
    
    #Pass Cx1-10 into a data.frame
    gen<-data.frame(cx1,cx2,cx3,cx4,cx5,cx6,cx7,cx8,cx9,cx10)
    
    sumx3<-data.frame(meanx3,stdx3)
    
    #Store ^^ Data.frame as a Global Matrix
    center<<-data.matrix(gen, rownames.force = NA)
    sumx3<<-data.matrix(sumx3, rownames.force = NA)
    xy<<-data.matrix(xy, rownames.force = NA)
  }
  
  #Begin looping through each iteration/simulation
  for(i in 1:nsims) {
    
    Generate(sampsize)
    
    cy = center %*% corr
    cx1 = cy[,1]
    cx2 = cy[,2]
    cx3 = cy[,3]
    cx4 = cy[,4]
    cx5 = cy[,5]
    cx6 = cy[,6]
    cx7 = cy[,7]
    cx8 = cy[,8]
    cx9 = cy[,9]
    cx10 = cy[,10]
    
    cx1 = (cx1*sqrt((pcx1)*(1-pcx1)))+pcx1
    cx2 = (cx2*sqrt((pcx2)*(1-pcx2)))+pcx2
    cx3 = (cx3*sumx3[1,2]) + sumx3[1,1] 
    cx4 = cx4*sqrt(3)
    cx5 = cx5*sqrt(3)
    cx6 = cx6*sqrt(3)
    
    multiple3 <- data.frame(cx1, cx2, cx3, cx4, cx5, cx6, cx7, cx8, cx9, cx10, xy)
    
    multiple3$logit <- with(multiple3, 
                            log(or0) + 
                              (log(or1)+f)*(cx1) +
                              (log(or2)+f)*(cx2) +
                              (log(or3)+f)*(cx3) +
                              (log(or4)+f)*(cx4) +
                              (log(or5)+f)*(cx5) +
                              (log(or6)+f)*(cx6) +
                              (log(or7)+f)*(cx7) +
                              (log(or8)+f)*(cx8) +
                              (log(or9)+f)*(cx9) +
                              (log(or10)+f)*(cx10)
    )
    
    multiple3$prob <- with(multiple3,
                           (exp(logit))/(1+exp(logit))
    )
    
    multiple3$y <- with(multiple3,
                        ifelse(xy<=prob,1,0)
    )
    
    #I may not be pulling the right number here
    fullLR <- deviance(glm(noquote(fullmodel), data=multiple3, family=binomial(logit)))
    reducLR <- deviance(glm(noquote(reducedmodel), data=multiple3, family=binomial(logit)))
    
    likelihoodratio = reducLR - fullLR
    critval = qchisq(1-alpha, dftest)
    both <- data.frame(likelihoodratio, critval)
    Reject <- with(both, ifelse(likelihoodratio >= critval, 1, 0))
    
    base <<- rbind(base, data.frame(Reject))
    
  }
  
  x = table(base)
  freq = cbind(x['1'])
  power <- freq/nsims
  se <- sqrt((power*(1-power))/nsims)
  LCL <- power - (1.96*se)
  UCL <- power + (1.96*se)
  
  cat("\nx: ", x)
  cat("\nfreq: ", freq)
  cat("\nPower: ", power)
  cat("\nse: ", se)
  cat("\nLCL: ", LCL)
  cat("\nUCL: ", UCL)
  
  cat("\n\nSample Size =", sampsize, "; Simulations =", nsims, "; Rho =", p, "; P(Y=1) =", or0,
      "\nOR1 =", or1, "; OR2 =", or2, "; OR3 =", or3, "; OR4 =", or4, "; OR5 =", or5, "; OR6 =", or6,
      "\nOR7 =", or7, "; OR8 =", or8, "; OR9 =", or9, "; OR10 =", or10,
      "\nFull Model: ", fullmodel,
      "\nReduced Model: ", reducedmodel,
      "\nPower:", round(power*100,2), "%      LCL:", round(LCL*100,2), "%      UCL:", round(UCL*100,2), "%")
  
}
