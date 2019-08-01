
#====================================================================

pimatfn = function(logitp0, act.center, trap.locations, sex, logsigmam, logsigmaf)
{
  logsigma = rep(logsigmaf, nrow(act.center))
  logsigma[sex == 1] = logsigmam
  # logsigma[sex == 0] = logsigmaf
  pimat = expit(logitp0) * exp(- e2dist(act.center, trap.locations) ^ 2 / (2 * (exp(logsigma) ^ 2))) # Division is always column wise
  # pimat = 1- exp(- expit(logitp0) * exp(- e2dist(act.center, trap.locations) ^ 2 / (2 * (exp(logsigma) ^ 2))))
  return(pimat)
}

#============================
# log - likelihood function
#============================
logLfn.da = function(ncap, n0mat, n0vec, logitphi, pimat, z, J) 
{
  phi = expit(logitphi)
  yyy1 = log( 1 + phi ^ ncap)
  A1 = log(exp(yyy1) - 1)
  yyy2 = log(1 + (1 - phi) ^ (2*n0vec - ncap))
  A2 = log(exp(yyy2) - 1)
  yyy3 = log(1 + pimat ^ n0mat)
  A3 = rowSums(log(exp(yyy3) - 1), na.rm = T)
  yyy4 = log(1 + (1 - phi* (2 - phi) * pimat) ^ (J-n0mat))
  A4 = rowSums(log(exp(yyy4) - 1), na.rm = T)
  A1[A1 == -Inf] = -10^200
  A2[A2 == -Inf] = -10^200
  A3[A3 == -Inf] = -10^200
  A4[A4 == -Inf] = -10^200
  return(sum(z*(A1 + A2 + A3 + A4)))
}

#====================================================================
# compute Monte Carlo standard error of a functional (FUNC) of a Markov chain (x) using the Subsampling Bootstrap Method (SBM)
mcse = function(x, FUNC,  batchsize = floor(sqrt(length(x)))) {
  n = length(x)
  nb = n - batchsize + 1
  bval = rep(NA,nb)
  for (j in 1:nb) {
    ind = j:(batchsize + j - 1)
    bval[j] = FUNC(x[ind])
  }
  ##  var.bval = var(bval)*(nb - 1) * n * batchsize / ( (n - batchsize) * nb )    #  OLBM method
  var.bval = var(bval)*(nb - 1) * batchsize / nb
  
  return(list(se = sqrt(var.bval / n), batchsize = batchsize))
}


#====================================================================
#  utility functions
#====================================================================
Mode <- function(x) 
{
  ux <- unique(x)
  ux[which.max(tabulate(match(x, ux)))]
}

logit = function(xx) { log(xx) - log(1 - xx) }

expit = function(xx) { ifelse(xx < 0, exp(xx) /(1 + exp(xx)),  1/(1 + exp(-xx))) }

lowerQuantile = function(z) {  quantile(z, probs = 0.025)  }

upperQuantile = function(z) {  quantile(z, probs = 0.975)  }

medQuantile = function(z) {  quantile(z, probs = 0.5)  }

bias = function(x){  mean((x[-1]-x[1]), na.rm =T)}

relbiaspercentage = function(x){  100*mean(abs(x[-1]-x[1]), na.rm =T)/x[1]}

mse = function(x){  mean((x[-1]-x[1])^2, na.rm =T)}

rmse = function(x){  sqrt(mean((x[-1]-x[1])^2, na.rm =T))}
sim.stats = function(x, hpdprob = 0.95)
{
  bias =  mean((x[-1]-x[1]), na.rm =T)
  relbiaspercentage = 100*mean(abs(x[-1]-x[1]), na.rm =T)/x[1]
  rmse =  sqrt(mean((x[-1]-x[1])^2, na.rm =T))
  hpd = HPDinterval(as.mcmc(x[-1]), prob = hpdprob)
  return(c(bias, relbiaspercentage, rmse, hpd))
}

e2dist = function(x,y) # x is matrix mx2, y is matrix nx2
{
  ind = sort(rep(1:nrow(y), nrow(x))) # [rep(1, nrow(x)), rep(2, nrow(x)), ..., rep(nrow(y), nrow(x))]'
  dvec = sqrt((x[,1] - y[ind,1]) ^ 2 + (x[,2] - y[ind,2]) ^ 2) # [d(1,1), d(2,1), ...,d(nrow(x), 1), ..., d(1,nrow(y)), d(2,nrow(y)),...,d(nrow(x), nrow(y))]
  dmat = matrix(dvec, nrow = nrow(x), ncol = nrow(y), byrow = F) 
  return(dmat)
}
makedata3d = function(capture, traploc)
{
  nID = max(capture[,2])
  nLOC = dim(traploc)[1]
  nSO = dim(traploc)[2] - 3
  
  traploc = traploc[order(traploc[,1]),]
  data3d = structure(rep(0, times = nID * nLOC * nSO), .Dim = c(nID, nLOC, nSO))
  len = length(capture[, 1])
  for (i in 1:len) 
  {
    id = capture[i, 2]
    so = capture[i, 3]
    loc = capture[i, 4]
    data3d[id, loc, so] = 1
  }
  return(data3d)
}
Lvec =  function(left, right, numl, numr, sexl, sexr, trap.locations, IDfixed, nloop)
{
  
  e2dist = function(x,y)
  {
    ind = sort(rep(1:nrow(y), nrow(x))) # [rep(1, nrow(x)), rep(2, nrow(x)), ..., rep(nrow(y), nrow(x))]'
    dvec = sqrt((x[,1] - y[ind,1]) ^ 2 + (x[,2] - y[ind,2]) ^ 2) # [d(1,1), d(2,1), ...,d(nrow(x), 1), ..., d(1,nrow(y)), d(2,nrow(y)),...,d(nrow(x), nrow(y))]
    matrix(dvec, nrow = nrow(x), ncol = nrow(y), byrow = F) 
  }
  
  # function needs to spit out initial L and activity centers
  M = dim(left)[1]
  K = dim(left)[2]
  
  idkp = (IDfixed + 1):M
  
  numl = numl - IDfixed
  numr = numr - IDfixed
  
  left2d = apply(left[idkp,,], c(1,2), sum, na.rm = T) # total left captures per trap for indinvidals (IDfixed+1):M
  right2d = apply(right[idkp,,],c(1,2), sum, na.rm = T) # total right captures per trap for indinvidals (IDfixed+1):M
  
  sexl = sexl[idkp]
  sexr = sexr[idkp]
  
  temp.lm = which(sexl == 1, arr.ind = T)
  temp.lf = which(sexl == 0, arr.ind = T)
  
  temp.rm = which(sexr == 1, arr.ind = T)
  temp.rf = which(sexr == 0, arr.ind = T)
  
  temp.lr = which(sexl == sexr, arr.ind = T)
  
  trap.locations = as.matrix(trap.locations)
  
  
  sbar.left = matrix(NA,nrow = (M - IDfixed),ncol = 2)
  sbar.right = matrix(NA,nrow = (M - IDfixed),ncol = 2)
  
  #------------------------------------------------------------------------
  # Getting an initial value for activity centres sbar.left and sbar.right
  #------------------------------------------------------------------------
  for (i in 1:(M - IDfixed))
  {
    #-----------  
    # FOR LSIs
    #-----------
    if (sum(left2d[i,]) > 0) # should always be satisfied i.e LSi i got captured at least once or not
    {  
      # get the traps where ith LSI got captured
      traps.loc = matrix(trap.locations[left2d[i,] > 0,], ncol = 2, byrow = F)
      # trps<-matrix(traps,ncol = 2,byrow = FALSE)
      sbar.left[i,] = apply(traps.loc, 2, weighted.mean, w = left2d[i,][left2d[i,] > 0]) #c(mean(trps[,1]),mean(trps[,2]))
    }
    if (sum(left2d[i,]) == 0) sbar.left[i,] = trap.locations[sample(1:K,1),] # if LSi i never got captured
    
    #-----------  
    # FOR RSIs
    #-----------
    
    if (sum(right2d[i,]) > 0) # should always be satisfied i.e RSi i got captured at least once or not
    { 
      # get the traps where ith LSI got captured
      traps.loc = matrix(trap.locations[right2d[i,] > 0,], ncol = 2, byrow = F)
      # trps<-matrix(traps,ncol = 2,byrow = FALSE)
      sbar.right[i,] = apply(traps.loc, 2, weighted.mean, w = right2d[i,][right2d[i,] > 0]) #c(mean(trps[,1]),mean(trps[,2]))
    }
    if (sum(right2d[i,]) == 0) sbar.right[i,] = trap.locations[sample(1:K,1),] # if RSi i never got captured
    
  } # end of for(i in 1:(M-IDfixed))
  
  #Getting euclidean distance between each right act centre sbar.right and each leftt act centre sbar.left
  D = e2dist(sbar.right,sbar.left) # dimension (M-IDfixed) x (M-IDfixed)
  
  L = sample(1:(M - IDfixed), (M - IDfixed), replace = F)
  
  Q = sum( D[cbind(1:(M - IDfixed),L)]  ) # sum of all the distances between pair sbar.right & sbar.left
  
  for (loop in 1:nloop)
  {
    for (i in 1:(M - IDfixed))
    {
      # Checking distance sum by linking RSI i to other LSI other than L[i]
      curr.spot = L[i]
      Qtmp = c()  
      
      ind1 = ifelse(sexr[i] == 1, setdiff(temp.lm, i), setdiff(temp.lf, i))
      # ind1 = setdiff(1:(M - IDfixed), i)
      
      for (k in ind1)
      {
        L[i] = L[k];  L[k] = curr.spot # swap
        Qtmp = c(Qtmp, sum( D[cbind(1:(M - IDfixed),L)] ))
        L[k] = L[i];  L[i] = curr.spot # swap back
      } # end of for(k in ind1)
      
      if (min(Qtmp, na.rm = T) < Q )
      {
        # Make the swap
        which_L = ind1[Qtmp == min(Qtmp)][1] # ind1[which(Qtmp == min(Qtmp), arr.ind = T)][1]
        L[i] = L[which_L]
        L[which_L] = curr.spot
        
        Q = min(Qtmp, na.rm = T)
      } # end of if(min(Qtmp, na.rm = T) < Q )
      
    } # close loop over 1:(M-IDfixed) [ =length of L ]
  } # close loop over 1:50
  
  if (IDfixed > 0) L = c((1:IDfixed) ,  IDfixed + L )
  if (IDfixed == 0) L = L
  
  return(L)
}

