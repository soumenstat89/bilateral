
PISCR.fn = function(
  capture_1, capture_2, gender, activity, traploc, 
  ndraws = 20, burnin = 10, M = 400, scale = 10^4,
  nloopL = 50, n.update = 20, batchsize = 1000, mindelta = 0.01, 
  sigmaOfProposal.logitphi = 0.05, sigmaOfProposal.logitp0 = 0.08,
  sigmaOfProposal.logsigmam = 0.02, sigmaOfProposal.logsigmaf = 0.02,
  sigmaOfProposal.L = 4, sigmaOfProposal.s = rep(3, 400), R = 5)
  {

  library(mvtnorm)
  library(MCMCpack)
  library(abind)
  library(ggplot2)
  library(raster)
  library(maptools)
  library(coda)
  library(spatstat) 
  library(gtools)
options(digits = 8)
source("utility.functions.R")

################################################

area <- activity@polygons[[1]]@Polygons[[1]]@area/1000000
activity1<-as.owin(activity)

trap = traploc[,'LOC_ID']
traploc.xx = traploc[,'X_COORD']
traploc.yy = traploc[,'Y_COORD']
K = dim(traploc)[1]
J = dim(traploc)[2] - 3 
dimnames(traploc)[[2]][4:dim(traploc)[2]] = c(1:J)
# convert units of locations (from meters to km) and adjust locations to have  x- and y- minima at zero

traploc.x = traploc.xx/scale
traploc.y = traploc.yy/scale
trap.locations=cbind(traploc.x, traploc.y)

left.obs = makedata3d(capture_1,traploc)
numl = dim(left.obs)[1] 
left = abind(left.obs, array(0, dim = c( M - numl, K, J)), along = 1)

right.obs = makedata3d(capture_2,traploc)
numr = dim(right.obs)[1] 
right = abind(right.obs, array(0, dim = c( M - numr, K, J)), along = 1)

active_trapso=as.matrix(traploc[,4:dim(traploc)[2]])
Jvec = apply(active_trapso, 1, sum) 

msk = active_trapso
msk3d = array(NA, c(M, K, J))
for (i in 1:M) msk3d[i, 1:K, 1:J] = msk[1:K, 1:J]
left = left * msk3d
right = right * msk3d


#==========================================================
# Initializing the parameters and processing the variables
logitphi = logit(0.8); logitp0 = logit(0.1); logsigmam = log(0.2); logsigmaf = log(0.15); psi = 0.5; theta = 0.5

## HANDLING THE SEX INFORMATION ##

sexl.obs = gender[1:numl, 2]  # length numl
sexr.obs = gender[1:numr, 4]  # length numr
sexl = sexr = rep(NA, M)

sexl[c(1:length(sexl.obs))[sexl.obs == 'Male']] = 1
sexl[c(1:length(sexl.obs))[sexl.obs == 'Female']] = 0

sexr[c(1:length(sexr.obs))[sexr.obs == 'Male']] = 1
sexr[c(1:length(sexr.obs))[sexr.obs == 'Female']] = 0

l.guys = levels(gender[,1])
IDfixed = suppressWarnings(sum(!is.na(as.numeric(l.guys))))

if(IDfixed == 0) known = 'none'
if(IDfixed != numl | IDfixed != numr) known = 'some'
if(IDfixed == numl & IDfixed == numr) known = 'ALL'

# Also re-label data sets if needed so that the left side always has more individuals
if (known!="ALL" & numl < numr){ 
  a = left; left = right; right = a
  b = numl; numl = numr; numr = b
  cc = sexl; sexl = sexr; sexr = cc
}


##################################################

numlm = sum(sexl == 1, na.rm = T); numlf = sum(sexl == 0, na.rm = T)
numrm = sum(sexr == 1, na.rm = T); numrf = sum(sexr == 0, na.rm = T)

sexl1 = c(); sexr1 = c()
if(numlm < numrm) sexl1 = rep(1, numrm - numlm)
if(numrm < numlm) sexr1 = rep(1, numlm - numrm) 
if(numlf < numrf) sexl1 = c(sexl1, rep(0, numrf - numlf))
if(numrf < numlf) sexr1 = c(sexr1, rep(0, numlf - numrf)) 


missing.sex.guys = is.na(sexl) 

nnn = max(numlm, numrm) + max(numlf, numrf) 

if(numl < nnn){ sexl[(numl+1) : nnn] = sample(sexl1, nnn - numl, replace = F)} 
if(numr < nnn){ sexr[(numr+1) : nnn] = sample(sexr1, nnn - numr, replace = F)} 

sexl[(nnn+1) : M] = sexr[(nnn+1) : M] = rbinom(M - nnn, 1, theta) 
if(sum(is.na(sexl)) > 0) {sexl[is.na(sexl)] = sexr[is.na(sexr)] = rbinom(sum(is.na(sexl)), 1, theta)}
# Latent variable vector for gender
sex = sexl # M x 1
#######################################################
if (known != 'ALL'){ L = Lvec(left, right, numl, numr, sexl, sexr, trap.locations, IDfixed, nloopL)}
s.left <- cbind((runifpoint(n = M , win = activity1))[[3]]/scale, 
                (runifpoint(n = M , win = activity1))[[4]]/scale)
# =============================

ncapl = rowSums(left) 
ncapr = rowSums(right)

if (known != 'ALL') {ncapr.star = ncapr[order(L)];  right.star = right[order(L),,]}
if (known == 'ALL') {ncapr.star = ncapr;            right.star = right}
ncap = ncapl + ncapr.star
n0mat = apply(left + right.star, c(1,2), function(a){sum(a > 0)})
n0vec = rowSums(n0mat) 
zero.guys = (ncapl + ncapr.star) == 0 
numcap = sum((ncapl + ncapr.star) > 0) 
z = c(rep(1, numcap), rbinom(M - numcap, 1, psi))
pimat = pimatfn(logitp0, s.left, trap.locations, sex, logsigmam, logsigmaf)
llik = logLfn.da(ncap, n0mat, n0vec, logitphi, pimat, z, J)

# .... assign values for adaptive Metropolis - Hastings sampling

batch = 0
delta = mindelta
continueMCMC = TRUE
draw = drawphi = drawp0 = drawsigmam = drawsigmaf = drawL = drawsex = drawact = 0
naccept.logitphi = naccept.logitp0 = naccept.logsigmam = naccept.logsigmaf = naccept.L = 0
naccept.s = rep(0,M)

ts = format(Sys.time(), "%d%m%y_%H%M%S")
folderName = paste("partiaID.da_", ts, sep = "")
dir.create(path = folderName)

start.time = Sys.time()
cat('Begin MCMC sampling:', '\n', '\n')

cat(c('N', 'psi', 'N.Male', 'theta', 'phi', 'p0', 'sigmam', 'sigmaf'), sep = ',', file = paste(folderName, '/markovchain.txt', sep = ""), append = TRUE)
cat('\n', file = paste(folderName, '/markovchain.txt', sep = ""), append = TRUE)
cat(c(paste('sex', 1:M, sep = '')), sep = ',', file = paste(folderName, '/markovchain.sex.txt', sep = ''), append = TRUE)
cat('\n', file = paste(folderName, '/markovchain.sex.txt', sep = ''), append = TRUE)
cat(c(paste('z', 1:M, sep = '')), sep = ',', file = paste(folderName, '/markovchain.z.txt', sep = ''), append = TRUE)
cat('\n', file = paste(folderName, '/markovchain.z.txt', sep = ''), append = TRUE)
cat(c(paste('sx', 1:M, sep = '')), sep = ',', file = paste(folderName, '/markovchain.sx.txt', sep = ''), append = TRUE)
cat('\n', file = paste(folderName, '/markovchain.sx.txt', sep = ''), append = TRUE)
cat(c(paste('sy', 1:M, sep = '')), sep = ',', file = paste(folderName, '/markovchain.sy.txt', sep = ''), append = TRUE)
cat('\n', file = paste(folderName, '/markovchain.sy.txt', sep = ''), append = TRUE)
cat(c(paste('loglik', 1:M, sep = '')), sep = ',', file = paste(folderName, '/markovchain.loglikelihood.txt', sep = ''), append = TRUE)
cat('\n', file = paste(folderName, '/markovchain.loglikelihood.txt', sep = ''), append = TRUE)
if (known != 'ALL') 
{
  cat(c(paste("L", 1:M, sep = '')), sep = ',', file = paste(folderName, '/markovchainL.txt', sep = ""), append = TRUE)
  cat('\n', file = paste(folderName, '/markovchainL.txt', sep = ""), append = TRUE)
}

while (continueMCMC) {
  
  draw = draw + 1
  drawinterval = 500
  if (draw == round(draw/drawinterval)*drawinterval)  cat('..... drawing sample #', draw, '\n')
  
  # update the increment/decrement for adaptive Metropolis - Hastings samplers
  if (floor(draw/batchsize)*batchsize == draw){
    batch = batch + 1
    if (1/sqrt(batch) < mindelta)  delta = 1/sqrt(batch)
  }
  
  # update phi
  drawphi = drawphi +1
  loglik.curr.phi = llik
  logitphi.cand = rnorm(1, logitphi, sigmaOfProposal.logitphi) 
  loglik.cand.phi = logLfn.da(ncap, n0mat, n0vec, logitphi.cand, pimat, z, J)
  lognum = loglik.cand.phi + log(exp(logitphi.cand)/((1+exp(logitphi.cand))^2))
  logden = loglik.curr.phi + log(exp(logitphi)/((1+exp(logitphi))^2)) 
  if (logden == -Inf){
    logitphi = logitphi.cand
    loglik.curr.phi = loglik.cand.phi
    naccept.logitphi = naccept.logitphi + 1
  }  
  if (logden != -Inf){
    logR = lognum - logden
    if (runif(1,0,1) <= exp(logR)){
      logitphi = logitphi.cand
      loglik.curr.phi = loglik.cand.phi
      naccept.logitphi = naccept.logitphi + 1
    }
  }
  if (floor(draw/batchsize)*batchsize == draw){
    SigmaDiff = ifelse(naccept.logitphi > 0.44*batchsize, exp(2*delta), exp(-2*delta))
    if(draw <= burnin){ sigmaOfProposal.logitphi = sigmaOfProposal.logitphi * SigmaDiff}
    cat("proposal sd of logitphi = ", sigmaOfProposal.logitphi, ' ')
    cat("naccept.logitphi = ", naccept.logitphi, '\n')
    naccept.logitphi = 0   
  }
  
  # update p0
  drawp0 = drawp0 +1
  loglik.curr.p0 = loglik.curr.phi
  logitp0.cand = rnorm(1, logitp0, sigmaOfProposal.logitp0) 
  pimat.cand = pimatfn(logitp0.cand, s.left, trap.locations, sex, logsigmam, logsigmaf)
  loglik.cand.p0 = logLfn.da(ncapl + ncapr.star, n0mat, n0vec, logitphi, pimat.cand, z, J)
  lognum = loglik.cand.p0 + log(exp(logitp0.cand)/((1+exp(logitp0.cand))^2))
  logden = loglik.curr.p0 + log(exp(logitp0)/((1+exp(logitp0))^2))
  if (logden == -Inf){
    logitp0 = logitp0.cand
    pimat = pimat.cand
    loglik.curr.p0 = loglik.cand.p0
    naccept.logitp0 = naccept.logitp0 + 1
  }
  if (logden != -Inf){
    logR = lognum - logden
    if (runif(1,0,1) <= exp(logR)){
      logitp0 = logitp0.cand
      pimat = pimat.cand
      loglik.curr.p0 = loglik.cand.p0
      naccept.logitp0 = naccept.logitp0 + 1
    }
  }
  if (floor(draw/batchsize)*batchsize == draw){
    SigmaDiff = ifelse(naccept.logitp0 > 0.44*batchsize, exp(2*delta), exp(-2*delta))
    if(draw <= burnin){ sigmaOfProposal.logitp0 = sigmaOfProposal.logitp0 * SigmaDiff}
    cat("proposal sd of logitp0 = ", sigmaOfProposal.logitp0, ' ')
    cat("naccept.logitp0 = ", naccept.logitp0, '\n')
    naccept.logitp0 = 0   
  }
  
  # update sigmam
  drawsigmam = drawsigmam +1
  loglik.curr.sigmam = loglik.curr.p0
  logsigmam.cand = rnorm(1, logsigmam, sigmaOfProposal.logsigmam)
  pimat.cand = pimatfn(logitp0, s.left, trap.locations, sex, logsigmam.cand, logsigmaf)
  loglik.cand.sigmam = logLfn.da(ncap, n0mat, n0vec, logitphi, pimat.cand, z, J)
  lognum = loglik.cand.sigmam + log(exp(logsigmam.cand)/R)
  logden = loglik.curr.sigmam + log(exp(logsigmam)/R) 
  if (logden == -Inf){
    logsigmam = logsigmam.cand
    pimat = pimat.cand
    loglik.curr.sigmam = loglik.cand.sigmam
    naccept.logsigmam = naccept.logsigmam + 1
  }
  if (logden != -Inf){
    logR = lognum - logden
    if (runif(1,0,1) <= exp(logR)){
      logsigmam = logsigmam.cand
      pimat = pimat.cand
      loglik.curr.sigmam = loglik.cand.sigmam
      naccept.logsigmam = naccept.logsigmam + 1
    }
  }
  if (floor(draw/batchsize)*batchsize == draw){
    SigmaDiff = ifelse(naccept.logsigmam > 0.44*batchsize, exp(2*delta), exp(-2*delta))
    if(draw <= burnin){ sigmaOfProposal.logsigmam = sigmaOfProposal.logsigmam * SigmaDiff}
    cat("proposal sd of logsigmam = ", sigmaOfProposal.logsigmam, ' ')
    cat("naccept.logsigmam = ", naccept.logsigmam, '\n')
    naccept.logsigmam = 0  
  }
  
  # update sigmaf
  drawsigmaf = drawsigmaf +1
  loglik.curr.sigmaf = loglik.curr.sigmam
  logsigmaf.cand = rnorm(1, logsigmaf, sigmaOfProposal.logsigmaf)
  pimat.cand = pimatfn(logitp0, s.left, trap.locations, sex, logsigmam, logsigmaf.cand)
  loglik.cand.sigmaf = logLfn.da(ncap, n0mat, n0vec, logitphi, pimat.cand, z, J)
  lognum = loglik.cand.sigmaf + log(exp(logsigmaf.cand)/R)
  logden = loglik.curr.sigmaf + log(exp(logsigmaf)/R)
  if (logden == -Inf){
    logsigmaf = logsigmaf.cand
    pimat = pimat.cand
    loglik.curr.sigmaf = loglik.cand.sigmaf
    naccept.logsigmaf = naccept.logsigmaf + 1
  }
  if (logden != -Inf){
    logR = lognum - logden
    if (runif(1,0,1) <= exp(logR)){
      logsigmaf = logsigmaf.cand
      pimat = pimat.cand
      loglik.curr.sigmaf = loglik.cand.sigmaf
      naccept.logsigmaf = naccept.logsigmaf + 1
    }
  }
  if (floor(draw/batchsize)*batchsize == draw){
    SigmaDiff = ifelse(naccept.logsigmaf > 0.44*batchsize, exp(2*delta), exp(-2*delta))
    if(draw <= burnin){ sigmaOfProposal.logsigmaf = sigmaOfProposal.logsigmaf * SigmaDiff}
    cat("proposal sd of logsigmaf = ", sigmaOfProposal.logsigmaf, ' ')
    cat("naccept.logsigmaf = ", naccept.logsigmaf, '\n')
    naccept.logsigmaf = 0   
  }
  
  #  update L
  if (known != "ALL"){
    loglik.curr.L = loglik.curr.sigmaf
    drawL = drawL +1
    indx1 = (IDfixed + 1):M
    if( IDfixed < numr) rightset1  = order(L)[z == 1 & order(L) >= (IDfixed+1) & order(L) <= numr]
    if( IDfixed >= numr) rightset1  = c()
    rightset2 = order(L)[z == 1 & order(L) > numr] # real right individuals who went uncaptured
    if(length(rightset2) >1) rightset2 = sample(rightset2, min(length(rightset2), n.update), replace = F)
    rightset3 = c()
    if(IDfixed < numl){
      leftset3 = c(1:M)[z == 1 & c(1:M)>IDfixed & c(1:M) <= numl]
      if(length(leftset3) > 0) for(ii in leftset3) rightset3 = c(rightset3, c(1:M)[L == ii]) 
      # right individuals to whom left captured individuals are linked with
    }
    rightset = sort(unique(c(rightset1, rightset2, rightset3)))
    for (r.guy1 in rightset){
      l.swap.out = L[r.guy1]
      possible.L = c(1:M)[z == 1 & c(1:M) > IDfixed & sex == sex[l.swap.out]]
      dv =  sqrt( (s.left[l.swap.out,1] - s.left[,1]) ^ 2 + (s.left[l.swap.out,2] - s.left[,2]) ^ 2 )
      wt.possible.L = exp( - (dv ^ 2) / sigmaOfProposal.L ^ 2)[z == 1 & c(1:M) > IDfixed  & sex == sex[l.swap.out]]
      if (length(possible.L) > 1)  l.swap.in =  sample( possible.L, 1, replace = F, prob = wt.possible.L)
      if (length(possible.L) == 1) next 
      if (length(possible.L) == 0) next # this case will never happen since l.swap.out is present there at the centre of the circle dv < 5
      if (l.swap.in == l.swap.out) next # saves computation time in for loop
      jump.prob.L =  wt.possible.L[which(possible.L == l.swap.in, arr.ind = T)] / sum(wt.possible.L) #  q(state.curr, state.cand)
      dv.back = sqrt( (s.left[l.swap.in,1] - s.left[,1]) ^ 2 + (s.left[l.swap.in,2] - s.left[,2]) ^ 2  )
      wt.possible.back.L = exp( - (dv.back ^ 2) / sigmaOfProposal.L ^ 2)[z == 1 & c(1:M) > IDfixed & sex == sex[l.swap.in]]
      # Note that, sex[l.swap.out] = sex[l.swap.in], and z[l.swap.out] may be != z[l.swap.in]
      jump.prob.back.L =  wt.possible.back.L[which(possible.L == l.swap.out, arr.ind = T)] / sum(wt.possible.back.L) #  q(state.cand, state.curr)
      ##  Which right encounter history is currently associated with left guy s.swap.in?
      r.guy2 = c(1:M)[L == l.swap.in] # which(L == l.swap.in, arr.ind = T) 
      # if (l.swap.in <=IDfixed) browser()
      L.cand = L
      L.cand[r.guy1] = l.swap.in
      L.cand[r.guy2] = l.swap.out
      right.star.cand = right[order(L.cand),,]
      ncapr.star.cand = ncapr[order(L.cand)]
      ncap.cand = ncapl + ncapr.star.cand
      n0mat.cand = apply(left + right.star.cand, c(1,2), function(a){sum(a > 0)})
      n0vec.cand = rowSums(n0mat.cand) # apply(n0mat.cand, 1, sum)
      loglik.cand.L = logLfn.da(ncap.cand, n0mat.cand, n0vec.cand, logitphi, pimat, z, J)
      lognum = loglik.cand.L + log(jump.prob.back.L)
      logden = loglik.curr.L + log(jump.prob.L)
      if (logden == -Inf){
        L = L.cand 
        loglik.curr.L = loglik.cand.L
        right.star = right.star.cand
        ncapr.star = ncapr.star.cand
        ncap = ncap.cand
        n0mat = n0mat.cand
        n0vec = n0vec.cand
        naccept.L = naccept.L + 1
      }
      if (logden != -Inf){
        logR = lognum - logden
        if (runif(1,0,1) <= exp(logR)){
          L = L.cand 
          loglik.curr.L = loglik.cand.L
          right.star = right.star.cand
          ncapr.star = ncapr.star.cand
          ncap = ncap.cand
          n0mat = n0mat.cand
          n0vec = n0vec.cand
          naccept.L = naccept.L + 1
        } 
      }
    } # end of for (r.guy1 in rightset)
    if (floor(draw/batchsize)*batchsize == draw){
      SigmaDiff = ifelse(naccept.L > 0.44*batchsize, exp(2*delta), exp(-2*delta))
      if(draw <= burnin){ sigmaOfProposal.L= sigmaOfProposal.L * SigmaDiff}
      cat(paste("proposal sd of L = ", sep = ""), sigmaOfProposal.L, ' ')
      cat(paste("naccept.L = ", sep = ""), naccept.L, '\n')
      naccept.L = 0   
    } 
  } # end of if (known!='ALL')
  
    # update z
    zero.guys = (ncapl + ncapr.star) == 0 
    ncprob = (1 - expit(logitphi) * (2 - expit(logitphi)) * pimat) ^ J 
    prob0 = exp(rowSums(log(ncprob))) 
    fc = prob0*psi / (prob0*psi + 1 - psi)
    z[zero.guys] = rbinom(sum(zero.guys), 1, fc[zero.guys])
    z[!zero.guys] = 1
    
    # update psi
    psi = rbeta(1, 1 + sum(z), 1 +  (M - sum(z) ) ) 
     
    # update sex
      drawsex = drawsex +1
      D1 = e2dist(s.left, trap.locations)
      pimat.m = expit(logitp0) * exp(- D1 * D1 / (2 * (exp(logsigmam) ^ 2))) # M x K
      Am1 = rowSums(log(pimat.m ^ n0mat), na.rm = T) # M x 1
      Am2 = rowSums(log((1 - expit(logitphi) * (2 - expit(logitphi)) * pimat.m) ^ (J-n0mat)), na.rm = T) # M x 1
      Am1[Am1 == -Inf] = -10^200; Am2[Am2 == -Inf] = -10^200; # .Machine$double.xmax
      AAm = exp(Am1 + Am2) # M x 1
      pimat.f = expit(logitp0) * exp(- D1 * D1 / (2 * (exp(logsigmaf) ^ 2))) # M x K
      Af1 = rowSums(log(pimat.f ^ n0mat), na.rm = T) # M x 1
      Af2 = rowSums(log((1 - expit(logitphi) * (2 - expit(logitphi)) * pimat.f) ^ (J-n0mat)), na.rm = T) # M x 1
      Af1[Af1 == -Inf] = -10^200; Af2[Af2 == -Inf] = -10^200; # .Machine$double.xmax
      AAf = exp(Af1 + Af2) # M x 1
      AA = theta * AAm + (1-theta) * AAf # M x 1
      theta0 = theta * AAm / AA
      mis = (z == 1) & missing.sex.guys & AA > 0 
      sex[mis] = rbinom(sum(mis), 1, theta0[mis])
      pimat = pimatfn(logitp0, s.left, trap.locations, sex, logsigmam, logsigmaf)
      
      # update theta 
      theta = rbeta( 1, 1 + sum(z*sex), 1 + sum(z*(1 - sex)) )
      
    # update the activity centers
    drawact = drawact +1
    
    s.left.cand = as.matrix(cbind(rnorm(M , s.left[, 1], sigmaOfProposal.s), rnorm(M , s.left[, 2], sigmaOfProposal.s)))
    Scoord = SpatialPoints(s.left.cand*scale, CRS("+proj=utm +zone=43 +datum=WGS84"))  
    SinPoly = over(Scoord, activity) # dim M x 1 
    # Sinpoly[i] equals the index among the mask points if the point i is inside the state space mask polygon, if it is not then Sinpoly = NA
    for (i in c(1:M)[z == 1 & !is.na(SinPoly[,1])]) {
      logsigma_i = ifelse(sex[i] == 1, logsigmam, logsigmaf)
      DD.cand = sqrt((s.left.cand[i,1] - trap.locations[,1]) ^ 2 + (s.left.cand[i,2] - trap.locations[,2]) ^ 2) # a vector Kx1
      pi.cand = expit(logitp0) * exp(-DD.cand * DD.cand / (2 * (exp(logsigma_i) ^ 2 )))
      yyy3 = log(1 + pi.cand^n0mat[i,])
      A3 = log(exp(yyy3) - 1)
      yyy4 = log(1 + (1 - expit(logitphi)* (2 - expit(logitphi)) * pi.cand)^(J - n0mat[i,]))
      A4 = log(exp(yyy4) - 1)
      A3[A3 == -Inf] = -10^200
      A4[A4 == -Inf] = -10^200
      loglik.cand.s = z[i] * (sum(A3+A4))
      yyy3 = log(1 + pimat[i,]^n0mat[i,])
      A3 = log(exp(yyy3) - 1)
      yyy4 = log(1 + (1 - expit(logitphi)* (2 - expit(logitphi)) * pimat[i,])^(J - n0mat[i,]))
      A4 = log(exp(yyy4) - 1)
      A3[A3 == -Inf] = -10^200
      A4[A4 == -Inf] = -10^200
      loglik.curr.s = z[i] * (sum(A3+A4))
      lognum = loglik.cand.s 
      logden = loglik.curr.s
      if (logden == -Inf){
        s.left[i, ] = s.left.cand[i, ]
        naccept.s[i] = naccept.s[i] +1
      }
      if (logden != -Inf){
        logR = lognum - logden
        if (runif(1,0,1) <= exp(logR)){
          s.left[i, ] = s.left.cand[i, ]
          naccept.s[i] = naccept.s[i] +1
        }
      }
      if (floor(draw/batchsize)*batchsize == draw){
        SigmaDiff = ifelse(naccept.s[i] > 0.44*batchsize, exp(2*delta), exp(-2*delta))
        if(draw <= burnin){ sigmaOfProposal.s[i] = sigmaOfProposal.s[i] * SigmaDiff}
        naccept.s[i] = 0 
      } 
  }
    pimat = pimatfn(logitp0, s.left, trap.locations, sex, logsigmam, logsigmaf)
    llik = logLfn.da(ncap, n0mat, n0vec, logitphi, pimat, z, J)
    cat(c(sum(z), psi, sum(sex*z), theta, 
          expit(logitphi), expit(logitp0), 
          exp(logsigmam)*(scale/1000), # output in kilometers
          exp(logsigmaf)*(scale/1000) # output in kilometers
    ), sep = ',', file = paste(folderName, '/markovchain.txt', sep = ""), append = TRUE)
    cat('\n', file = paste(folderName, '/markovchain.txt', sep = ""), append = TRUE)
    cat(c(sex), sep = ',', file = paste(folderName, '/markovchain.sex.txt', sep = ''), append = TRUE) 
    cat('\n', file = paste(folderName, '/markovchain.sex.txt', sep = ''), append = TRUE)
    cat(c(z), sep = ',', file = paste(folderName, '/markovchain.z.txt', sep = ''), append = TRUE)
    cat('\n', file = paste(folderName, '/markovchain.z.txt', sep = ''), append = TRUE)
    cat(c(s.left[,1]), sep = ',', file = paste(folderName, '/markovchain.sx.txt', sep = ''), append = TRUE)
    cat('\n', file = paste(folderName, '/markovchain.sx.txt', sep = ''), append = TRUE)
    cat(c(s.left[,2]), sep = ',', file = paste(folderName, '/markovchain.sy.txt', sep = ''), append = TRUE)
    cat('\n', file = paste(folderName, '/markovchain.sy.txt', sep = ''), append = TRUE)
    cat(llik, sep = ',', file = paste(folderName, '/markovchain.loglikelihood.txt', sep = ''), append = TRUE)
    cat('\n', file = paste(folderName, '/markovchain.loglikelihood.txt', sep = ''), append = TRUE)
    if (known != 'ALL') 
    {
      cat(c(L), sep = ',', file = paste(folderName, '/markovchainL.txt', sep = ""), append = TRUE)
      cat('\n', file = paste(folderName, '/markovchainL.txt', sep = ""), append = TRUE)
    }
  if (draw == ndraws){
    cat('Completed ', ndraws, ' draws of MCMC algorithm', '\n')
    numOfDraws = 0      
    if (numOfDraws == 0) continueMCMC = FALSE
  }
}  # end of MCMC loop

cat('MCMC sampling is completed!', '\n', '\n')
end.time = Sys.time()
(time.taken = end.time - start.time); print(time.taken)
#==================================================================================

post = read.csv(paste(folderName, '/markovchain.txt', sep = ""), sep = ",", header = T)
post = post[(burnin + 1):ndraws,]
geweke_diag = geweke.diag(post, frac1 = 0.2, frac2 = 0.4)
suppressWarnings(write.table(round(unlist(geweke_diag), 3), sep = '           ', 
                             col.names = F, append = F,
                             file = paste(folderName, '/geweke.diag.txt', sep = "")))
N.chain = post[, 'N']
Nvalues = 0:M
probN = rep(0, (M + 1))
for (i in 1:(M + 1)){ probN[i] = length(N.chain[N.chain == (i - 1)])/(ndraws - burnin)}
post.mode.N = Nvalues[probN == max(probN)][1]  
Bayes.Nmean = mean(N.chain, na.rm = T)
Bayes.Nvar = mean(N.chain ^ 2, na.rm = T) - (mean(N.chain, na.rm = T)) ^ 2
ind = cumsum(probN) >= 0.025 & cumsum(probN) <= 0.975
Bayes.Nlower = quantile(N.chain, 0.025, na.rm = T, type = 1) 
Bayes.Nupper = quantile(N.chain, 0.975, na.rm = T, type = 1)
fname5 = paste(folderName, "/mcmc plots of N.jpeg", sep = "")
ylimits = c(0, max(probN))
jpeg(fname5, width = 1000, height = 1000, units = "px", pointsize = 12, quality = 100)
plot(Nvalues[1:M], probN[1:M], type = 'h', ylim = ylimits, xlab = 'N', ylab = 'probability' )
dev.off()
out = as.matrix(c(Bayes.Nmean, sqrt(Bayes.Nvar), Bayes.Nlower, post.mode.N, Bayes.Nupper))
dimnames(out) = list(c('Mean.N', 'SE.N', '2.5%', 'post.mode', '97.5%'),c('HB Estimates of N'))
prob.quantiles = c(0.025, 0.5, 0.975)  
post.stats = cbind(apply(post,2,mean, na.rm = T), apply(post,2,sd, na.rm = T), t(apply(post, 2, quantile, probs = prob.quantiles, na.rm = T,type = 1)) )
prob.names = paste(as.character(100*prob.quantiles), '%', sep = '')
dimnames(post.stats)[2] = list(c('Mean.Chain', 'SD.Chain', prob.names))
mcse.mean.vec = mcse.sd.vec = rep(0, dim(post)[2])
mcse.lq.vec = mcse.uq.vec = mcse.med.vec = rep(0, dim(post)[2])
for (i in 1:dim(post)[2]){
  postvec = post[,i]
  mcse.mean.vec[i] = mcse(postvec[!is.na(postvec)], mean)$se
  mcse.sd.vec[i] = mcse(postvec[!is.na(postvec)], sd)$se
  mcse.med.vec[i] = mcse(postvec[!is.na(postvec)], medQuantile)$se
  mcse.lq.vec[i] = mcse(postvec[!is.na(postvec)], lowerQuantile)$se
  mcse.uq.vec[i] = mcse(postvec[!is.na(postvec)], upperQuantile)$se
}
mcse.mean.vec = unlist(mcse.mean.vec)
mcse.sd.vec = unlist(mcse.sd.vec)
mcse.med.vec = unlist(mcse.med.vec)
mcse.lq.vec = unlist(mcse.lq.vec)
mcse.uq.vec = unlist(mcse.uq.vec)
mcse.mat = cbind(mcse.mean.vec, mcse.sd.vec, mcse.lq.vec, mcse.med.vec, mcse.uq.vec)
HB_estimates = MCSE_estimates = c('', '', '', '', '')
dim.names = c('Mean.Chain', 'SD.Chain', prob.names)
dimnames(mcse.mat) = dimnames(post.stats)
information = as.data.frame(c(M, numl , numr, IDfixed, known, K, J, ndraws, burnin))
dimnames(information) = list(c('M', 'numl', 'numr', 'IDfixed', 'known', 'K', 'J', 'ndraws', 'burnin'), c('info'))
out.final = rbind(t(out), HB_estimates, dim.names, post.stats, MCSE_estimates, dim.names, mcse.mat)
write.csv(rbind(cbind(out.final, matrix('',nrow(out.final), nrow(information) - ncol(out.final))), 
                dimnames(information)[[1]], t(information)), 
          file = paste(folderName, '/EstimatesOfDerivedParam.csv', sep = ""), quote = F,row.names = T)
#============================================================================================================
post = read.csv(paste(folderName, '/markovchain.txt', sep = ""), sep = ",", header = T)
post.mcmc = as.mcmc(post)

fname6 = paste(folderName, "/traceplots_N.jpeg", sep = "")
jpeg(fname6, width = 1000, height = 1000, units = "px", pointsize = 12, quality = 100)
traceplot(post.mcmc[,'N'], xlab = "Iterations", ylab = "N", main = "Traceplot of N",
          cex.main = 2, cex.lab = 2, cex.axis = 2)
dev.off()

fname6 = paste(folderName, "/traceplots_psi.jpeg", sep = "")
jpeg(fname6, width = 1000, height = 1000, units = "px", pointsize = 12, quality = 100)
traceplot(post.mcmc[,'psi'], xlab = "Iterations", ylab = "psi", main = "Traceplot of psi",
          cex.main = 2, cex.lab = 2, cex.axis = 2)
dev.off()

fname6 = paste(folderName, "/traceplots_N.Male.jpeg", sep = "")
jpeg(fname6, width = 1000, height = 1000, units = "px", pointsize = 12, quality = 100)
traceplot(post.mcmc[,'N.Male'], xlab = "Iterations", ylab = "N.Males", main = "Traceplot of N.Males",
          cex.main = 2, cex.lab = 2, cex.axis = 2)
dev.off()

fname6 = paste(folderName, "/traceplots_theta.jpeg", sep = "")
jpeg(fname6, width = 1000, height = 1000, units = "px", pointsize = 12, quality = 100)
traceplot(post.mcmc[,'theta'], xlab = "Iterations", ylab = "theta", main = "Traceplot of theta",
          cex.main = 2, cex.lab = 2, cex.axis = 2)
dev.off()

fname6 = paste(folderName, "/traceplots_phi.jpeg", sep = "")
jpeg(fname6, width = 1000, height = 1000, units = "px", pointsize = 12, quality = 100)
traceplot(post.mcmc[,'phi'], xlab = "Iterations", ylab = "phi", main = "Traceplot of phi",
          cex.main = 2, cex.lab = 2, cex.axis = 2)
dev.off()

fname6 = paste(folderName, "/traceplots_p0.jpeg", sep = "")
jpeg(fname6, width = 1000, height = 1000, units = "px", pointsize = 12, quality = 100)
traceplot(post.mcmc[,'p0'], xlab = "Iterations", ylab = "p0", main = "Traceplot of p0",
          cex.main = 2, cex.lab = 2, cex.axis = 2)
dev.off()

fname6 = paste(folderName, "/traceplots_sigmam.jpeg", sep = "")
jpeg(fname6, width = 1000, height = 1000, units = "px", pointsize = 12, quality = 100)
traceplot(post.mcmc[,'sigmam'], xlab = "Iterations", ylab = "sigma.male", main = "Traceplot of sigma.male",
          cex.main = 2, cex.lab = 2, cex.axis = 2)
dev.off()

fname6 = paste(folderName, "/traceplots_sigmaf.jpeg", sep = "")
jpeg(fname6, width = 1000, height = 1000, units = "px", pointsize = 12, quality = 100)
traceplot(post.mcmc[,'sigmaf'], xlab = "Iterations", ylab = "sigma.female", main = "Traceplot of sigma.female",
          cex.main = 2, cex.lab = 2, cex.axis = 2)
dev.off()

#======================================================
post = read.csv(paste(folderName, '/markovchain.txt', sep = ""), sep = ",", header = T)
post = post[(burnin + 1):ndraws,]
seq1 = seq(1, ndraws - burnin, by = 10)
N1 = post[seq1,'N']; N1 = N1[!is.na(N1)]
psi1 = post[seq1,'psi']; psi1 = psi1[!is.na(psi1)]
N.Male1 = post[seq1,'N.Male']; N.Male1 = N.Male1[!is.na(N.Male1)]
theta1 = post[seq1,'theta']; theta1 = theta1[!is.na(theta1)]
phi1 = post[seq1,'phi']; phi1 = phi1[!is.na(phi1)]
p01 = post[seq1,'p0']; p01 = p01[!is.na(p01)]
sigmam1 = post[seq1,'sigmam']; sigmam1 = sigmam1[!is.na(sigmam1)]
sigmaf1 = post[seq1,'sigmaf']; sigmaf1 = sigmaf1[!is.na(sigmaf1)]
df = data.frame(N1, psi1, N.Male1, theta1, phi1, p01, sigmam1, sigmaf1)
fname6 = paste(folderName, "/Scatter Plots_", IDfixed,known, ".jpeg", sep = "")
jpeg(fname6, width = 1000, height = 1000, units = "px", pointsize = 12, quality = 100)
pairs(~N1 + psi1 + N.Male1 + theta1 + phi1 + p01 + sigmam1 + sigmaf1, data = df, main = "")
dev.off()
#========================================

save.image(paste(folderName, '/savingRimage.RData', sep = ""))

} # end of function PISCR.fn
