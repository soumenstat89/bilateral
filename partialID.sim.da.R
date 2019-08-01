library(MCMCpack)
library(abind)

source("partialID.simdata.R")
source("utility.functions.R")
source("Lvec.R")

ndrawss = rep(100, 70)
burninn = rep(10, 70)
MM = rep(400,70)
NN = rep(100, 70)
NN.Male = rep(40, 70)
phiphi = rep(sort(rep(seq(0.3, 0.9, by = 0.1), 5)), 2)
p0p0 = rep(c(0.005, seq(0.01, 0.07, by =0.02)), 14)

sigm = c(rep(0.3, 35), rep(0.4, 35)) # in km as state space
sigf = c(rep(0.15, 35), rep(0.2, 35)) # in km as state pace

sim.mat = cbind(ndrawss, burninn, MM, NN, NN.Male, phiphi, p0p0, sigm, sigf)

simdata = matrix(NA, 70, dim(sim.mat)[2] + 4)

# for(j in 1:70) {

j = 35
cat("j = ", j, '\n', sep = "")


M = sim.mat[j,3]

batchsize = 2000
mindeltaLogSigmaForRW = 0.01

sigmaOfProposal.logitphi = 0.02
sigmaOfProposal.logitp0 = 0.02
sigmaOfProposal.logsigmam = 0.01 
sigmaOfProposal.logsigmaf = 0.01 
sigmaOfProposal.L = 4
sigmaOfProposal.s = rep(3, M)

ndraws = sim.mat[j,1]; burnin = sim.mat[j,2]; 

nloopL = 20; n.update = 10
JJ = 50; plotitt = F;

nrow_trapp = 10; ncol_trapp = 16; xlimm = c(0,5); ylimm = c(0,7); bufferr = 1;


#=====================
# For simulated data
#====================


data = sim.gender.partial.data(N = sim.mat[j,4], N.Male = sim.mat[j,5], phi = sim.mat[j,6], 
                            p0 = sim.mat[j,7], sigmam = sim.mat[j,8], sigmaf = sim.mat[j,9], 
                            nrow_trap = nrow_trapp, ncol_trap = ncol_trapp, J = JJ, 
                            xlim = xlimm, ylim = ylimm, buffer = bufferr) 

plotit = plotitt


N.given = data$N; N.Male.given = data$N.Male; 

phi.given = data$phi; p0.given = data$p0

sigmam.given = data$sigmam; sigmaf.given = data$sigmaf

(numl = data$numl);  (numr = data$numr)

left.obs = data$left; right.obs = data$right

IDfixed = data$IDfixed; known = data$known
simdata[j,] = c(sim.mat[j,], IDfixed, numl,  numr,  known)

print(data.frame(N.given,  numl,  numr,  known,  IDfixed))

dimnames(simdata) = list(c(1:dim(simdata)[1]), c('ndraws', 'burnin', 'M', 'N', 'N.Male',
                                                 'phi', 'p0', 'sigmam', 'sigmaf', 'IDfixed', 'numl', 'numr', 'known'))

xlim = data$xlim 
ylim = data$ylim 
buffer = data$buffer 

trap.locations = data$trap.locations 
K = data$K; J = data$J


left = abind(left.obs, array(0, dim = c( M - numl, K, J)), along = 1)
right = abind(right.obs, array(0, dim = c( M - numr, K, J)), along = 1)


#######################################################
## HANDLING THE SEX INFORMATION ##

theta = 0.6  # Initialising theta

sexl.obs = sexl = data$sexl # length numl, but only first IDfixed no. of entries has values, rest NA
sexr.obs = sexr = data$sexr # length numr, but only first IDfixed no. of entries has values, rest NA

sex.true = data$sex.true # length N.given, the true sexes of the N.given individuals

# CHECKING
# Now no. of males and females in sexl and sexr are same.
# Hence no. of missing elements in sexl and sexr are also same.

unique(sexl[1:IDfixed] == sexr[1:IDfixed]) # TRUE

sum(sexl==1, na.rm = T) == sum(sexr==1, na.rm = T) # can be FALSE due to difference in ordering and number of left captured individuals
sum(sexl==0, na.rm = T) == sum(sexr==0, na.rm = T) # can be FALSE due to difference in ordering and number of right captured individuals


# Initializing the sex for (M - IDfixed) no. of non IDfixed guys
numlm = sum(sexl == 1, na.rm = T)
numlf = sum(sexl == 0, na.rm = T)

numrm = sum(sexr == 1, na.rm = T)
numrf = sum(sexr == 0, na.rm = T)


if(numlm < numrm) sexl = c(sexl, rep(1, numrm - numlm))
if(numrm < numlm) sexr = c(sexr, rep(1, numlm - numrm))
if(numlf < numrf) sexl = c(sexl, rep(0, numrf - numlf))
if(numrf < numlf) sexr = c(sexr, rep(0, numlf - numrf))

nnn = max(numlm, numrm) + max(numlf, numrf)

if(numl < nnn) sexl[(numl+1) : nnn] = sample(sexl[(numl+1) : nnn], length(nnn - numl), replace = F)

if(numr < nnn) sexr[(numr+1) : nnn] = sample(sexr[(numr+1) : nnn], length(nnn - numr), replace = F)

sex.allocation = rbinom(M - (max(numlm, numrm) + max(numlf, numrf)), 1, theta)

sexl = c(sexl, sex.allocation) # length M and have equal number of males and females
sexr = c(sexr, sex.allocation) # length M and have equal number of males and females

sex = sexl

# CHECKING
# Now no. of males and females in sexl and sexr are same.
# Hence no. of missing elements in sexl and sexr are also same.

sum(sexl==1) == sum(sexr==1) # TRUE
sum(sexl==0) == sum(sexr==0) # TRUE
sexl[1:nnn] == sexr[1:nnn] # not all true 

missing.sex.guys = c(rep(F, numl), rep(T, M - numl)) # the guys whose sex are not known

#######################################################

# Initial values for L vector of length M

if (known != 'ALL') L = Lvec(left, right, numl, numr, sexl, sexr, trap.locations, IDfixed, nloopL)

#==============================

s.left = cbind(runif(M, xlim[1], xlim[2]), runif(M, ylim[1], ylim[2]))

left2d = apply(left,c(1,2),sum) # no. of times individual i got captured at trap k on left side

for (i in 1:M)
{
  if (sum(left2d[i,]) == 0) next
  # get the traps where ith LSI got captured
  traps.loc = matrix(trap.locations[left2d[i,] > 0,], ncol = 2, byrow = F)
  # trps<-matrix(traps,ncol = 2,byrow = FALSE)
  s.left[i,] = apply(traps.loc, 2, weighted.mean, w = left2d[i,][left2d[i,] > 0]) #c(mean(trps[,1]),mean(trps[,2]))
  # s.left[i,] = apply(traps, 2, mean) #c(mean(trps[,1]),mean(trps[,2]))
}
#==============================


#==============================

# zero.guys are the (non - IDfixed) guys who never get captured either on left or right side.
# Some of them may be ghost some are not depending upon z_i = 0 or z_i = 1.
# The guys who are captured for them z_i = 1. (IDfixed guys have always z_i = 1)
# There are no zero guys among IDfixed guys.

ncapl = rowSums(left) # apply(left, 1, sum)
ncapr = rowSums(right) # apply(right, 1, sum)

if (known != 'ALL') {ncapr.star = ncapr[order(L)];  right.star = right[order(L),,]}
if (known == 'ALL') {ncapr.star = ncapr;            right.star = right}
ncap = ncapl + ncapr.star
n0mat = apply(left + right.star, c(1,2), function(a){sum(a > 0)})
n0vec = rowSums(n0mat) # apply(n0mat, 1, sum)

zero.guys = (ncapl + ncapr.star) == 0 # = (ncapl + ncapr.star) == 0 & c(1:M) > IDfixed
numcap = sum((ncapl + ncapr.star) > 0) # = M - sum(zero.guys) # = IDfixed + sum((ncapl + ncapr.star) > 0 & c(1:M) > IDfixed) 
# n.uncap = sum(z) - numcap

numcap >= max(numl, numr)   # max(numl, numr) <= numcap <= numl + numr + IDfixed
# if known == 'ALL', numcap = numlthis is true, if known != 'ALL', this is false

# Initial values
logitphi = logit(0.8); logitp0 = logit(0.05)
logsigmam = log(0.1); logsigmaf = log(0.5)
psi = 0.5

z = c(rep(1, numcap), rbinom(M - numcap, 1, psi))

# .... assign values for adaptive Metropolis - Hastings sampling
batch = 0
delta = mindeltaLogSigmaForRW

# Setting counters

naccept.logitphi = 0
naccept.logitp0 = 0
naccept.logsigmam = 0
naccept.logsigmaf = 0
naccept.L = 0
naccept.s = rep(0,M)

ts = format(Sys.time(), "%d%m%y_%H%M%S")
folderName = paste0("partialID_", ts, sep = "")
dir.create(path = folderName)

start.time = Sys.time()

# .... compute MCMC draws

continueGibbs = TRUE
draw = 0
drawphi = 0
drawp0 = 0
drawsigmam = 0
drawsigmaf = 0
drawL = 0
drawsex = 0
drawact = 0

cat('Begin MCMC sampling:', '\n', '\n')

fname1 = paste(folderName, '/sigma_proposals.txt', sep = "")
cat(c(paste("sigmaOfProposal.logitphi = ", sigmaOfProposal.logitphi, sep = ''), '\n'), sep = ',', file = fname1, append = TRUE)
cat(c(paste("sigmaOfProposal.logitp0 = ", sigmaOfProposal.logitp0, sep = ''), '\n'), sep = '', file = fname1, append = TRUE)
cat(c(paste("sigmaOfProposal.logsigmam = ", sigmaOfProposal.logsigmam, sep = ''), '\n'), sep = '', file = fname1, append = TRUE)
cat(c(paste("sigmaOfProposal.logsigmaf = ", sigmaOfProposal.logsigmaf, sep = ''), '\n'), sep = '', file = fname1, append = TRUE)
cat(c(paste("sigmaOfProposal.L = ", sigmaOfProposal.L, sep = ''), '\n'), sep = '', file = fname1, append = TRUE)
cat(c(paste("sigmaOfProposal.s = rep(", sigmaOfProposal.s[1], ",M)", sep = ''), '\n'), sep = '', file = fname1, append = TRUE)

if (known != 'ALL') 
{
  cat(c(paste("L", 1:M, sep = '')), sep = ',', file = paste(folderName, '/markovchainL.txt', sep = ""), append = TRUE)
  cat('\n', file = paste0(folderName, '/markovchainL.txt', sep = ""), append = TRUE)
}

cat(c('N', 'psi', 'numcap', 'N.Male', 'theta', 'phi', 'p0', 'sigmam', 'sigmaf'), sep = ',', file = paste0(folderName, '/markovchain.txt', sep = ""), append = TRUE)
cat('\n', file = paste0(folderName, '/markovchain.txt', sep = ""), append = TRUE)

while (continueGibbs) {
  
  draw = draw + 1
  drawinterval = 500
  if (draw == round(draw/drawinterval)*drawinterval)  cat('..... drawing sample #', draw, '\n')
  
  # update the increment/decrement for adaptive Metropolis - Hastings samplers
  if (floor(draw/batchsize)*batchsize == draw) 
  {
    batch = batch + 1
    if (1/sqrt(batch) < mindeltaLogSigmaForRW)  delta = 1/sqrt(batch)
  }
  
  
  # update phi
  
  drawphi = drawphi +1
  # prob = c(phi ^ 2, 2*phi*(1 - phi), (1 - phi)^2)
  pimat = pimatfn(logitp0, s.left, trap.locations, sex, logsigmam, logsigmaf)
  loglik.curr.phi = logLfn.da(ncap, n0mat, n0vec, logitphi, pimat, z, J)
 
  logitphi.cand = rnorm(1, logitphi, sigmaOfProposal.logitphi) 
  loglik.cand.phi = logLfn.da(ncap, n0mat, n0vec, logitphi.cand, pimat, z, J)
  
  lognum = loglik.cand.phi # + dunif(logit.phi.cand, log = T) + dnorm(logitphi, logitphi.cand, sigmaOfProposal.logitphi, log = T))  ## log-prior + log-proposal
  logden = loglik.curr.phi # + dunif(logit.phi     , log = T) + dnorm(logitphi.cand, logitphi, sigmaOfProposal.logitphi, log = T))  ## log-prior + log-proposal
  
  if (logden == -Inf)
  {
    logitphi = logitphi.cand
    loglik.curr.phi = loglik.cand.phi
    naccept.logitphi = naccept.logitphi + 1
  }  
  if (logden != -Inf)
  {
    logR = lognum - logden
    if (runif(1,0,1) <= exp(logR)) 
    {
      logitphi = logitphi.cand
      loglik.curr.phi = loglik.cand.phi
      naccept.logitphi = naccept.logitphi + 1
    }
  }
  
  if (floor(draw/batchsize)*batchsize == draw)
  {
    SigmaDiff = ifelse(naccept.logitphi > 0.44*batchsize, exp(2*delta), exp(-2*delta))
    if(draw <= burnin){ sigmaOfProposal.logitphi = sigmaOfProposal.logitphi * SigmaDiff}
    cat("proposal sd of logitphi = ", sigmaOfProposal.logitphi, ' ')
    cat("naccept.logitphi = ", naccept.logitphi, '\n')
    naccept.logitphi = 0   # reset counter for next batch
  }
  # update p0
  
  drawp0 = drawp0 +1
  # prob = c(p0 ^ 2, 2*p0*(1 - p0), (1 - p0)^2)
  loglik.curr.p0 = loglik.curr.phi
  
  logitp0.cand = rnorm(1, logitp0, sigmaOfProposal.logitp0) 
  # print(p0.cand)
  pimat.cand = pimatfn(logitp0.cand, s.left, trap.locations, sex, logsigmam, logsigmaf)
  loglik.cand.p0 = logLfn.da(ncap, n0mat, n0vec, logitphi, pimat.cand, z, J)
  
  lognum = loglik.cand.p0 # + dunif(logit.p0.cand, log = T) + dnorm(logitp0, logitp0.cand, sigmaOfProposal.logitp0, log = T))  ## log-prior + log-proposal
  logden = loglik.curr.p0 # + dunif(logit.p0     , log = T) + dnorm(logitp0.cand, logitp0, sigmaOfProposal.logitp0, log = T))  ## log-prior + log-proposal
  
  if (logden == -Inf)
  {
    logitp0 = logitp0.cand
    pimat = pimat.cand
    loglik.curr.p0 = loglik.cand.p0
    naccept.logitp0 = naccept.logitp0 + 1
  }
  if (logden != -Inf)
  {
    logR = lognum - logden
    if (runif(1,0,1) <= exp(logR)) 
    {
      logitp0 = logitp0.cand
      pimat = pimat.cand
      loglik.curr.p0 = loglik.cand.p0
      naccept.logitp0 = naccept.logitp0 + 1
    }
  }
  
  if (floor(draw/batchsize)*batchsize == draw)
  {
    SigmaDiff = ifelse(naccept.logitp0 > 0.44*batchsize, exp(2*delta), exp(-2*delta))
    if(draw <= burnin){ sigmaOfProposal.logitp0 = sigmaOfProposal.logitp0 * SigmaDiff}
    cat("proposal sd of logitp0 = ", sigmaOfProposal.logitp0, ' ')
    cat("naccept.logitp0 = ", naccept.logitp0, '\n')
    naccept.logitp0 = 0   # reset counter for next batch
  }
  
  
  # update sigmam
  loglik.curr.sigmam = loglik.curr.p0
  
  drawsigmam = drawsigmam +1
  
  logsigmam.cand = rnorm(1, logsigmam, sigmaOfProposal.logsigmam)
  													
  pimat.cand = pimatfn(logitp0, s.left, trap.locations, sex, logsigmam.cand, logsigmaf)
  
  loglik.cand.sigmam = logLfn.da(ncap, n0mat, n0vec, logitphi, pimat.cand, z, J)
  
  lognum = loglik.cand.sigmam # + dnorm(logsigmam, logsigmam.cand, sigmaOfProposal.logsigmam, log = T) ## proposal
  
  logden = loglik.curr.sigmam # + dnorm(logsigmam.cand, logsigmam, sigmaOfProposal.logsigmam, log = T) ## proposal
  
  if (logden == -Inf)
  {
    logsigmam = logsigmam.cand
    pimat = pimat.cand
    loglik.curr.sigmam = loglik.cand.sigmam
    naccept.logsigmam = naccept.logsigmam + 1
  }
  if (logden != -Inf)
  {
    logR = lognum - logden
    if (runif(1,0,1) <= exp(logR)) 
    {
      logsigmam = logsigmam.cand
      pimat = pimat.cand
      loglik.curr.sigmam = loglik.cand.sigmam
      naccept.logsigmam = naccept.logsigmam + 1
    }
  }
  # }
  if (floor(draw/batchsize)*batchsize == draw) 
  {
    SigmaDiff = ifelse(naccept.logsigmam > 0.44*batchsize, exp(2*delta), exp(-2*delta))
    if(draw <= burnin){ sigmaOfProposal.logsigmam = sigmaOfProposal.logsigmam * SigmaDiff}
    cat("proposal sd of logsigmam = ", sigmaOfProposal.logsigmam, ' ')
    cat("naccept.logsigmam = ", naccept.logsigmam, '\n')
    naccept.logsigmam = 0   # reset counter for next batch
  }
  
  # update sigmaf
  loglik.curr.sigmaf = loglik.curr.sigmam
  
  drawsigmaf = drawsigmaf +1
  
  logsigmaf.cand = rnorm(1, logsigmaf, sigmaOfProposal.logsigmaf)
  												
  pimat.cand = pimatfn(logitp0, s.left, trap.locations, sex, logsigmam, logsigmaf.cand)
  
  loglik.cand.sigmaf = logLfn.da(ncap, n0mat, n0vec, logitphi, pimat.cand, z, J)
  
  lognum = loglik.cand.sigmaf # + dnorm(logsigmaf, logsigmamf.cand, sigmaOfProposal.logsigmaf, log = T) ## proposal
  
  logden = loglik.curr.sigmaf # + dnorm(logsigmaf.cand, logsigmaf, sigmaOfProposal.logsigmaf, log = T) ## proposal
  
  if (logden == -Inf) 
  {
    logsigmaf = logsigmaf.cand
    pimat = pimat.cand
    loglik.curr.sigmaf = loglik.cand.sigmaf
    naccept.logsigmaf = naccept.logsigmaf + 1
  }
  if (logden != -Inf)
  {
    logR = lognum - logden
    if (runif(1,0,1) <= exp(logR)) 
    {
      logsigmaf = logsigmaf.cand
      pimat = pimat.cand
      loglik.curr.sigmaf = loglik.cand.sigmaf
      naccept.logsigmaf = naccept.logsigmaf + 1
    }
  }
  # }
  if (floor(draw/batchsize)*batchsize == draw) 
  {
    SigmaDiff = ifelse(naccept.logsigmaf > 0.44*batchsize, exp(2*delta), exp(-2*delta))
    if(draw <= burnin){ sigmaOfProposal.logsigmaf = sigmaOfProposal.logsigmaf * SigmaDiff}
    cat("proposal sd of logsigmaf = ", sigmaOfProposal.logsigmaf, ' ')
    cat("naccept.logsigmaf = ", naccept.logsigmaf, '\n')
    naccept.logsigmaf = 0   # reset counter for next batch
  }
  
  #  update L
  
  if (known != "ALL")
  {
    loglik.curr.L= loglik.curr.sigmaf
    drawL = drawL +1
    indx1 = (IDfixed + 1):M
    
    # n.update = 20 # 10, 20, 30
    if( IDfixed < numr) rightset1  = order(L)[z == 1 & order(L) >= (IDfixed+1) & order(L) <= numr]
    if( IDfixed >= numr) rightset1  = c()
    
    rightset2 = order(L)[z == 1 & order(L) > numr] # real right individuals wh went uncaptured
    
    if(length(rightset2) >1) rightset2 = sample(rightset2, min(length(rightset2), n.update), replace = F)
    # rightset2 = order(L)[z == 1]; rightset2 = rightset2[rightset2 > numr]
    rightset3 = c()
    if(IDfixed < numl )
    {
      leftset3 = c(1:M)[z == 1 & c(1:M)>IDfixed & c(1:M) <= numl]
      if(length(leftset3) > 0) for(ii in leftset3) rightset3 = c(rightset3, c(1:M)[L == ii]) 
      # right individuals to whom left captured individuals are linked with
    }
    
    rightset = sort(unique(c(rightset1, rightset2, rightset3)))
    
    for (r.guy1 in rightset)
    {
      
      l.swap.out = L[r.guy1]
      # if (z[l.swap.out]==0) next
      
      #  q(theta, theta_cand)
      possible.L = c(1:M)[z == 1 & c(1:M) > IDfixed & sex == sex[l.swap.out]]
      dv =  sqrt( (s.left[l.swap.out,1] - s.left[,1]) ^ 2 + (s.left[l.swap.out,2] - s.left[,2]) ^ 2 )
      wt.possible.L = exp( - (dv ^ 2) / sigmaOfProposal.L ^ 2)[z == 1 & c(1:M) > IDfixed  & sex == sex[l.swap.out]]
      
      if (length(possible.L) > 1)  l.swap.in =  sample( possible.L, 1, replace = F, prob = wt.possible.L)
      if (length(possible.L) == 1) next 
      if (length(possible.L) == 0) next # this case will never happen since l.swap.out is present there at the centre of the circle dv < 5
      if (l.swap.in == l.swap.out) next # saves computation time in for loop
      jump.prob.L =  wt.possible.L[which(possible.L == l.swap.in, arr.ind = T)] / sum(wt.possible.L) #  q(state.curr, state.cand)
      
      #  q(theta_cand, theta) 
      # possible.back = c(1:M)[ z == 1 & c(1:M) > IDfixed] 
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
      
      if (plotit)
      { 
        #         plot(s.left[z==1])
        #         if (sex[l.swap.out]==1) points(s.left[c(l.swap.in,l.swap.out),], pch = c(24, 22), bg = 'blue', xlab = 'Easting', ylab = 'Northing')
        #         if (sex[l.swap.out]==0) points(s.left[c(l.swap.in,l.swap.out),], pch = c(24, 22), bg = 'red', xlab = 'Easting', ylab = 'Northing')
        #         
        # pch = 22 is square, pch = 24 is triangle.
        plot(s.left[z == 1])
        points(s.left[c(l.swap.in,l.swap.out),], pch = c(24, 22), bg = 'blue', xlab = 'Easting', ylab = 'Northing')
        
      }
      
      # loglik.curr = logLfn.da(ncapl, ncapr.star, n0mat, n0vec, phi, pimat, z, J)
      right.star.cand = right[order(L.cand),,]
      ncapr.star.cand = ncapr[order(L.cand)]
      ncap.cand = ncapl + ncapr.star.cand
      n0mat.cand = apply(left + right.star.cand, c(1,2), function(a){sum(a > 0)})
      n0vec.cand = rowSums(n0mat.cand) # apply(n0mat.cand, 1, sum)
      loglik.cand.L = logLfn.da(ncap.cand, n0mat.cand, n0vec.cand, logitphi, pimat, z, J)
      
      lognum = loglik.cand.L + log(jump.prob.back.L)
      
      logden = loglik.curr.L + log(jump.prob.L)
      
      if (logden == -Inf) 
      {
        L = L.cand 
        loglik.curr.L = loglik.cand.L
        right.star = right.star.cand
        ncapr.star = ncapr.star.cand
        ncap = ncap.cand
        n0mat = n0mat.cand
        n0vec = n0vec.cand
        naccept.L = naccept.L + 1
      }
      if (logden != -Inf)
      {
        logR = lognum - logden
        if (runif(1,0,1) <= exp(logR)) 
        {
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
    
    
    
    if (floor(draw/batchsize)*batchsize == draw) 
    {
      SigmaDiff = ifelse(naccept.L > 0.44*batchsize, exp(2*delta), exp(-2*delta))
      if(draw <= burnin){ sigmaOfProposal.L= sigmaOfProposal.L * SigmaDiff}
      cat(paste("proposal sd of L = ", sep = ""), sigmaOfProposal.L, ' ')
      cat(paste("naccept.L = ", sep = ""), naccept.L, '\n')
      naccept.L = 0   # reset counter for next batch
    } 
    
    #======================================
    # update z
    
    zero.guys = (ncapl + ncapr.star) == 0 # = (ncapl + ncapr.star) == 0 & c(1:M) > IDfixed
    numcap = sum((ncapl + ncapr.star) > 0) # = M - sum(zero.guys) # = IDfixed + sum((ncapl + ncapr.star) > 0 & c(1:M) > IDfixed) 
    # n.uncap = sum(z) - numcap
    
    ncprob = (1 - expit(logitphi) * (2 - expit(logitphi)) * pimat) ^ J # = ((1 - phi)^2)*pimat + (1 - pimat)
    # prob0 = exp(apply(log(ncprob), 1, sum))  # *(theta ^ sex)*((1 - theta)^(1 - sex))
    prob0 = exp(rowSums(log(ncprob))) # *(theta ^ sex)*((1 - theta)^(1 - sex))
    
    fc = prob0*psi / (prob0*psi + 1 - psi)
    z[zero.guys] = rbinom(sum(zero.guys), 1, fc[zero.guys])
    z[!zero.guys] = 1
    
    # update psi
    
    # beta(1, 1) prior for psi # beta(0.01, 1) prior can also be used
    psi = rbeta(1, 1 + sum(z), 1 +  (M - sum(z) ) )  # count IDfixed guys
    
    
  } # end of if (known!='ALL')
  
  if (known == 'ALL')
  {
    # update z
    
    zero.guys = (ncapl + ncapr.star) == 0 # == (ncapl + ncapr.star) == 0 & c(1:M) > IDfixed
    # numcap does not change in known  = "ALL" situation
    # numcap = sum((ncapl + ncapr.star) > 0) # = M - sum(zero.guys) # = IDfixed + sum((ncapl + ncapr.star) > 0 & c(1:M) > IDfixed) 
    # n.uncap = sum(z) - numcap
    
    ncprob = (1 - expit(logitphi) * (2 - expit(logitphi)) * pimat) ^ J # = (1 - pimat)+ ((1 - phi)^2)*pimat
    # prob0 = exp(apply(log(ncprob), 1, sum))  #*(theta ^ sex)*((1 - theta)^(1 - sex))
    prob0 = exp(rowSums(log(ncprob))) # *(theta ^ sex)*((1 - theta)^(1 - sex))
    
    fc = prob0*psi / (prob0*psi + 1 - psi)
    z[zero.guys] = rbinom(sum(zero.guys), 1, fc[zero.guys])
    z[!zero.guys] = 1
    
    # update psi
    
    # beta(1, 1) prior for psi # beta(0.01, 1) prior can also be used
    psi = rbeta(1, 1 + sum(z), 1 +  (M - sum(z) ) )
    
  } # end of if (known=='ALL')
  
  # update sex
  
  drawsex = drawsex +1
  sex.cand = sex
  for(i in  c(1:M)[missing.sex.guys & z == 1] )
  {
    sex.cand[i]= rbinom(1, 1, theta) # proposal
    
    DD = sqrt((s.left[i,1] - trap.locations[,1]) ^ 2 + (s.left[i,2] - trap.locations[,2]) ^ 2) # a vector Kx1
    
    logsigma.cand = ifelse(sex.cand[i] == 1, logsigmam, logsigmaf)
    pi.cand = expit(logitp0) * exp(-DD * DD / (2 * (exp(logsigma.cand) ^ 2)))
    
    yyy3 = log(1 + pi.cand^n0mat[i,])
    A3 = log(exp(yyy3) - 1)
    yyy4 = log(1 + (1 - expit(logitphi)* (2 - expit(logitphi)) * pi.cand)^(J - n0mat[i,]))
    A4 = log(exp(yyy4) - 1)
    A3[A3 == -Inf] = -10^200
    A4[A4 == -Inf] = -10^200
    loglik.cand.sex = z[i] * (sum(A3+A4))
    
    yyy3 = log(1 + pimat[i,]^n0mat[i,])
    A3 = log(exp(yyy3) - 1)
    yyy4 = log(1 + (1 - expit(logitphi)* (2 - expit(logitphi)) * pimat[i,])^(J - n0mat[i,]))
    A4 = log(exp(yyy4) - 1)
    A3[A3 == -Inf] = -10^200
    A4[A4 == -Inf] = -10^200
    loglik.curr.sex = z[i] * (sum(A3+A4))
    
    lognum = loglik.cand.sex # + dbinom(sex.cand[i], theta, log = T) + ## prior
    logden = loglik.curr.sex # + dbinom(sex[i], theta, log = T) + ## prior
    
    if (logden == -Inf)
    {
      sex[i] = sex.cand[i]
    }
    if (logden != -Inf)
    {
      logR = lognum - logden
      if (runif(1,0,1) <= exp(logR)) 
      {
        sex[i] = sex.cand[i]
      }
    }
  }
 														
  pimat = pimatfn(logitp0, s.left, trap.locations, sex, logsigmam, logsigmaf)
  
  # update theta
  
  theta = rbeta( 1, 1 + sum(z*sex), 1 + sum(z*(1 - sex)) )
  
  # update the activity centers
  drawact = drawact +1
  
  s.left.cand = matrix(c(rnorm(M, s.left[, 1], sigmaOfProposal.s), rnorm(M, s.left[, 2], sigmaOfProposal.s)), nrow = M, ncol = 2 , byrow = F)
  inbox = s.left.cand[,1] < xlim[2] & s.left.cand[,1] > xlim[1] & s.left.cand[,2] < ylim[2] & s.left.cand[,2] > ylim[1]
  
  for (i in c(1:M)[z == 1 & inbox])  
  {
   # if (!inbox[i]) next 
    
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
      
      lognum = loglik.cand.s #+ log(jump.prob.back.s)
      
      logden = loglik.curr.s #+ log(jump.prob.s)
      
      if (logden == -Inf)
      {
        s.left[i, ] = s.left.cand[i, ]
        naccept.s[i] = naccept.s[i] +1
      }
      if (logden != -Inf)
      {
        logR = lognum - logden
        if (runif(1,0,1) <= exp(logR)) 
        {
          s.left[i, ] = s.left.cand[i, ]
          naccept.s[i] = naccept.s[i] +1
        }
      }
      
      if (floor(draw/batchsize)*batchsize == draw) 
      {
        SigmaDiff = ifelse(naccept.s[i] > 0.44*batchsize, exp(2*delta), exp(-2*delta))
        if(draw <= burnin){ sigmaOfProposal.s[i] = sigmaOfProposal.s[i] * SigmaDiff}
        # cat(paste("naccept.s", i, " = ", sep = ""), naccept.s[i], '\n')
        # cat(paste("proposal sd of s", i, " = ", sep = ""), sigmaOfProposal.s[i], '\n')
        naccept.s[i] = 0   # reset counter for next batch
      } 
  } # end of for (i in c(1:M)[z == 1])
  
  if (known != 'ALL') 
  {
    cat(c(L), sep = ',', file = paste(folderName, '/markovchainL.txt', sep = ""), append = TRUE)
    cat('\n', file = paste(folderName, '/markovchainL.txt', sep = ""), append = TRUE)
  }
  
  cat(c(sum(z), psi, numcap, sum(sex*z), theta, 
        expit(logitphi), expit(logitp0), 
        exp(logsigmam), 
        exp(logsigmaf)
  ), sep = ',', file = paste(folderName, '/markovchain.txt', sep = ""), append = TRUE)
  cat('\n', file = paste(folderName, '/markovchain.txt', sep = ""), append = TRUE)
  
  if (draw == ndraws) 
  {
    cat('Completed ', ndraws, ' draws of MCMC algorithm', '\n')
    # numOfDraws = as.integer(readline(prompt = 'Enter additional number of MCMC draws -> '))
    numOfDraws = 0      
    if (numOfDraws == 0) continueGibbs = FALSE
    # else ndraws = ndraws + numOfDraws
    cat(c(paste('\n', 'Completed ', ndraws, ' draws of MCMC algorithm', sep = ''), '\n \n'), sep = ',', file = fname1, append = TRUE)
    cat(c(paste("sigmaOfProposal.logitphi = ", sigmaOfProposal.logitphi, sep = ''), '\n'), sep = ',', file = fname1, append = TRUE)
    cat(c(paste("sigmaOfProposal.logitp0 = ", sigmaOfProposal.logitp0, sep = ''), '\n'), sep = '', file = fname1, append = TRUE)
    cat(c(paste("sigmaOfProposal.logsigmam = ", sigmaOfProposal.logsigmam, sep = ''), '\n'), sep = '', file = fname1, append = TRUE)
    cat(c(paste("sigmaOfProposal.logsigmaf = ", sigmaOfProposal.logsigmaf, sep = ''), '\n'), sep = '', file = fname1, append = TRUE)
    cat(c(paste("sigmaOfProposal.L = ", sigmaOfProposal.L, sep = ''), '\n'), sep = '', file = fname1, append = TRUE)
    cat(c(paste("sigmaOfProposal.s = rep(", sigmaOfProposal.s[1], ",M)", sep = ''), '\n'), sep = '', file = fname1, append = TRUE)
    
  }
  
}  # end of while loop

cat('Gibbs sampling is completed!', '\n', '\n')
end.time = Sys.time()
(time.taken = end.time - start.time); print(time.taken)
#==================================================================================

# Compute estimates of population abundance for a specific region within S

# .... compute Bayes estimate of posterior mean N and estimate of posterior probabilities

post = read.csv(paste0(folderName, '/markovchain.txt', sep = ""), sep = ",", header = T)
post = post[(burnin + 1):ndraws,]
geweke_diag = geweke.diag(post, frac1 = 0.2, frac2 = 0.4)
print(geweke_diag)
suppressWarnings(write.table(round(unlist(geweke_diag), 3), sep = '           ', 
                             col.names = F, append = F,
                             file = paste0(folderName, '/geweke.diag.txt', sep = "")))

N.chain = post[, 'N']

Nvalues = 0:M
probN = rep(0, (M + 1))
for (i in 1:(M + 1)) probN[i] = length(N.chain[N.chain == (i - 1)])/(ndraws - burnin)

post.mode.N = Nvalues[probN == max(probN)][1]  # which(probN==max(probN), arr.ind = T)-1

Bayes.Nmean = mean(N.chain, na.rm = T)
Bayes.Nvar = mean(N.chain ^ 2, na.rm = T) - (mean(N.chain, na.rm = T)) ^ 2
ind = cumsum(probN) >= 0.025 & cumsum(probN) <= 0.975
Bayes.Nlower = quantile(N.chain, 0.025, na.rm = T, type = 1) # min(Nvalues[ind])
Bayes.Nupper = quantile(N.chain, 0.975, na.rm = T, type = 1) # max(Nvalues[ind])

# # # .... plot histograms of empirical Bayes and Bayes estimates
fname5 = paste0(folderName, "/mcmc plots of N.jpeg", sep = "")
ylimits = c(0, max(probN))
jpeg(fname5, width = 1000, height = 1000, units = "px", pointsize = 12, quality = 100)
# # plot(n.uncap + n, Bayes.probN, type = 'h', ylim = ylimits); lines(n.uncap + n, dlnorm(n.uncap, meanlog = 3.13, sdlog = 0.315), col = 'red', lwd = 2)
plot(Nvalues[1:M], probN[1:M], type = 'h', ylim = ylimits, xlab = 'N', ylab = 'probability' )
# lines(n.uncap + n, dlnorm(n.uncap, meanlog = log(Bayes.n0Mean), sdlog = sqrt(Bayes.n0Var)), col = 'red', lwd = 2)
dev.off()


# Print results of Bayesian analysis

out = as.matrix(c(Bayes.Nmean, sqrt(Bayes.Nvar), Bayes.Nlower, post.mode.N, Bayes.Nupper))
dimnames(out) = list(c('Mean.N', 'SE.N', '2.5%', 'post.mode', '97.5%'),c('HB Estimates of N'))

# Compute posterior means and quantiles
if(known!='ALL') numcap.given = length(unique(c(data$which.left.obsd, data$which.right.obsd)))
if(known=='ALL') numcap.given = numl

par_true = c(N.given, N.given / M, numcap.given, N.Male.given, N.Male.given / N.given, phi.given, p0.given, sigmam.given, sigmaf.given)

prob.quantiles = c(0.025, 0.5, 0.975)  # for credible limits
post.stats = cbind(apply(post,2,mean, na.rm = T), apply(post,2,sd, na.rm = T), 
                   t(apply(post, 2, quantile, probs = prob.quantiles, na.rm = T,type = 1)))

post.stats2 = t(cbind(apply(rbind(par_true, post),2,sim.stats, hpdprob = 0.95)))
dimnames(post.stats2)[2] = list(c('bias', 'relative_bias_percentage', 'RMSE', 'HPD_lower', 'HPD_upper'))
prob.names = paste(as.character(100*prob.quantiles), '%', sep = '')
dimnames(post.stats)[2] = list(c('Mean.Chain', 'SD.Chain', prob.names))

mcse.mean.vec = mcse.sd.vec = rep(0, dim(post)[2])
mcse.lq.vec = mcse.uq.vec = mcse.med.vec = rep(0, dim(post)[2])
for (i in 1:dim(post)[2])
{
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
# Print summaries of posterior
# Print results of classical and Bayesian analyses
# cat ('\n', '\n', 'Bayesian estimates of model parameters', '\n', '\n', sep=',', file = paste0(folderName, '/EstimatesOfDerivedParam.txt', sep=""), append = TRUE)
# cat(post.stats, sep=',', file = paste0(folderName, '/EstimatesOfDerivedParam.txt', sep=""), append = TRUE)
# post.stats = cbind(post.stats, rep('', nrow(post.stats)))
HB_estimates = sim_estimates = MCSE_estimates = c('', '', '', '', '')
dim.names1 = c('Mean.Chain', 'SD.Chain', prob.names)
dim.names2 = c('bias', 'relative_bias_percentage', 'RMSE', 'HPD_lower', 'HPD_upper')
dimnames(mcse.mat) = dimnames(post.stats)

# information = as.data.frame(c(M, N.given, N.Male.given, Ninexpo, numl , numr, IDfixed, known, rad.given, phi.given, 
#                             sigmam.given, sigmaf.given, K, J, ndraws, burnin))
# dimnames(information)=list(c('M', 'N.given', 'N.Male.given','Ninexpo', 'numl', 'numr', 'IDfixed', 'known', 'rad.given', 'phi.given', 
#                              'sigmam.given', 'sigmaf.given', 'K', 'J', 'ndraws', 'burnin'),
#                            c('info'))


information = as.data.frame(c(M, N.given, N.Male.given, numl , numr, IDfixed, known, phi.given, p0.given, 
                              sigmam.given, sigmaf.given, K, J, ndraws, burnin, buffer))
dimnames(information) = list(c('M', 'N.given', 'N.male.Given', 'numl', 'numr', 'IDfixed', 'known', 'phi.given', 'p0.given', 
                               'sigmam.given', 'sigmaf.given', 'K', 'J', 'ndraws', 'burnin', 'buffer'),
                             c('info'))

# write.csv(sim.settings, file = paste0(folderName, '/simulation settings.csv', sep=""), quote = F,row.names = T)

out.final = rbind(t(out), HB_estimates, dim.names1, post.stats, sim_estimates, dim.names2, post.stats2, MCSE_estimates, dim.names1, mcse.mat)
print(post.stats)
print(post.stats2)
print(mcse.mat)

# write.csv(out.final, file = paste0(folderName, '/EstimatesOfDerivedParam.csv', sep=""), quote = F,row.names = T)
write.csv(rbind(cbind(out.final, matrix('',nrow(out.final), nrow(information) - ncol(out.final))), 
                dimnames(information)[[1]], t(information)), 
          file = paste0(folderName, '/EstimatesOfDerivedParam.csv', sep = ""), quote = F,row.names = T)


#============================================================================================================
post = read.csv(paste0(folderName, '/markovchain.txt', sep = ""), sep = ",", header = T)

post.mcmc = as.mcmc(post)

fname6 = paste0(folderName, "/traceplots_N.jpeg", sep = "")
jpeg(fname6, width = 1000, height = 1000, units = "px", pointsize = 12, quality = 100)
traceplot(post.mcmc[,1], xlab = "Iterations", ylab = "N", main = "Traceplot of N",
          cex.main = 2, cex.lab = 2, cex.axis = 2)
dev.off()

fname6 = paste0(folderName, "/traceplots_psi.jpeg", sep = "")
jpeg(fname6, width = 1000, height = 1000, units = "px", pointsize = 12, quality = 100)
traceplot(post.mcmc[,2], xlab = "Iterations", ylab = "psi", main = "Traceplot of psi",
          cex.main = 2, cex.lab = 2, cex.axis = 2)
dev.off()

fname6 = paste0(folderName, "/traceplots_numcap.jpeg", sep = "")
jpeg(fname6, width = 1000, height = 1000, units = "px", pointsize = 12, quality = 100)
traceplot(post.mcmc[,3], xlab = "Iterations", ylab = "numcap", main = "Traceplot of numcap",
          cex.main = 2, cex.lab = 2, cex.axis = 2)
dev.off()

fname6 = paste0(folderName, "/traceplots_N.Male.jpeg", sep = "")
jpeg(fname6, width = 1000, height = 1000, units = "px", pointsize = 12, quality = 100)
traceplot(post.mcmc[,4], xlab = "Iterations", ylab = "N.Males", main = "Traceplot of N.Males",
          cex.main = 2, cex.lab = 2, cex.axis = 2)
dev.off()

fname6 = paste0(folderName, "/traceplots_theta.jpeg", sep = "")
jpeg(fname6, width = 1000, height = 1000, units = "px", pointsize = 12, quality = 100)
traceplot(post.mcmc[,5], xlab = "Iterations", ylab = "theta", main = "Traceplot of theta",
          cex.main = 2, cex.lab = 2, cex.axis = 2)
dev.off()

fname6 = paste0(folderName, "/traceplots_phi.jpeg", sep = "")
jpeg(fname6, width = 1000, height = 1000, units = "px", pointsize = 12, quality = 100)
traceplot(post.mcmc[,6], xlab = "Iterations", ylab = "phi", main = "Traceplot of phi",
          cex.main = 2, cex.lab = 2, cex.axis = 2)
dev.off()

fname6 = paste0(folderName, "/traceplots_p0.jpeg", sep = "")
jpeg(fname6, width = 1000, height = 1000, units = "px", pointsize = 12, quality = 100)
traceplot(post.mcmc[,7], xlab = "Iterations", ylab = "p0", main = "Traceplot of p0",
          cex.main = 2, cex.lab = 2, cex.axis = 2)
dev.off()

fname6 = paste0(folderName, "/traceplots_sigmam.jpeg", sep = "")
jpeg(fname6, width = 1000, height = 1000, units = "px", pointsize = 12, quality = 100)
traceplot(post.mcmc[,8], xlab = "Iterations", ylab = "sigma.male", main = "Traceplot of sigma.male",
          cex.main = 2, cex.lab = 2, cex.axis = 2)
dev.off()

fname6 = paste0(folderName, "/traceplots_sigmaf.jpeg", sep = "")
jpeg(fname6, width = 1000, height = 1000, units = "px", pointsize = 12, quality = 100)
traceplot(post.mcmc[,9], xlab = "Iterations", ylab = "sigma.female", main = "Traceplot of sigma.female",
          cex.main = 2, cex.lab = 2, cex.axis = 2)
dev.off()

#======================================================
post = read.csv(paste0(folderName, '/markovchain.txt', sep = ""), sep = ",", header = T)
post = post[(burnin + 1):ndraws,]
seq1 = seq(1, ndraws - burnin, by = 10)
# post.mcmc = as.mcmc(post)
N1 = post[seq1,1]; N1 = N1[!is.na(N1)]
psi1 = post[seq1,2]; psi1 = psi1[!is.na(psi1)]
numcap1 = post[seq1,3]; numcap1 = numcap1[!is.na(numcap1)]
N.Male1 = post[seq1,4]; N.Male1 = N.Male1[!is.na(N.Male1)]
theta1 = post[seq1,5]; theta1 = theta1[!is.na(theta1)]
phi1 = post[seq1,6]; phi1 = phi1[!is.na(phi1)]
p01 = post[seq1,7]; p01 = p01[!is.na(p01)]
sigmam1 = post[seq1,8]; sigmam1 = sigmam1[!is.na(sigmam1)]
sigmaf1 = post[seq1,9]; sigmaf1 = sigmaf1[!is.na(sigmaf1)]

df = data.frame(N1, psi1, N.Male1, theta1, numcap1, phi1, p01, sigmam1, sigmaf1)

fname6 = paste0(folderName, "/Scatter Plots_", IDfixed,known, ".jpeg", sep = "")
jpeg(fname6, width = 1000, height = 1000, units = "px", pointsize = 12, quality = 100)
# pairs(~N1+ psi1+ N.Male1+ theta1+ phi1+ sigmam1+ sigmaf1+ rad1, data = df, main="")
pairs(~N1 + psi1 + N.Male1 + theta1 + phi1 + p01 + sigmam1 + sigmaf1, data = df, main = "")

dev.off()
#========================================

if(known!='ALL')
{
  which.left.obsd = data$which.left.obsd
  which.right.obsd = data$which.right.obsd
  known.ID.indx = data$known.ID.indx
}

if(known =='ALL')
{
  which.left.obsd = 1:numl
  which.right.obsd = 1:numl
  known.ID.indx = 1:numl
}

tdf = cbind(1:K, trap.locations, matrix(1, K, J))
dimnames(tdf)[2] = list(c("", "Easting", "Northing", 1:J))
write.csv(tdf, 
          file = paste0(folderName, '/trap_deployment_file_right.csv', sep = ""), 
          quote = F,row.names = F)

edf.left = NULL
for(i in 1:dim(left.obs)[1]){ for(k in 1:K){ for(t in 1:J){
  if(left.obs[i,k,t] == 1 & tdf[k, t] == 1){ edf.left = rbind(edf.left, c(1, i, t, k))} # ; print(c(i, k, t))}
  
}}}
dimnames(edf.left)[2] = list(c("session", "individual", "occasion", "trap"))
write.csv(edf.left, 
          file = paste0(folderName, '/encounter_data_file_left.csv', sep = ""), 
          quote = F,row.names = F)
# left3d = SCR23darray(edf.left, tdf) # nind x ntraps x noccassions.

edf.right = NULL
for(i in 1:dim(right.obs)[1]){ for(k in 1:K){ for(t in 1:J){
  if(right.obs[i,k,t] == 1 & tdf[k, t] == 1){ edf.right = rbind(edf.right, c(1, i, t, k))} # ; print(c(i, k, t))}
  
}}}
dimnames(edf.right)[2] = list(c("session", "individual", "occasion", "trap"))
write.csv(edf.right, 
          file = paste0(folderName, '/encounter_data_file_right.csv', sep = ""), 
          quote = F,row.names = F)

#====================================================
save.image(paste0(folderName, '/savingRimage.RData', sep = ""))

# } # loop for j

#=================================================

