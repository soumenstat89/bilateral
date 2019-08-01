sim.gender.partial.data = function(N, N.Male, phi, p0, sigmam, sigmaf, nrow_trap, ncol_trap, J, xlim, ylim, buffer)
{
 
  
# N = 100; N.Male = 40; phi = 0.8; p0 = 0.07; sigmam = 2; sigmaf = 1; 
# nrow_trap = 10; ncol_trap = 16; J = 30; 
# xlim = c(0,29); ylim = c(0,35); buffer = 10; 
  #======================================================================
  
  e2dist = function(x,y) # x is matrix mx2, y is matrix nx2
  {
    ind = sort(rep(1:nrow(y), nrow(x))) # [rep(1, nrow(x)), rep(2, nrow(x)), ..., rep(nrow(y), nrow(x))]'
    dvec = sqrt((x[,1] - y[ind,1]) ^ 2 + (x[,2] - y[ind,2]) ^ 2) # [d(1,1), d(2,1), ...,d(nrow(x), 1), ..., d(1,nrow(y)), d(2,nrow(y)),...,d(nrow(x), nrow(y))]
    dmat = matrix(dvec, nrow = nrow(x), ncol = nrow(y), byrow = F) 
    return(dmat)
  }
  
  
  pimatfn = function(logitp0, act.center, trap.locations, sex, logsigmam, logsigmaf)
  {
    logsigma = rep(logsigmaf, nrow(act.center))
    logsigma[sex == 1] = logsigmam
    pimat = expit(logitp0) * exp(- e2dist(act.center, trap.locations) ^ 2 / (2 * (exp(logsigma) ^ 2))) # Division is always column wise
    # pimat = 1- exp(- expit(logitp0) * exp(- e2dist(act.center, trap.locations) ^ 2 / (2 * (exp(logsigma) ^ 2))))
    return(pimat)
  }
  
  act.center = cbind(runif(N, xlim[1], xlim[2]), runif(N, ylim[1], ylim[2]))
  act.x = act.center[,1]
  act.y = act.center[,2]
  sgrid.area = (xlim[2] - xlim[1])*(ylim[2] - ylim[1])
  trap.locations = as.matrix(expand.grid(seq(xlim[1] + buffer, xlim[2] - buffer, length.out = nrow_trap), 
                                         seq(ylim[1] + buffer, ylim[2] - buffer, length.out = ncol_trap)))
  traploc.x = trap.locations[,1]
  traploc.y = trap.locations[,2]
  K = nrow_trap * ncol_trap

  
  # Getting the captures
  ind=sample(1:N, N.Male, replace=F)
  sex.true=rep(0,N)
  sex.true[ind]=1
 
  pimat = pimatfn(logit(p0), act.center, trap.locations, sex.true, log(sigmam), log(sigmaf))
  
  left = right = array(0,dim = c(N,K,J))
  
  # Now loop over individuals and traps and fill-in the encounters
  for (i in 1:N){ for (k in 1:K){ for (t in 1:J){
    u = runif(1,0,1)
    if (u <= pimat[i,k])
    {
      left[i,k,t] = rbinom(1,1,phi)
      right[i,k,t] = rbinom(1,1,phi)
    } 
  }}}
  
  #   Put known IDfixed guys first.
  known.ID = apply(left+right, 1, function(amat){sum(amat>1)>0}) # individuals with at least one simultaneous capture (l,r)=(1,1)
  known.ID.indx = c(1:N)[known.ID] # which(known.ID, arr.ind = T)
  IDfixed = sum(known.ID)
  
  kp.left = apply(left,1,sum) > 0 # captured left guys
  kp.right = apply(right,1,sum) > 0 # captured right guys
  
  # Now we have a simulated data set
  if (IDfixed == 0) known = 'none'
  if (IDfixed>0)
  {
    known = 'some'
    logic1 = all.equal(kp.left, known.ID)
    logic2 = all.equal(kp.right, known.ID)
    if (logic1 == T & logic2 == T) known = 'ALL'
  }

  original.left = left
  original.right = right
  
  # Determine which ones were captured OR have known ID (which includes some 0 guys)
  # The IDfixed guys are given first in the output data set followed by the captured guys
  if (known != 'ALL')
  {
    which.left.obsd = c(1:N)[kp.left] # which(kp.left, arr.ind = T)
    left = left[c(known.ID.indx, setdiff(which.left.obsd, known.ID.indx)),,]
    numl = sum(kp.left)
    sexl = sex.true[c(known.ID.indx, setdiff(which.left.obsd, known.ID.indx))] # first IDfixed number of elements of sexl are known
    which.left.obsd = c(known.ID.indx, setdiff(which.left.obsd, known.ID.indx))
    
    which.right.obsd = c(1:N)[kp.right] # which(kp.right, arr.ind = T)
    right = right[c(known.ID.indx, setdiff(which.right.obsd, known.ID.indx)),,]
    numr = sum(kp.right)
    sexr = sex.true[c(known.ID.indx, setdiff(which.right.obsd, known.ID.indx))] # first IDfixed number of elements of sexr are known
    which.right.obsd = c(known.ID.indx, setdiff(which.right.obsd, known.ID.indx))
    
  }
  
  if (known == 'ALL')
  {
    cap = apply(left+right,1,sum)>0
    left = left[cap,,]
    right = right[cap,,]
    sexl = sexr = sex.true[cap]
    numl = numr = sum(cap)
  }
  
  if (known!='ALL')
  {
    out = list(left = left, right = right, IDfixed = IDfixed, known = known, 
             sgrid.area = sgrid.area, trap.locations = trap.locations, 
             K = K, J = J, xlim = xlim, ylim = ylim, buffer = buffer, act.center = act.center,
             which.left.obsd = which.left.obsd, which.right.obsd = which.right.obsd, known.ID.indx = known.ID.indx,
             original.left = original.left, original.right = original.right, kp.left = kp.left, kp.right = kp.right,
             N = N, N.Male = N.Male, 
             phi = phi, p0 = p0, sigmam = sigmam, sigmaf = sigmaf, 
             sexl = sexl, sexr = sexr, sex.true = sex.true, 
             numl = numl, numr = numr)
  }
  if (known=='ALL')
  {
    out = list(left = left, right = right, IDfixed = IDfixed, known = known, 
             sgrid.area = sgrid.area, trap.locations = trap.locations, 
             K = K, J = J, xlim = xlim, ylim = ylim, buffer = buffer, act.center = act.center,
             original.left = original.left, original.right = original.right, kp.left = kp.left, kp.right = kp.right,
             N = N, N.Male = N.Male, 
             phi = phi, p0 = p0, sigmam = sigmam, sigmaf = sigmaf, 
             sexl = sexl, sexr = sexr, sex.true = sex.true, 
             numl = numl, numr = numr)
  }
  
  return(out)
} # end of sim.gender.data.partial

