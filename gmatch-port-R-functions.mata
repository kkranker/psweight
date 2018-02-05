// put into class definitions
//    real rowvector cbpsloss(), cbpsloss2(), cbpslossOID()
//    void           cbpseval_port()
//
// add to CDPS_evaleval
//   if(fctn=="cbpseval_port") M.cbpseval_port(todo,beta,est,fctn,denominator,lnf,g,H)
//   else

mata:

real rowvector gmatch::cbpsloss2()
{
//          abs( (this.W :* X) * invsym(quadcross(this.X,this.W,this.X)) * (this.W :* X)' )
//  return( abs( (this.W :* X) * invsym(quadcross(this.X,this.W,this.X)) * (this.W :* X)' ) )
          abs( this.W' * X * invsym(quadcross(this.X,this.W,this.X)) * X' * this.W )
  return( abs( this.W' * X * invsym(quadcross(this.X,this.W,this.X)) * X' * this.W )  )
}

real rowvector gmatch::cbpsloss(real colvector probs, real colvector w, | string scalar est)
{
  real rowvector loss1
  if (!length(this.covariancesP)) this.calccovariances()
  loss1 = abs( w' * (this.W :* X) * invsym(quadcross(this.X,this.W,this.X)) * (this.W :* X)' * w )
  return(loss1)
}

real rowvector gmatch::cbpslossOID(real colvector probs, real colvector w, | string scalar est)
{
  real colvector w_del, gbar, wx0, wx1, wxP, loss
  real matrix V, invV
  if (args()<3) est="ate"

  // Generate the vector of mean imbalance by weights.
  w_del = 1/this.N :* (this.W :* this.X)' * (w)

  // Generate g-bar, as in the paper.
  gbar = 1/this.N :* (this.W :* this.X)' * (this.T-probs)
  gbar = (gbar \ w_del )

  // Generate the covariance matrix used in the GMM estimate.
  // Was for the initial version that calculates the analytic variances.
  if (strlower(est)=="atet") {
    wx0   = this.W:^.5:*this.X:*((1:-probs):*probs):^.5
    wx1   = this.W:^.5:*this.X:*(probs:/(1:-probs)):^.5
    wxP   = this.W:^.5:*this.X:*(probs):^.5
  }
  else if (strlower(est)=="ate") {
    wx0   = this.W:^.5:*this.X:*((1:-probs):*probs):^.5
    wx1   = this.W:^.5:*this.X:*(probs:*(1:-probs)):^-.5
    wxP   = this.W:^.5:*this.X
  }
  else _error(est + " is not allowed.")

  if (strlower(est)=="atet") {
    V = ( (1/this.N) :* ((wx0' * wx0) , (wxP' * wxP)) :* (this.N / this.N1)  \
          (1/this.N) :* ((wxP' * wxP) :* (this.N/this.N1) , (wx1' * wx1) :* (this.N^2/this.N1^2) ) )
  }
  else {
    V = ( (1/this.N) :* ((wx0' * wx0) , (wxP' * wxP)) \
          (1/this.N) :* ((wxP' * wxP) , (wx1' * wx1)) )
  }

  invV = invsym(V)

  // Calculate the GMM loss
  loss = gbar' * invV * gbar
  return(loss)

}

// Port of the gmm.func()  function from CBPS.Binary.R (version 0.17)
void gmatch::cbpseval_port(real    /*scalar        */ todo,
                       real    /*rowvector     */ beta,
                       string  /*scalar        */ est,
                       string  /*scalar        */ fctn,
                       real    /*scalar        */ denominator,
                       real    /*colvector     */ lnf,
                       real    /*rowvector     */ g,
                       real    /*matrix        */ H)
{
    real colvector probs, w, w_del, gbar, wx0, wx1, wxP, loss
    real matrix V, invV
    // Designate sample size, number of treated and control observations,
    // theta, which are used to generate probabilities.
    // Trim probabilities, and generate weights.
    probs = this.logitpredict(this.X, beta)
    probs = this.trim(probs)
    w = this.logitweights(probs, est)
    // Generate the vector of mean imbalance by weights.
    w_del = 1/this.N :* (this.W :* this.X)' * (w)
    // Generate g-bar, as in the paper.
    gbar = 1/this.N :* (this.W :* this.X)' * (this.T-probs)
    gbar = (gbar \ w_del )
    // Generate the covariance matrix used in the GMM estimate.
    if (strlower(est)=="atet") {
      wx0   = this.W:^.5:*this.X:*((1:-probs):*probs):^.5
      wx1   = this.W:^.5:*this.X:*(probs:/(1:-probs)):^.5
      wxP   = this.W:^.5:*this.X:*(probs):^.5
    }
    else if (strlower(est)=="ate") {
      wx0   = this.W:^.5:*this.X:*((1:-probs):*probs):^.5
      wx1   = this.W:^.5:*this.X:*(probs:*(1:-probs)):^-.5
      wxP   = this.W:^.5:*this.X
    }
    else _error(est + " is not allowed.")
    if (strlower(est)=="atet") {
      V = ( (1/this.N) :* ((wx0' * wx0) , (wxP' * wxP)) :* (this.N / this.N1)  \
            (1/this.N) :* ((wxP' * wxP) :* (this.N/this.N1) , (wx1' * wx1) :* (this.N^2/this.N1^2) ) )
    }
    else {
      V = ( (1/this.N) :* ((wx0' * wx0) , (wxP' * wxP)) \
            (1/this.N) :* ((wxP' * wxP) , (wx1' * wx1)) )
    }

    invV = invsym(V)

    // Calculate the GMM loss
    lnf = gbar' * invV * gbar

}

end

