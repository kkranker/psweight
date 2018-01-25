
mata:

class gmatch
{
  private:
    real colvector   T, W, sel0, sel1, Y0
    real matrix      X
    string scalar    treatvar, depvar, wgtvar
    string rowvector varlist
    real scalar      N1, N0, N1_raw, N0_raw, havewgts
    void             calcmeans(), calcvariances()
    real rowvector   means1, means0, meansP, variances0, variances1, variancesP
    real rowvector   diagvariance()

  public:
    void             new(), set(), set_X(), set_T(), set_W(), set_Y()
    real rowvector   diff()
    real rowvector   stddiff()
}

// The following functions read data into the instance of the class
// set_W() needs to be called after set_T()
void gmatch::new()
{
  this.wgtvar = this.depvar = ""
  this.W = 1
  this.havewgts = 0
  this.Y0 = .
  /* */ "By default, data is unweighed and there is no depvar"
}

void gmatch::set_X(string scalar varlist, string scalar tousevar)
{
  this.varlist = tokens(varlist)
  st_view(this.X, ., this.varlist, tousevar)
  /* */  "X is " + strofreal(rows(this.X)) + " by " + strofreal(cols(this.X))
}

void gmatch::set_T(string scalar treatvar, string scalar tousevar)
{
  this.treatvar = treatvar
  st_view(this.T, ., treatvar, tousevar)
  this.sel0 = selectindex(!T)
  this.sel1 = selectindex( T)
  this.N0 = this.N0_raw = rows(this.sel0)
  this.N1 = this.N1_raw = rows(this.sel1)
  if (!this.N0 | !this.N1) _error("At least one treatment and control observation required.")
  /* */  "T is " + strofreal(rows(this.T)) + " by " + strofreal(cols(this.T))
  /* */  strofreal(this.N0) + " control obs"
  /* */  strofreal(this.N1) + " treatment obs"
}

void gmatch::set_W(string scalar wgtvar, string scalar tousevar, | real scalar normalize)
{
  this.havewgts = 1
  this.wgtvar = wgtvar
  st_view(this.W, ., wgtvar, tousevar)
  if (normalize) this.W = this.W :/ (rows(this.W) / quadcolsum(this.W))
  this.N0 = quadcolsum(this.W[this.sel0])
  this.N1 = quadcolsum(this.W[this.sel1])
  if (!this.N0 | !this.N1) _error("Sum of weights is 0 in the treatment or control group.")
  /* */ "W is " + strofreal(rows(this.T)) + " by " + strofreal(cols(this.T))
  /* */ strofreal(this.N0) + " weighted control obs"
  /* */ strofreal(this.N1) + " weighted treatment obs"
}

void gmatch::set_Y(string scalar depvar, string scalar tousevar)
{
  real colvector Y
  this.depvar = depvar
  st_view(Y, ., this.depvar, tousevar)
  st_select(this.Y0, Y, !this.T)
  /* */  "Y0 is " + strofreal(rows(this.Y0)) + " by " + strofreal(cols(this.Y0))
}

// This function calculates the mean for the T and C group, and fills them into the private class variables
void gmatch::calcmeans()
{
  if (this.havewgts) {
    this.means0 = mean(this.X[this.sel0, .], this.W[this.sel0])
    this.means1 = mean(this.X[this.sel1, .], this.W[this.sel1])
    this.meansP = mean(this.X              , this.W           )
  }
  else {
    this.means0 = mean(this.X[this.sel0, .])
    this.means1 = mean(this.X[this.sel1, .])
    this.meansP = mean(this.X              )
  }
  /* */ "Control   group means:"; this.means0
  /* */ "Treatment group means:"; this.means1
  /* */ "Pooled          means:"; this.meansP
}

// This function calculates the difference in means between the T and C groups
real rowvector gmatch::diff()
{
  if (!rows(this.means1)) this.calcmeans()
  /* */ "Difference in means: "
  return(this.means1 :- this.means0)
}

// Define function diagvariance(x, w) == diagonal(quadvariance(x, w))'
// This function can be a lot faster than quadvariance, especially when you have lots of columns.
// Optionally, you can provide weights and/or provide a rowvector with the column means.
// For testing, mreldif(diagvariance(X, w), diagonal(quadvariance(X, w))') should be small
real rowvector gmatch::diagvariance(real matrix X, | real colvector w, real rowvector Xmean)
{
  real rowvector v
  if (args()<2) w = 1
  if (args()<3) Xmean = mean(X, w)

  if (w==1) v = quadcolsum( (X:-Xmean):^2)     / (rows(X)-1)
  else      v = quadcolsum(((X:-Xmean):^2):*w) / (quadcolsum(w)-1)
  return(v)
}

// This function calculates the variancess for the T and C group, and fills them into the private class variables
void gmatch::calcvariances()
{
  if (!rows(this.means1)) this.calcmeans()
  if (this.havewgts) {
    this.variances0 = this.diagvariance(this.X[this.sel0, .], this.W[this.sel0], this.means0)
    this.variances1 = this.diagvariance(this.X[this.sel1, .], this.W[this.sel1], this.means1)
    this.variancesP = this.diagvariance(this.X              , this.W           , this.meansP)
  }
  else {
    this.variances0 = this.diagvariance(this.X[this.sel0, .], 1, this.means0)
    this.variances1 = this.diagvariance(this.X[this.sel1, .], 1, this.means1)
    this.variancesP = this.diagvariance(this.X              , 1, this.meansP)
  }
  "Control   group variances:"; this.variances0
  "Treatment group variances:"; this.variances1
  "Pooled          variances:"; this.variancesP
}

// This function calculates standardized differences in means between the T and C groups
// The function argument is optional.  
//    If missing, it uses the pooled 
real rowvector gmatch::stddiff(| real scalar denominator)
{
  if (!rows(this.variances1))   this.calcvariances()
  
  /* */ "Standardized differences: "
  if (!args())             return(this.diff() :/ this.variancesP)
  else if (denominator==0) return(this.diff() :/ this.variances0)
  else if (denominator==1) return(this.diff() :/ this.variances1)
  else _error(strofreal(denominator)+ " is an invalid arguement")
}


end
