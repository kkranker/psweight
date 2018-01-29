
mata:

mata set matastrict on

class gmatch
{
  private:
    real colvector   T, W, sel0, sel1, Y0
    real matrix      X
    string scalar    treatvar, depvar, wgtvar
    string rowvector varlist
    real scalar      N1, N0, N1_raw, N0_raw, havewgts
    void             calcmeans(), calcvariances(), calccovariances(), logit_eval()
    real rowvector   olsbeta(), diagvariance(), logitbeta()
    real rowvector   means1, means0, meansP, variances0, variances1, variancesP, covariances0, covariances1, covariancesP

  public:
    void             new(), set(), set_W(), set_Y()
    real rowvector   diff(), stddiff(), varratio()
    real rowvector   stddiff, varratio
    real scalar      mean_asd, min_asd, max_asd
    real colvector   logitpscores(), ipwweights()
}

// The following functions read data into the instance of the class
// set_W() needs to be called after set_T()
void gmatch::new()
{
  /* */  "New instance of gmatch() created"
  /* */  "gmatch::new() doesn't do anything"
  /* */  "T is " + strofreal(rows(this.T)) + " by " + strofreal(cols(this.T))

}

void gmatch::set(string scalar treatvar, string scalar varlist, string scalar tousevar, | string scalar wgtvar)
{
  // Define treatment dummy
  this.treatvar = treatvar
  st_view(this.T, ., treatvar, tousevar)
  this.sel0 = selectindex(!T)
  this.sel1 = selectindex( T)
  /* */  "T is " + strofreal(rows(this.T)) + " by " + strofreal(cols(this.T))

  // Save number of observations in private macros
  this.N0_raw = rows(this.sel0)
  this.N1_raw = rows(this.sel1)
  if (!this.N0_raw | !this.N1_raw) _error("At least one treatment and control observation required.")
  /* */  strofreal(this.N0_raw) + " control obs"
  /* */  strofreal(this.N1_raw) + " treatment obs"

  // Define covariates
  this.varlist = tokens(varlist)
  st_view(this.X, ., this.varlist, tousevar)
  /* */  "X is " + strofreal(rows(this.X)) + " by " + strofreal(cols(this.X))

  // Define weights
  if (args()>3) {
    this.havewgts = 1
    this.wgtvar = wgtvar
    st_view(this.W, ., wgtvar, tousevar)
    // This code assumes weights are **already** normalized. Here's code to normalize if I change my mind:
    //    this.W = this.W :/ (rows(this.W) / quadcolsum(this.W))
    this.N0 = quadcolsum(this.W[this.sel0])
    this.N1 = quadcolsum(this.W[this.sel1])
    if (!this.N0 | !this.N1) _error("Sum of weights is 0 in the treatment or control group.")
    /* */ "W is " + strofreal(rows(this.T)) + " by " + strofreal(cols(this.T))
    /* */ strofreal(this.N0) + " weighted control obs"
    /* */ strofreal(this.N1) + " weighted treatment obs"
  }
  else {
    this.havewgts = 0
    this.wgtvar = ""
    this.W = 1
    this.N0 = this.N0_raw
    this.N1 = this.N1_raw
    /* */ "Data are unweighted"
  }
}

void gmatch::set_Y(string scalar depvar, string scalar tousevar)
{
  real colvector Y
  this.depvar = depvar
  Y=.
  st_view(Y, ., this.depvar, tousevar)
  st_select(this.Y0, Y, !this.T)
  /* */  "Y0 is " + strofreal(rows(this.Y0)) + " by " + strofreal(cols(this.Y0))
}

// This function calculates the means for the T and C groups,
// and saves the results in private variables
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
  /* */ "Control group means:"; this.means0
  /* */ "Treatment group means:"; this.means1
  /* */ "Pooled means:"; this.meansP
}

// This function calculates the difference in means between the T and C groups
real rowvector gmatch::diff()
{
  if (!length(this.means1)) this.calcmeans()
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

// This function calculates the variances for the T and C group,
// and saves the results in private variables
void gmatch::calcvariances()
{
  if (!length(this.means1)) this.calcmeans()
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
  /* */ "Control group variances:"; this.variances0
  /* */ "Treatment group variances:"; this.variances1
  /* */ "Pooled variances:"; this.variancesP
}

// This function calculates the variances for the T and C group,
// and saves the results in private variables
void gmatch::calccovariances()
{
  if (this.havewgts) {
    this.covariances0 = quadvariance(this.X[this.sel0, .], this.W[this.sel0])
    this.covariances1 = quadvariance(this.X[this.sel1, .], this.W[this.sel1])
    this.covariancesP = quadvariance(this.X              , this.W           )
  }
  else {
    this.covariances0 = quadvariance(this.X[this.sel0, .])
    this.covariances1 = quadvariance(this.X[this.sel1, .])
    this.covariancesP = quadvariance(this.X              )
  }
  /* */ "Control group covariances:"; this.covariances0
  /* */ "Treatment group covariances:"; this.covariances1
  /* */ "Pooled covariances:"; this.covariancesP
  this.variances0 = diagonal(this.covariances0)'
  this.variances1 = diagonal(this.covariances1)'
  this.variancesP = diagonal(this.covariancesP)'
}

// This function calculates standardized differences in means between the T and C groups
// The first argument is optional, and tells the function which variance to use in the denominator
//    If first argument is 0, it uses the control groups' variances
//    If first argument is 1, it uses the treatment groups' variances (this is the default)
//    If first argument is 2, it uses the pooled variances
// The second argument is optional, and tells the function whether to calculate the full covariance matrix.
//    By default, the function only calculates the variances (the diagonal of the covariance matrix, without the off-diagonal elements)
// The standardized differences and mean/min/max absolute standardized differences are all saved in class variables
real rowvector gmatch::stddiff(| real scalar denominator, real scalar cov)
{
  if (args()<1) denominator=1
  if (args()<2) cov=0

  if (cov & !length(this.covariances1)) this.calccovariances()
  else if (!length(this.variances1))    this.calcvariances()

  if      (denominator==0) this.stddiff = (this.diff() :/ this.variances0)
  else if (denominator==1) this.stddiff = (this.diff() :/ this.variances1)
  else if (denominator==2) this.stddiff = (this.diff() :/ this.variancesP)
  else _error(strofreal(denominator)+ " is an invalid arguement for gmatch::stddiff()")

  this.mean_asd = mean(abs(this.stddiff'))
  this.min_asd  =  min(abs(this.stddiff ))
  this.max_asd  =  max(abs(this.stddiff ))

  return(this.stddiff)
}

// This function calculates ratio of variances between the T and C groups
// The first argument is optional, and tells the function whether to calculate the full covariance matrix.
//    By default, the function only calculates the variances (the diagonal of the covariance matrix, without the off-diagonal elements)
// The ratio of variances are all saved in class variables
real rowvector gmatch::varratio(| real scalar cov)
{
  if (args()<1) cov=0

  if (cov & !length(this.covariances1)) this.calccovariances()
  else if (!length(this.variances1))    this.calcvariances()

  this.varratio = this.variances1 :/ this.variances0
  return(this.varratio)
}

// Define function to calculate coefficients for an OLS regression model
// A contant term is included in the regression, but its coefficient
// is dropped from output (so number of coefficients matches number of columns)
real rowvector gmatch::olsbeta(real colvector y, real matrix X, | real colvector w)
{
  real colvector beta
  real scalar C
  real matrix XX, Xy
  if (args()<3) w=1
  C = cols(X)
  XX = quadcross(X, 1, w, X, 1)
  Xy = quadcross(X, 1, w, y, 0)
  beta  = invsym(XX,++C)*Xy
  return(beta[1..--C]')
}

// function that returns IPW weights (based on logit model that regresses T on X)
real colvector gmatch::ipwweights(| real scalar ate)
{
  real colvector pm, pscore, ipwwgt
  if (args()<1) ate=0

  pm      = 1 :- (!this.T)
  pscore  = this.logitpscores()

  if (ate) ipwwgt = (pm :* (1:/pscore)) :+ (!pm :* (1:/(1:-pscore)))
  else     ipwwgt = pm :+ (!pm :* (pscore:/(1:-pscore)))

  // normalizelize the weights
  ipwwgt = ipwwgt :/ mean(ipwwgt,this.W)

  /* */ "pm[1..20,.]"; pm[1..20,.]
  /* */ "pscore[1..20,.]"; pscore[1..20,.]
  /* */ "ipwwgt[1..20,.]"; ipwwgt[1..20,.]
  /* */ "POmean: "; mean(this.Y0, ipwwgt[this.sel0]:*this.W)

  return(ipwwgt)
}

// function that returns IPW weights (based on logit model that regresses T on X)
real colvector gmatch::logitpscores()
{
  real matrix    X1
  real rowvector beta
  real colvector pscore

  X1     = (this.X,J(rows(this.X),1,1))
  beta   = this.logitbeta(this.T, X1)
  pscore = invlogit(X1*beta')
  return(pscore)
}

// Define function to calculate coefficients for a logit regression model
// A contant term is included in the regression and its coefficient is included in the vector of betas
// Source: https://www.stata.com/statalist/archive/2010-10/msg01188.html
//    From   jpitblado@stata.com (Jeff Pitblado, StataCorp LP)
//    To   statalist@hsphsun2.harvard.edu
//    Subject   Re: st: pointing to a class member function (likelihood) with optimize() in mata
//    Date   Thu, 28 Oct 2010 10:19:51 -0500
real rowvector gmatch::logitbeta(real colvector y, real matrix X)
{
  real rowvector beta
  transmorphic S

  S = optimize_init()
  optimize_init_argument(S, 1, y)
  optimize_init_argument(S, 2, X)
  optimize_init_evaluator(S, &logit_eval())
  optimize_init_evaluatortype(S, "d2")
  optimize_init_params(S, J(1,cols(X),0))

  beta = optimize(S)

  /* */ "Logit model coefficients & std. errors"; (optimize_result_params(S) \ sqrt(diagonal(optimize_result_V_oim(S)))')'

  return(beta)
}

void logit_eval( real scalar    todo,
                 real rowvector beta,
                 real colvector y,
                 real matrix    X,
                 real scalar    lnf,
                 real rowvector g,
                 real matrix    H)
{
  real colvector  pm
  real colvector  xb
  real colvector  lj
  real colvector  dllj
  real colvector  d2llj

  pm      = 2*(y :!= 0) :- 1
  xb      = X*beta'

  lj      = invlogit(pm:*xb)
  if (any(lj :== 0)) {
          lnf = .
          return
  }
  lnf = quadcolsum(ln(lj))
  if (todo == 0) return

  dllj    = pm :* invlogit(-pm:*xb)
  if (missing(dllj)) {
          lnf = .
          return
  }
  g       = quadcross(dllj, X)
  if (todo == 1) return

  d2llj   = abs(dllj) :* lj
  if (missing(d2llj)) {
          lnf = .
          return
  }
  H       = - quadcross(X, d2llj, X)
}

end
