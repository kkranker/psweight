clear all
cls

mata:

mata set matastrict on

class gmatch
{
  private:
    real colvector   T, W, sel1, sel0, Y0, W_orig
    real matrix      X
    string scalar    treatvar, depvars, wgtvar
    string rowvector varlist
    void             logit_eval(), calcmeans(), calcvariances(), calcN()
    real scalar      N1, N0, N, N1_raw, N0_raw, N_raw
    real rowvector   olsbeta(), diagvariance(), logitfit()
    real colvector   olspredict(), logitpredict(), logitweights()
    real rowvector   means1, means0, variances0, variances1, variancesP, variancesA

  // real matrix variances0, variances1, variancesP
  // void        calccovariances()

  public:
    void             new(), set(), set_W(), set_Y(), clone(), multweight()
    real rowvector   diff(), stddiff(), varratio(), prognosticdiff(), pomean()
    real scalar      mean_asd(), min_asd(), max_asd()
    real colvector   ipw()
    real matrix      balancetable()
}

// The following functions read data into the instance of the class
// set_W() needs to be called after set_T()
void gmatch::new()
{
  // /* */ "New instance of gmatch() created"
  // /* */ "gmatch::new() doesn't do anything"
  // /* */ "T is " + strofreal(rows(this.T)) + " by " + strofreal(cols(this.T))
  this.depvars = ""
}

void gmatch::clone(class gmatch scalar origcopy)
{
  this.N          =  origcopy.N
  this.N_raw      =  origcopy.N_raw
  this.N0         =  origcopy.N0
  this.N0_raw     =  origcopy.N0_raw
  this.N1         =  origcopy.N1
  this.N1_raw     =  origcopy.N1_raw
  this.T          =  origcopy.T
  this.W          =  origcopy.W
  this.W_orig     =  origcopy.W_orig
  this.X          =  origcopy.X
  this.Y0         =  origcopy.Y0
  this.depvars    =  origcopy.depvars
  this.sel0       =  origcopy.sel0
  this.sel1       =  origcopy.sel1
  this.treatvar   =  origcopy.treatvar
  this.varlist    =  origcopy.varlist
  this.wgtvar     =  origcopy.wgtvar
}

void gmatch::set(string scalar treatvar, string scalar varlist, string scalar tousevar, | string scalar wgtvar)
{
  // Define treatment dummy
  this.treatvar = treatvar
  st_view(this.T, ., treatvar, tousevar)
  // /* */  "T is " + strofreal(rows(this.T)) + " by " + strofreal(cols(this.T))

  // Define covariates
  this.varlist = tokens(varlist)
  st_view(this.X, ., this.varlist, tousevar)
  // /* */  "X is " + strofreal(rows(this.X)) + " by " + strofreal(cols(this.X))
  // /* */  "X contains" ; this.varlist

  // Define weights
  // This code assumes weights are **already** normalized. Here's code to normalize: this.W = this.W :/ (rows(this.W) / quadcolsum(this.W))
  if (args()>3) {
    this.wgtvar = wgtvar
    st_view(this.W, ., this.wgtvar, tousevar)
    st_view(this.W_orig, ., this.wgtvar, tousevar) // an extra copy of the weight variable that can only be set via this function. Useful for reweighting/matching situations.
    // /* */  "W is " + strofreal(rows(this.W)) + " by " + strofreal(cols(this.W))
  }
  else {
    this.W  = this.W_orig = J(rows(T),1,1)
    "Data are unweighted."
  }
  this.calcN()
}
  
void gmatch::calcN() {

  // Index to select observations in control and treatment groups  
  this.sel0 = selectindex(!this.T :& this.W)
  this.sel1 = selectindex( this.T :& this.W)
  
  // Save number of observations
  this.N0_raw = rows(this.sel0)
  this.N1_raw = rows(this.sel1)
  this.N_raw = this.N0_raw + this.N1_raw
  if (min((this.N0_raw,this.N1_raw)==0)) _error("At least one treatment and control observation required.")
  
  // Save weighted number of observations
  this.N0 = quadcolsum(this.W[this.sel0])
  this.N1 = quadcolsum(this.W[this.sel1])
  this.N = this.N0 + this.N1 
  if (min((this.N0    ,this.N1    )==0)) _error("Sum of weights is 0 in the treatment or control group.")
  strofreal(this.N0_raw) + " control obs   (sum of weights = " + strofreal(this.N0) + ")"
  strofreal(this.N1_raw) + " treatment obs (sum of weights = " + strofreal(this.N1) + ")"
  
  // these means/varinaces are saved internally in the class (to avoid computing them over and over).  
  // They need to be reset because we just reweighted the sample.
  // If I'm re-calcuating sample sizes, this is probably the case.  Set to missing here just to be safe.
  this.means0 = this.means1 = this.variances0 = this.variances1 = this.variancesP = this.variancesA = J(1,0,.)  
  
}

// Note: this function doesn't allow the class to touch the treatment group's outcome data
void gmatch::set_Y(string scalar depvarnames, string scalar tousevar)
{
  real colvector Y
  this.depvars = tokens(depvarnames)
  Y=.
  st_view(Y, ., this.depvars, tousevar)
  st_select(this.Y0, Y, !this.T)
  // /* */  "Y0 is " + strofreal(rows(this.Y0)) + " by " + strofreal(cols(this.Y0))
}

// multipy the original weights by something
void gmatch::multweight(real colvector newweight)
{
  // Define weights
  this.W = this.W_orig :* newweight
  
  // some new weights could be zero.  recalculate N and set means/variances to missing.
  this.calcN()

  // these means/varinaces are saved internally in the class (to avoid computing them over and over).  They need to be reset because we just reweighted the sample.
  this.means0 = this.means1 = this.variances0 = this.variances1 = this.variancesP = this.variancesA = J(1,0,.)  
  
}


// This function makes a balance table and prints it to the screen
// The argument is the same as their definition in stddiff() and varratio()
real matrix gmatch::balancetable(| real scalar denominator)
{
  real matrix table, sd
  if (args()<1) denominator=1

  if (!length(this.means1))     this.calcmeans()
	if (!length(this.variances1)) this.calcvariances()

  table = ( this.means1
          \ this.means0
          \ this.diff()
          \ this.stddiff(denominator)
          \ (denominator==0 ? sqrt(this.variances0) : (denominator==1 ? sqrt(this.variances1) : (denominator==2 ? sqrt(this.variancesP) : (denominator==3 ? sqrt(this.variancesA) : _error("denominator argument invalid")))))
          \ this.varratio())'

  // print to screen with labels
  (("Variable" \ varlist'),
  (( "means1",
     "means0",
     "diff()",
     "stddiff()",
     (denominator==0 ? "sd0" : (denominator==1 ? "sd1" : (denominator==2 ? "sdP" : (denominator==3 ? "sdA" : "")))),
     "varratio()")
     \ strofreal(table) ))

  return(table)
}

// This function calculates the means for the T and C groups
// These means are saved internally in the class (to avoid computing them over and over)
// Call this function whenever sample or weights change
void gmatch::calcmeans()
{
  this.means0 = mean(this.X[this.sel0, .], this.W[this.sel0])
  this.means1 = mean(this.X[this.sel1, .], this.W[this.sel1])
  // /* */ "Control group means:"  ; this.means0
  // /* */ "Treatment group means:"; this.means1
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
real rowvector gmatch::diagvariance(real matrix x, | real colvector w, real rowvector xmean)
{
  real rowvector v
  if (args()<2) w = 1
  if (args()<3) xmean = mean(x, w)

  if (all(w:==1)) v = quadcolsum( (x:-xmean):^2)     / (rows(x)-1)
  else            v = quadcolsum(((x:-xmean):^2):*w) / (quadcolsum(w)-1)
  return(v)
}

// This function calculates the variances for the T and C group,
// These variances are saved internally in the class (to avoid computing them over and over)
// Call this function whenever sample or weights change
void gmatch::calcvariances()
{
  if (!length(this.means1)) this.calcmeans()
  this.variances0 = this.diagvariance(this.X[this.sel0, .], this.W[this.sel0], this.means0)
  this.variances1 = this.diagvariance(this.X[this.sel1, .], this.W[this.sel1], this.means1)
  this.variancesP = this.diagvariance(this.X, this.W)
  this.variancesA = (this.variances0 :+ this.variances1) :/ 2
  // /* */ "Control group variances:"; this.variances0
  // /* */ "Treatment group variances:"; this.variances1
  // /* */ "Pooled variances:"; this.variancesP
  // /* */ "Average of variances from treatment and control groups"; this.variancesA
}

//  // This function calculates the variances for the T and C group,
//  // and saves the results in private variables
//  void gmatch::calccovariances()
//  {
//    if (all(this.W:==1)) {
//      this.covariances0 = quadvariance(this.X[this.sel0, .])
//      this.covariances1 = quadvariance(this.X[this.sel1, .])
//      this.covariancesP = quadvariance(this.X)
//    }
//    else {
//      this.covariances0 = quadvariance(this.X[this.sel0, .], this.W[this.sel0])
//      this.covariances1 = quadvariance(this.X[this.sel1, .], this.W[this.sel1])
//      this.covariancesP = quadvariance(this.X, this.W)
//    }
//    this.covariancesA = (this.covariances0 :+ cothis.variances1) :/ 2
//    // /* */ "Control group covariances:"; this.covariances0
//    // /* */ "Treatment group covariances:"; this.covariances1
//    // /* */ "Pooled covariances:"; this.covariancesP
//    // /* */ "Average of covariances from treatment and control groups"; this.covariancesA
//    this.variances0 = diagonal(this.covariances0)'
//    this.variances1 = diagonal(this.covariances1)'
//    this.variancesP = diagonal(this.covariancesP)'
//    this.variancesA = diagonal(this.covariancesA)'
//    // /* */ "Average of variances from treatment and control groups"; this.variancesA
//  }

// This function calculates standardized differences in means between the T and C groups
// The first argument is optional, and tells the function which variance to use in the denominator
//    = 0, it uses the control groups' variances
//    = 1, it uses the treatment groups' variances (this is the default)
//    = 2, it uses the pooled variances
//    = 3, it uses (control groups' variances + treatment groups' variances)/2  (the definition from Stata's tbalance command)
real rowvector gmatch::stddiff(| real scalar denominator)
{
  if (args()<1) denominator=1
  real rowvector stddiff
  if (!length(this.variances1)) this.calcvariances()
  if      (denominator==0) stddiff = (this.diff() :/ sqrt(this.variances0))
  else if (denominator==1) stddiff = (this.diff() :/ sqrt(this.variances1))
  else if (denominator==2) stddiff = (this.diff() :/ sqrt(this.variancesP))
  else if (denominator==3) stddiff = (this.diff() :/ sqrt(this.variancesA))
  else _error(strofreal(denominator)+ " is an invalid argument for gmatch::stddiff()")
  return(stddiff)
}

// functions to return mean/min/max absolute standardized differences
real scalar gmatch::mean_asd(| real scalar denominator) {
  if (args()<1) denominator=1
  return(mean(abs(this.stddiff()')))
}
real scalar gmatch::min_asd(| real scalar denominator)  {
  if (args()<1) denominator=1
  return(min(abs(this.stddiff())))
}
real scalar gmatch::max_asd(| real scalar denominator)  {
  if (args()<1) denominator=1
  return(max(abs(this.stddiff())))
}

// This function calculates ratio of variances between the T and C groups
real rowvector gmatch::varratio()
{
  if  (!length(this.variances1)) this.calcvariances()
  return((this.variances1 :/ this.variances0))
}


// function that returns difference in y_hat, where y_hat is generated using a
// OLS regression of y on X using the control group data
real rowvector gmatch::prognosticdiff()
{
  real rowvector beta, progdiff
  real colvector yhat
  real scalar yhat_bar_0, yhat_bar_1, c
  if (!length(this.depvars)) _error("Dependent variable is undefined.  Use gmatch::set_Y().")

  yhat = J(rows(this.X), cols(this.Y0), .)
  for (c=1; c<=cols(this.Y0); c++) {
    beta = this.olsbeta(this.Y0[.,c], this.X[this.sel0,.], this.W[this.sel0])
    yhat[.,c] = this.olspredict(this.X, beta, 1)
  }

  yhat_bar_0 = mean(yhat[this.sel0,.], this.W[this.sel0])
  yhat_bar_1 = mean(yhat[this.sel1,.], this.W[this.sel1])
  progdiff = yhat_bar_1 :- yhat_bar_0
  
  // /* */ "Control group mean prognostic score: "  ; yhat_bar_0
  // /* */ "Treatment group mean prognostic score: "; yhat_bar_1
  // /* */ "Difference:"; progdiff
  return(progdiff)
}


// Define function to calculate coefficients for an OLS regression model
// A contant term is included in the regression.
real rowvector gmatch::olsbeta(real matrix y, real matrix X, | real colvector w)
{
  real colvector beta
  real scalar C
  real matrix XX, Xy

  if (args()<3) w=1
  C = cols(X)
  XX = quadcross(X, 1, w, X, 1)
  Xy = quadcross(X, 1, w, y, 0)
  beta  = invsym(XX,++C)*Xy
  return(beta')
}

// Function that returns predicted values (e.g., propensity scores) if given the X's and betas, using the logit model functional form
// If addconst==1, the function assumes the last coefficient corresponds to the constant term, and X doesn't have a constant term
// Warning: this function doesn't check the conformability; I rely on Stata to produce errors with invalid arguments
real colvector gmatch::olspredict(real matrix X, real rowvector beta, | real scalar addconst)
{
  real colvector yhat
  if (args()<3) addconst = 1
  if (addconst) yhat = (X*beta[1..(cols(beta)-1)]') :+ beta[cols(beta)]
  else          yhat =  X*beta'
  return(yhat)
}


// function that returns IPW weights
//    est corresponds to the options in gmatch::logitweights()
//    est = "ate"  computes weights for average treatment effect (the default)
//        = "atet" computes weights for average treatment effect on the treated
//        = "ateu" computes weights for average treatment effect on the untreated
real colvector gmatch::ipw(string scalar est)
{
  real rowvector beta
  real colvector pscore, ipwwgt
  if (args()<1) est="ate"
  beta   = this.logitfit(this.T, this.X, this.W)
  pscore = this.logitpredict(this.X, beta, 1)
  ipwwgt = this.logitweights(pscore, est)
  return(ipwwgt)
}


// function that returns (weighted) mean of the dependent variable(s)
real rowvector gmatch::pomean()
{
  if (this.depvars=="") _error("dependent variable not defined. use gmatch::set_Y()")
  return( mean(this.Y0, this.W[this.sel0]) )
}

// Function that returns predicted values (e.g., propensity scores) if given the X's and betas, using the logit model functional form
// If addconst==1, the function assumes the last coefficient corresponds to the constant term, and X doesn't have a constant term
// Warning: this function doesn't check the conformability; I assume Stata will produce an error with invalid arguments
real colvector gmatch::logitpredict(real matrix X, real rowvector beta, | real scalar addconst)
{
  real colvector pr
  if (args()<3) addconst = 1
  if (addconst) pr = invlogit((X*beta[1..(cols(beta)-1)]') :+ beta[cols(beta)])
  else          pr = invlogit( X*beta')
  return(pr)
}

// This turns a vector of pscores into IPW weights. this assumes a logit setup.
// Formulas match the normalized weights in Stata's teffects IPW command
//    pscore is a vector of propensity scores
//    est = "ate"  computes weights for average treatment effect (the default)
//        = "atet" computes weights for average treatment effect on the treated
//        = "ateu" computes weights for average treatment effect on the untreated
real colvector gmatch::logitweights(real colvector pscore, | string scalar est)
{
  real colvector pm, ipwwgt, minmax
  if (args()<2) est="ate"

  minmax = minmax(pscore)
  if (minmax[1,1]<=0 :| minmax[1,2]>=1) _error("Propensity scores need to be greater than 0 and less than 1.")
  if (minmax[1,1]<=0.03 & (est=="ate" | est=="ateu")) errprintf("Warning: minimum propensity score is %12.0g \n", minmax[1,1])
  if (minmax[1,2]>=0.97 & (est=="ate" | est=="atet")) errprintf("Warning: maximum propensity score is %12.0g \n", minmax[1,2])

  pm = 1 :- (!this.T)
  if      (est=="ate")   ipwwgt = (pm :* (1:/pscore)) :+ (!pm :* (1:/(1:-pscore)))
  else if (est=="atet")  ipwwgt =  pm :+ (!pm :* (pscore:/(1:-pscore)))
  else if (est=="ateu")  ipwwgt = !pm :+ ( pm :* ((1:-pscore):/pscore))
  else _error(est + " is an invalid argument for gmatch::logitweights()")

  // normalize the weights to have mean 1 in each group
  if (est=="ate" | est=="atet") ipwwgt[this.sel0] = ipwwgt[this.sel0] :/ mean(ipwwgt[this.sel0], this.W[this.sel0])
  if (est=="ate" | est=="ateu") ipwwgt[this.sel1] = ipwwgt[this.sel1] :/ mean(ipwwgt[this.sel1], this.W[this.sel1])

  return(ipwwgt)
}

// Define function to calculate coefficients for a logit regression model
// A contant term is added to the model and its coefficient is included in the vector of betas
real rowvector gmatch::logitfit(real colvector Y, real matrix X , | real colvector W)
{
  transmorphic M
  M=moptimize_init()
  moptimize_init_evaluator(M,&logit_eval())
  moptimize_init_evaluatortype(M,"lf")
  moptimize_init_eq_cons(M, 1, "on")
  moptimize_init_depvar(M,1,Y)
  moptimize_init_eq_indepvars(M,1,X)
  moptimize_init_eq_colnames(M,1,(J(1,cols(X),"x") + strofreal((1..cols(X)))))
  moptimize_init_vcetype(M, "robust")
  if (args()>2 & W!=1) moptimize_init_weight(M, W)

  moptimize(M)
  // /* */ "Logit model coefficients and robust standard errors:"; moptimize_result_display(M)
  return(moptimize_result_coefs(M))
}

void logit_eval(transmorphic S, real rowvector beta, real colvector lnf)
{
  real colvector  Y, pm, xb, lj
  Y  = moptimize_util_depvar(S, 1)
  xb = moptimize_util_xb(S, beta, 1)
  pm = 2*(Y :!= 0) :- 1
  lj = invlogit(pm:*xb)

  if (any(lj :== 0)) {
    lnf = .
    return
  }
  lnf  = ln(lj)
}


end

include C:\Users\kkranker\Documents\Stata\Ado\Devel\gmatch\gmatch.ado
