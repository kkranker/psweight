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
    void             calcmeans(), calcvariances(), calcN(), calccovariances()
    real scalar      N1, N0, N, N1_raw, N0_raw, N_raw
    real scalar      mean_sd_sq(),  entropydistance()
    real rowvector   olsbeta(), diagvariance(), logitfit(), sd_sq(), asd()
    real colvector   olspredict(), logitpredict(), logitweights(), cbps_moments(), trim()
    real matrix      cbps_wgt_matrix(), cbpsgradient()
    real rowvector   means1, means0, variances0, variances1, variancesP, variancesA
    real matrix      covariances0, covariances1, covariancesP, covariancesA

  public:
    void             new(), set(), set_W(), set_Y(), clone(), multweight(), cbpseval()
    real rowvector   diff(), stddiff(), varratio(), prognosticdiff(), pomean()
    real scalar      mean_asd(), max_asd()
    real colvector   ipw(), cbps()
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
  /* */  "X is " + strofreal(rows(this.X)) + " by " + strofreal(cols(this.X))
  /* */  "X contains" ; this.varlist

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
  strofreal(this.N0_raw) + " control obs (sum of weights = " + strofreal(this.N0) + ")"
  strofreal(this.N1_raw) + " treatment obs (sum of weights = " + strofreal(this.N1) + ")"
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
  if (min((this.N0,this.N1)==0)) _error("Sum of weights is 0 in the treatment or control group.")

  // these means/varinaces are saved internally in the class (to avoid computing them over and over).
  // They need to be reset because we just reweighted the sample.
  // If I'm re-calcuating sample sizes, this is probably the case.  Set to missing here just to be safe.
  this.means0 = this.means1 = this.variances0 = this.variances1 = this.variancesP = this.variancesA = J(1,0,.)
  this.covariances0 = this.covariances1 = this.covariancesP = this.covariancesA = J(0,0,.)

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
void gmatch::multweight(|real colvector newweight)
{
  if (args()<1) newweight=1
  this.W = this.W_orig :* newweight
  // some new weights could be zero.  recalculate N and set means/variances to missing.
  this.calcN()
  this.means0 = this.means1 = this.variances0 = this.variances1 = this.variancesP = this.variancesA = J(1,0,.)
  this.covariances0 = this.covariances1 = this.covariancesP = this.covariancesA = J(0,0,.)
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
  (( "mean (T)",
     "mean (C)",
     "diff()",
     "stddiff()",
     (denominator==0 ? "sd (C)" : (denominator==1 ? "sd (T)" : (denominator==2 ? "sd (Pooled)" : (denominator==3 ? "sd (Avg)" : "")))),
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


real scalar gmatch::entropydistance(real colvector x, | real colvector w) {
  real colvector e
  if (args()<2) w=1
  e = x :* ln( x :* length(x) )
  return( mean(e, 1) )
  //return( quadcolsum( e ) )
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

// This function calculates the variances for the T and C group,
// and saves the results in private variables
void gmatch::calccovariances()
{
  if (all(this.W:==1)) {
    this.covariances0 = quadvariance(this.X[this.sel0, .])
    this.covariances1 = quadvariance(this.X[this.sel1, .])
    this.covariancesP = quadvariance(this.X)
  }
  else {
    this.covariances0 = quadvariance(this.X[this.sel0, .], this.W[this.sel0])
    this.covariances1 = quadvariance(this.X[this.sel1, .], this.W[this.sel1])
    this.covariancesP = quadvariance(this.X, this.W)
  }
  this.covariancesA = (this.covariances0 :+ this.covariances1) :/ 2
  // /* */ "Control group covariances:"; this.covariances0
  // /* */ "Treatment group covariances:"; this.covariances1
  // /* */ "Pooled covariances:"; this.covariancesP
  // /* */ "Average of covariances from treatment and control groups"; this.covariancesA
  this.variances0 = diagonal(this.covariances0)'
  this.variances1 = diagonal(this.covariances1)'
  this.variancesP = diagonal(this.covariancesP)'
  this.variancesA = diagonal(this.covariancesA)'
  // /* */ "Average of variances from treatment and control groups"; this.variancesA
}

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

// functions to return mean/max absolute standardized differences
real rowvector gmatch::asd(| real scalar denominator)
{
  if (args()<1) denominator=1
  return(abs(this.stddiff(denominator)))
}
real rowvector gmatch::sd_sq(| real scalar denominator)
{
  if (args()<1) denominator=1
  return(this.stddiff(denominator):^2)
}
real scalar gmatch::mean_asd(| real scalar denominator)
{
  if (args()<1) denominator=1
  return(mean(this.asd()'))
}
real scalar gmatch::max_asd(| real scalar denominator)
{
  if (args()<1) denominator=1
  return(max(this.asd()))
}
real scalar gmatch::mean_sd_sq(| real scalar denominator)
{
  if (args()<1) denominator=1
  return(mean(this.stddiff(denominator)')^2)
}


//  I don't know what this is anymore
// real scalar gmatch::cbps_moments(| real scalar denominator) {
//   if (args()<1) denominator=1
// // raw or weighted N?
//  return(sqrt(quadcolsum(this.X:*this.W)*invsym((this.X:*this.T)'this.X)*quadcolsum(this.X:*this.W)':/this.N_raw^2:*this.N1_raw))
// 	return(sqrt(quadcolsum(this.X:*this.W)*invsym((this.X:*this.T)'this.X)*quadcolsum(this.X:*this.W)':/this.N^2:*this.N1))
// }


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
    yhat[.,c] = this.olspredict(this.X, beta)
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
// If cols(X)+1==cols(beta), the function assumes the last coefficient corresponds to the constant term, and X just doesn't have a constant term
// Warning: this function doesn't check the conformability; I rely on Stata to produce errors with invalid arguments
real colvector gmatch::olspredict(real matrix X, real rowvector beta)
{
  if (cols(X)+1==cols(beta)) {
    return((X*beta[1..(cols(beta)-1)]') :+ beta[cols(beta)])
  }
  else {
    return(X*beta')
  }
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
  /* */ "logit model beta:"; beta
  pscore = this.logitpredict(this.X, beta)
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
// If cols(X)+1==cols(beta), the function assumes the last coefficient corresponds to a constant term, and X just doesn't include it
// Warning: this function doesn't check the conformability; I assume Stata will produce an error with invalid arguments
real colvector gmatch::logitpredict(real matrix X, real rowvector beta)
{
 if (cols(X)+1==cols(beta)) {
    return(invlogit((X*beta[1..(cols(beta)-1)]') :+ beta[cols(beta)]))
  }
  else if (cols(X)==cols(beta)) {
    return(invlogit(X*beta'))
  }
  else _error("X and beta have non-conformable dimensions")
}

// trims X at val and (1-val)
real colvector gmatch::trim(real colvector x, | real scalar val)
{
  real colvector xt
  if (args()<2) val = 1e-6
  xt = rowmin((J(rows(x),1,1-val),x))
  xt = rowmax((J(rows(x),1,  val),xt))
  return(xt)
}


// This turns a vector of pscores into IPW weights. this assumes a logit setup.
// Formulas match the normalized weights in Stata's teffects IPW command
//    pscore is a vector of propensity scores
//    est = "ate"  computes weights for average treatment effect (the default)
//        = "atet" computes weights for average treatment effect on the treated
//        = "ateu" computes weights for average treatment effect on the untreated
real colvector gmatch::logitweights(real colvector pscore, | string scalar est)
{
  real colvector pm,  minmax
  real matrix ipwwgt
  if (args()<2) est="ate"

  if (any(pscore:<=0) | any(pscore:>=1)) _error("Propensity scores need to be greater than 0 and less than 1.")
//  /* */ if (minmax[1,1]<=0.03 & (strlower(est)=="ate" | strlower(est)=="ateu")) errprintf("Warning: minimum propensity score is %12.0g \n", minmax[1,1])
//  /* */ if (minmax[1,2]>=0.97 & (strlower(est)=="ate" | strlower(est)=="atet")) errprintf("Warning: maximum propensity score is %12.0g \n", minmax[1,2])

  pm = 1 :- (!this.T)
  if      (strlower(est)=="ate")   ipwwgt = (pm :* (1:/pscore)) :+ (!pm :* (1:/(1:-pscore)))
  else if (strlower(est)=="atet")  ipwwgt =  pm :+ (!pm :* (pscore:/(1:-pscore)))
  else if (strlower(est)=="ateu")  ipwwgt = !pm :+ ( pm :* ((1:-pscore):/pscore))
  else _error(est + " is an invalid argument for gmatch::logitweights()")

  // normalize the weights to have mean 1 in each group
  if (strlower(est)=="ate" | strlower(est)=="atet") ipwwgt[this.sel0] = ipwwgt[this.sel0] :/ mean(ipwwgt[this.sel0], this.W[this.sel0])
  if (strlower(est)=="ate" | strlower(est)=="ateu") ipwwgt[this.sel1] = ipwwgt[this.sel1] :/ mean(ipwwgt[this.sel1], this.W[this.sel1])
  return(ipwwgt)
}

// Define function to calculate coefficients for a logit regression model
// A contant term is added to the model and its coefficient is included in the vector of betas
real rowvector gmatch::logitfit(real colvector Y, real matrix X, | real colvector W)
{
  transmorphic S
  S=moptimize_init()
  moptimize_init_evaluator(S, &logit_eval())
  moptimize_init_evaluatortype(S,"lf")
  moptimize_init_eq_cons(S, 1, "on")
  moptimize_init_depvar(S,1,Y)
  moptimize_init_eq_indepvars(S,1,X)
  moptimize_init_eq_colnames(S,1,(J(1,cols(X),"x") + strofreal((1..cols(X)))))
  moptimize_init_vcetype(S, "robust")
  if (args()>2 & W!=1) moptimize_init_weight(S, W)

  moptimize(S)
  // /* */ "Logit model coefficients and robust standard errors:"; moptimize_result_display(S)
  return(moptimize_result_coefs(S))
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

// function that returns CBPS weights
//    est corresponds to the options in gmatch::logitweights()
//        "ate"  computes weights for average treatment effect (the default)
//        "atet" computes weights for average treatment effect on the treated
//        "ateu" computes weights for average treatment effect on the untreated
//    fctn corresponds to the balance measure
//        "mean_sd_sq" minimizes the mean standardized difference squared
//    denominator is passed to stddiff() and related functions
real colvector gmatch::cbps(| string scalar est, string scalar fctn, real scalar denominator)
{
  real rowvector beta
  real colvector pscore, cbpswgt
  class gmatch scalar M
  if (args()<1) est="ate"
  if (args()<2) fctn="sd_sq"
  if (args()<3) denominator=1

  M.clone(this)

  transmorphic S
  S=optimize_init()
  optimize_init_evaluator(S, &cbps_eval())
  optimize_init_which(S, "min")
  optimize_init_argument(S, 1, M)
  optimize_init_argument(S, 2, est)
  optimize_init_argument(S, 3, fctn)
  optimize_init_argument(S, 4, denominator)
	optimize_init_singularHmethod(S,"hybrid")  // equivalent to ml's "difficult" option
  optimize_init_technique(S, "bfgs 15 nr 15")
	optimize_init_conv_ptol(S, 1e-7)
	optimize_init_conv_vtol(S, 1e-8)
	optimize_init_conv_nrtol(S, 1e-6)
  optimize_init_tracelevel(S, "value" )  // "none", "value", "params"


  if (fctn=="cbpsloss" | fctn=="sd_sq" | fctn=="asd" | fctn=="mean_sd_sq_ent" | fctn=="mean_sd_sq_cv") {
    optimize_init_evaluatortype(S,"gf0")
	  optimize_init_conv_ignorenrtol(S, "on")
  }
  else if (fctn=="cbps_moments") {
    // optimize_init_evaluatortype(S,"d1")
    optimize_init_conv_ptol(S, 1e-11)
  	optimize_init_conv_vtol(S, 1e-12)
    optimize_init_tracelevel(S, "none" )  // "none", "value", "params"
// this should be automatic
    optimize_init_evaluatortype(S,"d1")   // d1 if I'm running plain vanilla. otherwise just use "do"
  optimize_init_evaluatortype(S,"d0")    // numerical gradient
  optimize_init_evaluatortype(S,"gf1")  // for overidentified version
	  optimize_init_conv_ignorenrtol(S, "on")
  }
  else {
    optimize_init_evaluatortype(S,"d0")
	  optimize_init_conv_ignorenrtol(S, "off")
  }


  "Step 1 (initial values from logit model):"
  real rowvector betalogit
  betalogit = M.logitfit(M.T, M.X, M.W)
  optimize_init_params(S, betalogit)
  // /* */ "  optimize_init_params(S)";   optimize_init_params(S)
  // /* */ "optimize_result_value0(S)"; optimize_result_value0(S)

  
  if (fctn=="cbps_moments") {
    external real matrix ww
    ww = M.cbps_wgt_matrix(betalogit, 0, est)
    ww = invsym(ww)
  }
  else ww = .
  optimize_init_argument(S, 5, denominator)
  ""

  "Step 2 (CBPS) :"
  /* */ // This temp code keeps optimize() from producing an error
  /* */ // Once it's working switch back to
  /* */ // (void) optimize(S)
  /* */
  /* */
  /* */ (void) _optimize(S)
  /* */               if (optimize_result_returncode(S)!=0) {
  /* */                       errprintf("{p}\n")
  /* */                       errprintf("%s\n", optimize_result_errortext(S))
  /* */                       "current beta"; optimize_result_params(S)
  /* */                       errprintf("\n \n Exiting the function early.\n")
  /* */                       errprintf("{p_end}\n")
  /* */                       // exit(optimize_result_returncode(S))
  /* */                       return(J(this.N,1,.))
  /* */               }


  beta    = optimize_result_params(S)
  /* */ "CBPS beta"; ((M.varlist,"_cons")', strofreal(beta)')
  // /* */ "optimize_result_value(S)" ; optimize_result_value(S)
  // /* */ "optimize_result_scores(S)" ; optimize_result_scores(S)
  // /* */ optimize_query(S)

  pscore  = M.logitpredict(M.X, beta)
  pscore  = M.trim(pscore)
  cbpswgt = M.logitweights(pscore, est)
  M.multweight(cbpswgt)

  "Balance after CBPS (" + fctn + "):"
  /* */   if      (fctn=="mean_sd_sq")      M.mean_sd_sq(denominator)
  /* */   else if (fctn=="mean_asd")        M.mean_asd(denominator)
  /* */   else if (fctn=="max_asd")         M.max_asd(denominator)
//  /* */   else if (fctn=="cbps_moments")      M.cbps_moments(denominator)
  /* */   else if (fctn=="sd_sq")           M.sd_sq(denominator)'
  /* */   else if (fctn=="asd")             M.asd(denominator)'
//  /* */   else if (fctn=="cbpsloss2")       M.cbpsloss2()'

  /* */ "optimize_result_value(S)" ; optimize_result_value(S)
  /* */ "balance table after matching (" + fctn + "):"; real matrix temp; temp = M.balancetable(denominator)
  /* */ "entropydistance of control weights (" + fctn + "):"; (this.entropydistance(cbpswgt[this.sel0], this.W_orig[this.sel0]))
  /* */ "cv of control weights (" + fctn + "):"; (sqrt(this.diagvariance(cbpswgt[this.sel0], this.W_orig[this.sel0])) :/  mean(cbpswgt[this.sel0], this.W_orig[this.sel0]))


  return(cbpswgt)
}



void gmatch::cbpseval( real    /* scalar      */   todo,
                       real    /* rowvector   */   beta,
                       string  /* scalar      */   est,
                       string  /* scalar      */   fctn,
                       real    /* scalar      */   denominator,
                       real    /* scalar      */   ww,
                       real    /* colvector   */   lnf,
                       real    /* rowvector   */   g,
                       real    /* matrix      */   H)
{
   real colvector  pscore, cbpswgt
   pscore = this.logitpredict(this.X, beta)
   pscore = this.trim(pscore)
   cbpswgt = this.logitweights(pscore, est)

    /* if (fctn=="cbpsloss") {
                                     lnf = this.cbpsloss(pscore,cbpswgt)'
     return
   }
   else if (fctn=="cbpslossOID") {
                                     lnf = this.cbpslossOID(pscore,cbpswgt)'
     return
   }
   else */  if (fctn=="cbps_moments") {
     real matrix dpscore, gg, grad, oid
     oid = 1
     dpscore = pscore:*(1:-pscore)
     gg = this.cbps_moments(pscore, dpscore, oid, est)
     lnf = gg' * ww * gg
 ///* */  "lnf is " + strofreal(rows(lnf)) + " by " + strofreal(cols(lnf))
     if (todo==0) return
     grad = this.cbpsgradient(pscore, oid , est)
 ///* */  "grad is " + strofreal(rows(grad)) + " by " + strofreal(cols(grad))
     g = grad' * ww * gg :* (2:*this.N)
     g = g'

 ///* */  "g is " + strofreal(rows(g)) + " by " + strofreal(cols(g))
     return
   }
   this.multweight(cbpswgt)
   if      (fctn=="mean_sd_sq")      lnf = this.mean_sd_sq(denominator)
   else if (fctn=="mean_asd")        lnf = this.mean_asd(denominator)
   else if (fctn=="max_asd")         lnf = this.max_asd(denominator)
//   else if (fctn=="cbps_moments")      lnf = this.cbps_moments(denominator)
   else if (fctn=="sd_sq")           lnf = this.sd_sq(denominator)'
   else if (fctn=="mean_sd_sq_ent")  lnf = this.mean_sd_sq(denominator) + .25*(this.entropydistance(cbpswgt[this.sel0], this.W_orig[this.sel0]))
   else if (fctn=="mean_sd_sq_cv")   lnf = this.mean_sd_sq(denominator) + .25*(sqrt(this.diagvariance(cbpswgt[this.sel0], this.W_orig[this.sel0])) :/  mean(cbpswgt[this.sel0], this.W_orig[this.sel0]))
   else if (fctn=="asd")             lnf = this.asd(denominator)'
//   else if (fctn=="cbpsloss2")       lnf = this.cbpsloss2()'
   else                              _error(fctn + " is invalid with gmatch::cbpseval()")

}


void cbps_eval(real todo, real beta,
               class gmatch scalar M,
               string est, string fctn, real denominator, real ww,
               real lnf, real g, real H)
{
  M.cbpseval(todo,beta,est,fctn,denominator,ww,lnf,g,H)
}


// Port of the moment function from the Stata verion of CBPS
real colvector gmatch::cbps_moments(real colvector pscore, real matrix dpscore, real scalar overid, string scalar est)
{
  real matrix X_const
  real colvector gg
  if (any(this.W:!=1)) _error("gmatch::cbpslossS() does not yet accomodate weighted samples")
  
	X_const = (this.X, J(this.N,1,1))

  if (strlower(est)=="ate") {
			gg=X_const'*((this.T-pscore):/pscore:/(1:-pscore))
	}
  else if (strlower(est)=="atet") {
			gg=X_const'*((this.T-pscore):/(1:-pscore):*(this.N_raw/this.N1_raw))
	}
  else _error(est + " is invalid with gmatch::cbps_moments()")

  if(overid) {
    gg = (X_const'*(dpscore:*(this.T-pscore):/pscore:/(1:-pscore)) \ gg)
	}

  gg = gg:/this.N_raw

  return(gg)
}


// Port of the gradient function from the Stata verion of CBPS
real matrix gmatch::cbpsgradient(real colvector pscore, real scalar overid, string scalar est)
{
  real matrix G, dw,  X_const

	X_const = (this.X, J(this.N,1,1))

  if (strlower(est)=="ate") {
    G = -(X_const:*((this.T:-pscore):^2):/pscore:/(1:-pscore))'X_const
  }
  else if (strlower(est)=="atet") {
    dw=(pscore:*(this.T:-1)):/(1:-pscore):*(this.N/this.N1)
    G = ((X_const:*dw)') * X_const
  }
  if (overid) {
    G = -(X_const:*(pscore:*(1:-pscore)))' X_const \ G
  }
  G = G :/ this.N
  return(G)
}


// Port of the weighting matrix function from the Stata verion of CBPS
real matrix gmatch::cbps_wgt_matrix(real rowvector beta, real scalar overid, string scalar est)
{
  real matrix X_const, ww
  real colvector pscore, dpscore

  pscore  = this.logitpredict(this.X, beta)
  pscore  = this.trim(pscore)
  dpscore = pscore:*(1:-pscore)

  X_const = (this.X, J(this.N_raw,1,1))

	if (!overid) {
		if (strlower(est)=="ate") {
      /* */ // Can I use quadcross with dpscore as a weight
			ww = (X_const:/(pscore:*(1:-pscore)))'*X_const  
		}
		else if (strlower(est)=="atet") {
			ww = (X_const:*(pscore:/(1:-pscore)):*(this.N_raw/this.N1_raw):^2)'*X_const
		}
	}
  else {
		if (strlower(est)=="ate") {
			ww = (      (X_const:*(dpscore:^2:/pscore:/(1:-pscore)))'*X_const, 
                  (X_const:*(dpscore:/pscore:/(1:-pscore)))'*X_const)
			ww = ( ww \ (X_const:*(dpscore:/pscore:/(1:-pscore)))'*X_const, 
                  (X_const:*(1:/pscore:/(1:-pscore)))'*X_const)
	  }
	  else if (strlower(est)=="atet") {
			ww = (      (X_const:*(pscore:/(1:-pscore):*dpscore:^2:/pscore:^2))'*X_const, 
                  (X_const:*(pscore:/(1:-pscore):*dpscore:/pscore):*(this.N_raw/this.N1_raw))'*X_const)
			ww = ( ww \ (X_const:*(pscore:/(1:-pscore):*dpscore:/pscore):*(this.N_raw/this.N1_raw))'*X_const,  
                  (X_const:*(pscore:/(1:-pscore)):*((this.N_raw/this.N1_raw)^2))'*X_const )
		}
	}
	ww=ww:/this.N_raw
  return(ww)
}


end
// include C:\Users\kkranker\Documents\Stata\Ado\Devel\gmatch\gmatch-port-R-functions.mata
include C:\Users\kkranker\Documents\Stata\Ado\Devel\gmatch\gmatch.ado

