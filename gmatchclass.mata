//****************************************************************************/
*! $Id: gmatch.ado,v dbef04cb535b 2018/02/14 06:08:14 kkranker $
*! Generalization of IPW and CBPS estimators
*! Class defintion: gmatch()
//
*! By Keith Kranker
// Last updated $Date: 2018/02/14 06:08:14 $
//
// Copyright (C) Mathematica Policy Research, Inc.
// This code cannot be copied, distributed or used without the express written
// permission of Mathematica Policy Research, Inc.
//*****************************************************************************/

version 15.1
mata:

mata set matastrict on
mata set matafavor speed

class gmatch
{
  private:
    real colvector   T, W, sel1, sel0, Y0, W_orig, W_mtch, PS_mtch
    real matrix      X, XC, Xstd
    string scalar    treatvar, depvars, wgtvar
    string rowvector varlist
    real rowvector   means1, means0, meansP, variances0, variances1, variancesP, variancesA
    real matrix      covariances0, covariances1, covariancesP, covariancesA
    void             calcN(), cbps_port_stata(), cbps_port_r(), postbeta()
    real scalar      K, N1, N0, N, N1_raw, N0_raw, N_raw
    real scalar      mean_sd_sq(), entropydistance()
    real rowvector   olsbeta(), diagvariance(), logitbeta(), sd_sq(), asd(), wgt_moments()
    real colvector   olspredict(), logitpredict(), logitweights(), cbps_port_stata_moments(), trim()
    real matrix      cbps_port_stata_wgt_matrix(), cbps_port_stata_gradient(), Ct()

  public:
    void             new(), set(), set_Y(), reweight(), get_scores()
    void             cbpseval(), balanceresults()
    real rowvector   gmatch(), ipw(), cbps(), cbpsoid()
    real rowvector   diff(), stddiff(), varratio(), progdiff(), stdprogdiff(), pomean()
    real rowvector   means1(), means0(), meansP(), variances0(), variances1(), variancesP(), variancesA()
    real matrix      covariances0(), covariances1(), covariancesP(), covariancesA()
    real scalar      mean_asd(), max_asd(), wgt_cv(), wgt_sd(), wgt_skewness(), wgt_kurtosis(), wgt_max()
    real matrix      balancetable()
}


// The following function is called every time a class instance is created
void gmatch::new()
{
  this.depvars = ""
}

// loads the main data into the class, using views wherever possible
void gmatch::set(string scalar treatvar, string scalar varlist, string scalar tousevar, | string scalar wgtvar)
{
  // Define treatment dummy
  this.treatvar = treatvar
  st_view(this.T, ., treatvar, tousevar)
  // /* */  "T is " + strofreal(rows(this.T)) + " by " + strofreal(cols(this.T))

  // Define covariates
  this.varlist  = tokens(varlist)
  st_view(this.X , .,    this.varlist                , tousevar)
  this.K = cols(this.X)
  // /* */  "X contains" ; this.varlist
  // /* */  "X is " + strofreal(rows(this.X)) + " by " + strofreal(cols(this.X))

  // Define weights
  // This code assumes weights are **already** normalized. Here's code to normalize: this.W = this.W :/ (rows(this.W) / quadcolsum(this.W))
  if (args()>=4) {
    this.wgtvar = wgtvar
    st_view(this.W_orig, ., this.wgtvar, tousevar) // an extra copy of the weight variable that can only be set via this function. Useful for reweighting/matching situations.
  }
  else this.W_orig = J(rows(this.T),1,1)
  // /* */  "W is " + strofreal(rows(this.W)) + " by " + strofreal(cols(this.W))
  // initialize W_mtch=1 and W=W_orig
  // W_orig is a view, but W and W_mtch are not
  this.reweight()
  // /* */  "W_orig is " + strofreal(rows(this.W_orig )) + " by " + strofreal(cols(this.W_orig))
  // /* */  "W is " + strofreal(rows(this.W)) + " by " + strofreal(cols(this.W))

  // Index to select observations in control and treatment groups
  this.sel0 = selectindex(!this.T :& this.W_orig)
  this.sel1 = selectindex( this.T :& this.W_orig)

  // Save raw number of observations
  this.N0_raw = rows(this.sel0)
  this.N1_raw = rows(this.sel1)
  this.N_raw = this.N0_raw + this.N1_raw
  if (min((this.N0_raw,this.N1_raw)==0)) _error("At least one treatment and control observation required.")

  this.calcN()
  if (all(this.W_orig:==1)) {
    strofreal(this.N0_raw) + " control obs"
    strofreal(this.N1_raw) + " treatment obs"
    "(Data are unweighted.)"
  }
  else {
    strofreal(this.N0_raw) + " control obs (sum of weights = " + strofreal(this.N0) + ")"
    strofreal(this.N1_raw) + " treatment obs (sum of weights = " + strofreal(this.N1) + ")"
  }
}

// flags observations with weights!=0
// calculates wegihted sample sizes for treatment and control groupss
void gmatch::calcN() {
  this.N0 = quadcolsum(this.W[this.sel0])
  this.N1 = quadcolsum(this.W[this.sel1])
  this.N = this.N0 + this.N1
  if (min((this.N0,this.N1)==0)) _error("Sum of weights is 0 in the treatment or control group.")
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

// multiply the original weights by a matching weight and store in this.W
//   without the first argument, the function will reset the weights back to the original.
// The second (optional) argument will store IPW weights in this.PS_mtch
// Afterward, N's are re-calculated and means/(co)variances are re-set to missing.
void gmatch::reweight(|real colvector newweight, real colvector newpscores)
{
  // weights
  if (args()<1) this.W_mtch = J(rows(this.T),1,1)
  else          this.W_mtch = newweight
  this.W = this.W_orig :* this.W_mtch
  // p-scores
  if (args()>1) this.PS_mtch = newpscores
  else if (args()<1) this.PS_mtch = .
  // recalculate N and set means/variances to missing.
  this.calcN()
  this.means0 = this.means1 = this.meansP = this.variances0 = this.variances1 = this.variancesP = this.variancesA = J(1,0,.)
  this.covariances0 = this.covariances1 = this.covariancesP = this.covariancesA = J(0,0,.)
}

// used to push the resulting weights and propensity scores back into Stata.
// newvarnames needs to have four new variable names
//   - the first  variable is filled with this.W
//   - the second variable is filled with this.W_mtch
//   - the third  variable is filled with this.PS_mtch
//   - the fourth variable is filled with this.T
// The touse variable is required.
// This function was written for the .ado program.
// Use this with caution in other contexts; it asssumes the data have not been edited
// or re-sorted since the data were originally loaded into the class instance using gmatch::set().
void gmatch::get_scores(string rowvector newvarnames, string scalar tousevar)
{
  real matrix thisview
  newvarnames = tokens(newvarnames)
  if (length(newvarnames)!=4) _error("gmatch::get_scores() requires four variable names")
  st_view(thisview, ., newvarnames, tousevar)

  if (rows(thisview)==rows(this.W)) thisview[.,1] = this.W
  else thisview[.,1] = J(rows(thisview),1,.)

  if (rows(thisview)==rows(this.W_mtch)) thisview[.,2] = this.W_mtch
  else thisview[.,2] = J(rows(thisview),1,.)

  if (rows(thisview)==rows(this.PS_mtch)) thisview[.,3] = this.PS_mtch
  else thisview[.,3] = J(rows(thisview),1,.)

  if (rows(thisview)==rows(this.T)) thisview[.,4] = this.T
  else thisview[.,4] = J(rows(thisview),1,.)
}

// This function makes a balance table and prints it to the screen
// The argument is the same as the definition in stddiff() and varratio()
// results are also saved in stata in r(bal)
real matrix gmatch::balancetable(| real scalar denominator)
{
  real matrix table
  string rowvector colstripe, tmp
  if (args()<1) denominator=1

  table = ( this.means1()
          \ this.means0()
          \ this.diff()
          \ this.stddiff(denominator)
          \ (denominator==0 ? sqrt(this.variances0()) : (denominator==1 ? sqrt(this.variances1()) : (denominator==2 ? sqrt(this.variancesP()) : (denominator==3 ? sqrt(this.variancesA()) : _error("denominator argument invalid")))))
          \ this.varratio())'

  colstripe = ("mean_T",
               "mean_C",
               "diff",
               "std_diff",
               (denominator==0 ? "sd_C" : (denominator==1 ? "sd_T" : (denominator==2 ? "sd_pool" : (denominator==3 ? "sd_avg" : "")))),
               "var_ratio")

  st_matrix("r(bal)", table)
  st_matrixcolstripe("r(bal)", (J(length(colstripe),1,""), colstripe'))
  st_matrixrowstripe("r(bal)", (J(length(varlist)  ,1,""), varlist'  ))

  tmp=st_tempname()
  stata("matrix "+tmp+"=r(bal)")
  stata("_matrix_table "+tmp+","+st_local("diopts")); ""

  return(table)
}

// This function prints the balance table to the screen
// The argument is the same as their definition in stddiff() and varratio()
void gmatch::balanceresults(| string scalar est, real scalar denominator)
{
  if (args()<1) est="ate"
  if (args()<2) denominator=1
  transmorphic temp
  if (all(this.W_mtch:==1)) "Unmatched data"
  st_rclear()
  "Balance:";                            temp = this.balancetable(denominator)
  "Mean standardized diff., squared";    this.mean_sd_sq(denominator)
  "Mean absolute standardized diff.";    this.mean_asd(denominator)
  "Maximum absolute standardized diff."; this.max_asd(denominator)
  "C.V. of matching weights:";           this.wgt_cv(est)
  "S.D. of matching weights:";           this.wgt_sd(est)
  "Skewness of matching weights:";       this.wgt_skewness(est)
  "Kurtosis of matching weights:";       this.wgt_kurtosis(est)
  "Maximum matching weight:";            this.wgt_max(est)
  if (this.depvars!="") {
    "Prognostic scores:";                temp = this.progdiff()
  }
}

// These functions calculate the (weighted) means for the C group, the T group, and the pooled sample
// These means are saved internally in the class (to avoid computing them over and over).
// If means have previously been calculated, they are simply returned.
real rowvector gmatch::means0()
{
  if (!length(this.means0)) this.means0 = mean(this.X[this.sel0, .], this.W[this.sel0])
  return(this.means0)
}
real rowvector gmatch::means1()
{
  if (!length(this.means1)) this.means1 = mean(this.X[this.sel1, .], this.W[this.sel1])
  return(this.means1)
}
real rowvector gmatch::meansP()
{
  if (!length(this.meansP)) this.meansP = mean(this.X, this.W)
  return(this.meansP)
}

// This function calculates the difference in means between the T and C groups
real rowvector gmatch::diff()
{
  return(this.means1() :- this.means0())
}

// This function calculates entropy for a vector x
real scalar gmatch::entropydistance(real colvector x, | real colvector w) {
  real colvector e
  real scalar sumw
  if (args()<2) {
    w=1
    sumw=rows(x)
  }
  else {
    sumw=quadcolsum(w)
  }

  // for me, the sum of weights is rows(x) (unweighted) or sum(w) (weighted)
  // in entropy balancing from Hainmueller et al., the sum of the weights = 1.
  e = x :* ln( x :* sumw )
  return( quadcolsum( e ) )
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

// These functions calculate the (weighted) variances for the C group, the T group, and the pooled sample, and the average of the T and C group variabnces
// These variances are saved internally in the class (to avoid computing them over and over).
// If variances have previously been calculated, they are simply returned.
real rowvector gmatch::variances0()
{
  if (!length(this.variances0)) this.variances0 = this.diagvariance(this.X[this.sel0, .], this.W[this.sel0], this.means0())
  return(this.variances0)
}
real rowvector gmatch::variances1()
{
  if (!length(this.variances1)) this.variances1 = this.diagvariance(this.X[this.sel1, .], this.W[this.sel1], this.means1())
  return(this.variances1)
}
real rowvector gmatch::variancesP()
{
  if (!length(this.variancesP)) this.variancesP = this.variancesP = this.diagvariance(this.X, this.W)
  return(this.variancesP)
}
real rowvector gmatch::variancesA()
{
  if (!length(this.variancesA)) this.variancesA = (this.variances0() :+ this.variances1()) :/ 2
  return(this.variancesA)
}

// These functions calculate the (weighted) covariance matrices for the C group, the T group, and the pooled sample, and the average of the T and C group variabnces
// These covariance matrices are saved internally in the class (to avoid computing them over and over).
// If covariance matrices have previously been calculated, they are simply returned.
real matrix gmatch::covariances0()
{
  if (!length(this.covariances0)) {
    if (allof(this.W,1)) this.covariances0 = quadvariance(this.X[this.sel0, .])
    else                 this.covariances0 = quadvariance(this.X[this.sel0, .], this.W[this.sel0])
    this.variances0 = diagonal(this.covariances0)'
  }
  return(this.covariances0)
}
real matrix gmatch::covariances1()
{
  if (!length(this.covariances1)) {
    if (allof(this.W,1)) this.covariances1 = quadvariance(this.X[this.sel1, .])
    else                 this.covariances1 = quadvariance(this.X[this.sel1, .], this.W[this.sel1])
    this.variances1 = diagonal(this.covariances1)'
  }
  return(this.covariances1)
}
real matrix gmatch::covariancesP()
{
  if (!length(this.covariancesP)) {
    if (allof(this.W,1)) this.covariancesP = quadvariance(this.X)
    else                 this.covariancesP = quadvariance(this.X, this.W)
    this.variancesP = diagonal(this.covariancesP)'
  }
  return(this.covariancesP)
}
real matrix gmatch::covariancesA()
{
  if (!length(this.covariancesA)) {
    this.covariancesA = (this.covariances0() :+ this.covariances1()) :/ 2
    this.variancesA = diagonal(this.covariancesA)'
  }
  return(this.covariancesA)
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
  if      (denominator==0) stddiff = (this.diff() :/ sqrt(this.variances0()))
  else if (denominator==1) stddiff = (this.diff() :/ sqrt(this.variances1()))
  else if (denominator==2) stddiff = (this.diff() :/ sqrt(this.variancesP()))
  else if (denominator==3) stddiff = (this.diff() :/ sqrt(this.variancesA()))
  else _error(strofreal(denominator)+ " is an invalid argument for gmatch::stddiff()")
  return(stddiff)
}

// The following functions calculate UNWEIGHTED means, variance, CV, SD, skewness, kurtosis, and higher moments of the matching weights
// As usual, est controls whether we do with with all observations, the treatment group, or the control group
real rowvector gmatch::wgt_moments(real scalar r, string scalar est)
{
  real scalar v, m
  real colvector W_sel
  if      (strlower(est)=="ate" ) W_sel=this.W_mtch
  else if (strlower(est)=="atet") W_sel=this.W_mtch[this.sel0]
  else if (strlower(est)=="ateu") W_sel=this.W_mtch[this.sel1]
  else _error(est + " is an invalid argument for gmatch::wgt_moments()")
  m = mean(W_sel)
  if (r==0) { // the only exception is that r==0 gives the sd
    v = sqrt(quadcolsum((W_sel:-m):^2) / (rows(W_sel)-1))
  }
  else v = quadcolsum((W_sel:-m):^r)
  return((v,m))
}
real scalar gmatch::wgt_cv(string scalar est)
{
  real rowvector vm
  real scalar cv
  vm = this.wgt_moments(0,est)
  cv = vm[1]/vm[2]
  st_numscalar("r(wgt_cv)",cv)
  return(cv)
}
real scalar gmatch::wgt_sd(string scalar est)
{
  real scalar sd
  sd = this.wgt_moments(0,est)[1]
  st_numscalar("r(wgt_sd)",sd)
  return(sd)
}
real scalar gmatch::wgt_skewness(string scalar est)
{
  real scalar skew
  skew = (this.wgt_moments(3,est)[1]) * (this.wgt_moments(2,est)[1])^(-3/2)
  st_numscalar("r(wgt_skewness)",skew)
  return(skew)
}
real scalar gmatch::wgt_kurtosis(string scalar est)
{
  real scalar kurt
  kurt = (this.wgt_moments(4,est)[1]) * (this.wgt_moments(2,est)[1])^(-2)
  st_numscalar("r(wgt_kurtosis)",kurt)
  return(kurt)
}
real scalar gmatch::wgt_max(string scalar est)
{
  real scalar mx
  if      (strlower(est)=="ate" ) mx = max(this.W_mtch)
  else if (strlower(est)=="atet") mx = max(this.W_mtch[this.sel0])
  else if (strlower(est)=="ateu") mx = max(this.W_mtch[this.sel1])
  else _error(est + " is an invalid argument for gmatch::wgt_moments()")
  st_numscalar("r(wgt_max)",mx)
  return(mx)
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
  real scalar out
  if (args()<1) denominator=1
  out = mean(this.asd(denominator)')
  st_numscalar("r(mean_asd)",out)
  return(out)
}
real scalar gmatch::max_asd(| real scalar denominator)
{
  real scalar out
  if (args()<1) denominator=1
  out = max(this.asd(denominator))
  st_numscalar("r(max_asd)", out)
  return(out)
}
real scalar gmatch::mean_sd_sq(| real scalar denominator)
{
  real scalar out
  if (args()<1) denominator=1
  out = mean(this.stddiff(denominator)')^2
  st_numscalar("r(mean_sd_sq)", out)
  return(out)
}


// This function calculates ratio of variances between the T and C groups
real rowvector gmatch::varratio()
{
  return((this.variances1() :/ this.variances0()))
}


// function that returns difference in y_hat, where y_hat is generated using a
// OLS regression of y on X using the control group data
// Denominator is defined the same as in stddiff(), and is passed to stdprogdiff()
real rowvector gmatch::progdiff(| real scalar denominator)
{
  real rowvector beta, progdiff, yhat_bar_0, yhat_bar_1, y_bar_0, stdprogdiff
  real colvector yhat
  real scalar c
  real matrix table
  string rowvector colstripe,tmp
  if (!length(this.depvars)) _error("Dependent variable is undefined.  Use gmatch::set_Y().")
  if (args()<1) denominator=1

  yhat = J(rows(this.X), cols(this.Y0), .)
  for (c=1; c<=cols(this.Y0); c++) {
    beta = this.olsbeta(this.Y0[.,c], this.X[this.sel0,.], this.W[this.sel0])
    yhat[.,c] = this.olspredict(this.X, beta)
  }

  yhat_bar_0  = mean(yhat[this.sel0,.], this.W[this.sel0])
  yhat_bar_1  = mean(yhat[this.sel1,.], this.W[this.sel1])
  progdiff    = yhat_bar_1 :- yhat_bar_0
  stdprogdiff = stdprogdiff(denominator, yhat, progdiff)
  y_bar_0     = mean(this.Y0, this.W[this.sel0])

  table = (yhat_bar_1 \ yhat_bar_0 \ progdiff \ stdprogdiff \ y_bar_0)'

  colstripe = ("mean_yhat_T",
               "mean_yhat_C",
               "diff",
               "std_diff",
               "mean_y_C")

  st_matrix("r(prog)", table)
  st_matrixcolstripe("r(prog)", (J(length(colstripe),1,""), colstripe'))
  st_matrixrowstripe("r(prog)", (J(length(depvars)  ,1,""), depvars'  ))

  tmp=st_tempname()
  stata("matrix "+tmp+"=r(prog)")
  stata("_matrix_table "+tmp+","+st_local("diopts"))
  "Note: The std_diff column does not account for the standard error of the linear predictions."

  return(progdiff)
}

// This function calculates standardized differences in prognostic scores
// Variances do not account for the OLS modeling; they are just the variance of the y_hat variable
// Denominator is defined the same as in stddiff()
// Note: when this is used in the optimization program, the OLS model is re-estimated
// each iteration.  The reason is that it allows the prognostic scores to be estiamted
// with a reweighted comparison group that "looks like" the treatment group.
real rowvector gmatch::stdprogdiff(| real scalar denominator, real matrix yhat, real rowvector progdiff)
{
  real rowvector beta, yhat_bar_0, yhat_bar_1, stddiff
  real scalar c
  if (args()<1) denominator=1
  if (args()<2) {
    if (!length(this.depvars)) _error("Dependent variable is undefined.  Use gmatch::set_Y().")
    yhat = J(rows(this.X), cols(this.Y0), .)
    for (c=1; c<=cols(this.Y0); c++) {
      beta = this.olsbeta(this.Y0[.,c], this.X[this.sel0,.], this.W[this.sel0])
      yhat[.,c] = this.olspredict(this.X, beta)
    }
  }
  if (args()<3) {
    yhat_bar_0 = mean(yhat[this.sel0,.], this.W[this.sel0])
    yhat_bar_1 = mean(yhat[this.sel1,.], this.W[this.sel1])
    progdiff = yhat_bar_1 :- yhat_bar_0
  }
  if      (denominator==0) stddiff = (progdiff :/ sqrt(this.diagvariance(yhat[this.sel0, .], this.W[this.sel0])))
  else if (denominator==1) stddiff = (progdiff :/ sqrt(this.diagvariance(yhat[this.sel1, .], this.W[this.sel1])))
  else if (denominator==2) stddiff = (progdiff :/ sqrt(this.diagvariance(yhat, this.W)))
  else if (denominator==3) stddiff = (progdiff :/ sqrt((this.diagvariance(yhat[this.sel0, .], this.W[this.sel0]) :+ this.diagvariance(yhat[this.sel1, .], this.W[this.sel1])) :/ 2))
  else _error(strofreal(denominator)+ " is an invalid argument for gmatch::stddiff()")
  return(stddiff)
}

// Define function to calculate coefficients for an OLS regression model
// A contant term is included in the regression.
real rowvector gmatch::olsbeta(real matrix y, real matrix X, | real colvector w, real scalar addconst)
{
  real colvector beta
  real matrix XX, Xy
  if (args()<3) w=1
  if (args()<4) addconst=1

  if (addconst) {
    XX = quadcross(X, 1, w, X, 1)
    Xy = quadcross(X, 1, w, y, 0)
    beta = invsym(XX,(cols(X)+1))*Xy
  }
  else {
    XX = quadcross(X, 0, w, X, 0)
    Xy = quadcross(X, 0, w, y, 0)
    beta = invsym(XX)*Xy
  }
  return(beta')
}

// Function that returns predicted values, X*beta'
// If cols(X)+1==cols(beta), the function assumes the last coefficient corresponds to the constant term, and X just doesn't have a constant term
// Warning: this function doesn't check the conformability; I rely on Stata to produce errors with invalid arguments
real colvector gmatch::olspredict(real matrix X, real rowvector beta)
{
  if ((cols(X)==cols(beta)) & cols(beta)) {
    return(X*beta')
  }
  else if ((cols(X)==cols(beta)-1) & cols(beta)) {
    return((X*beta[1..(cols(beta)-1)]') :+ beta[cols(beta)])
  }
  else _error("X and beta are not conformable.")
}


// Function that performs IPW
//    est corresponds to the options in gmatch::logitweights()
//    est = "ate"  computes weights for average treatment effect (the default)
//        = "atet" computes weights for average treatment effect on the treated
//        = "ateu" computes weights for average treatment effect on the untreated
real rowvector gmatch::ipw(string scalar est)
{
  real rowvector beta
  real colvector pscore, ipwwgt
  real matrix Ct
  this.reweight()
  if (args()<1) est="ate"
  Ct = this.Ct((this.varlist,"_cons"))
  beta   = this.logitbeta(this.T, this.X, this.W_orig, 1, Ct)
  // /* */ "propensity score (logit) model beta:"; beta
  pscore = this.logitpredict(this.X, beta)
  ipwwgt = this.logitweights(pscore, est)
  this.postbeta(beta)
  this.reweight(ipwwgt, pscore)
  return(beta)
}

// function that returns (weighted) mean of the dependent variable(s) in the control group
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
  if ((cols(X)==cols(beta)) & cols(beta)) {
    return(invlogit(X*beta'))
  }
  else if ((cols(X)==cols(beta)-1) & cols(beta)) {
    return(invlogit((X*beta[1..(cols(beta)-1)]') :+ beta[cols(beta)]))
  }
  else _error("X and beta are not conformable.")
}

// trims a generic column vector, x
// by default, trimming is at 1e-6 and 1-1e-6, which is useful for trimming propensity scores very close to 0 or 1
real colvector gmatch::trim(real colvector x, | real scalar minval, real scalar maxval)
{
  real colvector out
  if (args()<2) minval = 1e-6
  if (args()<3) maxval = 1-minval
  out = rowmax((J(rows(x),1,minval),rowmin((J(rows(x),1,maxval),x))))
  return(out)
}

// Function that turns a vector of pscores into IPW weights.
// Formulas match the normalized weights in Stata's teffects IPW command
//    pscore is a vector of propensity scores
//    est = "ate"  computes weights for average treatment effect (the default)
//        = "atet" computes weights for average treatment effect on the treated
//        = "ateu" computes weights for average treatment effect on the untreated
real colvector gmatch::logitweights(real colvector pscore, | string scalar est)
{
  real colvector pm
  real matrix ipwwgt
  if (args()<2) est="ate"

  if (any(pscore:<=0) | any(pscore:>=1)) _error("Propensity scores need to be greater than 0 and less than 1.")
//  /* */ if (minmax[1,1]<=0.03 & (strlower(est)=="ate" | strlower(est)=="ateu")) errprintf("Warning: minimum propensity score is %12.0g \n", minmax[1,1])
//  /* */ if (minmax[1,2]>=0.97 & (strlower(est)=="ate" | strlower(est)=="atet")) errprintf("Warning: maximum propensity score is %12.0g \n", minmax[1,2])

  pm = 1 :- (!this.T)
  if      (strlower(est)=="ate")   ipwwgt = (pm :/pscore) :+ (!pm:/(1:-pscore))
  else if (strlower(est)=="atet")  ipwwgt =  pm :+ (!pm :* (pscore:/(1:-pscore)))
  else if (strlower(est)=="ateu")  ipwwgt = !pm :+ ( pm :* ((1:-pscore):/pscore))
  else _error(est + " is an invalid argument for gmatch::logitweights()")

  // normalize the weights to have mean 1 in each group
  if (strlower(est)=="ate" | strlower(est)=="atet") ipwwgt[this.sel0] = ipwwgt[this.sel0] :/ mean(ipwwgt[this.sel0], this.W_orig[this.sel0])
  if (strlower(est)=="ate" | strlower(est)=="ateu") ipwwgt[this.sel1] = ipwwgt[this.sel1] :/ mean(ipwwgt[this.sel1], this.W_orig[this.sel1])
  return(ipwwgt)
}

// Define function to calculate coefficients for a logit regression model
// A constant term is added to the model and its coefficient is included in the vector of betas
// The program looks at Stata local mlopts for options related to controlling maximization
real rowvector gmatch::logitbeta(real colvector Ymat, real matrix Xmat, | real colvector Wmat, real scalar addconst, real matrix Ct)
{
  transmorphic S
  if (args()<4) addconst=1

  S=moptimize_init()
  moptimize_init_evaluator(S, &gmatch_logit_eval())
  moptimize_init_evaluatortype(S,"lf")
  moptimize_init_depvar(S,1,Ymat)
  moptimize_init_eq_indepvars(S,1,Xmat)
  if (!addconst) moptimize_init_eq_cons(S, 1, "off")
  if (args()>=3 & any(Wmat:!=1)) moptimize_init_weight(S, Wmat)
  moptimize_init_eq_colnames(S, 1, (J(1,cols(Xmat),"x") + strofreal((1..cols(Xmat)))))
  moptimize_init_vcetype(S, "robust")
  if (args()>=5) moptimize_init_constraints(S, Ct)
  if (st_local("mlopts")!="") moptimize_init_mlopts(S, st_local("mlopts"))

  moptimize(S)
  // /* */ "Logit model coefficients and robust standard errors:"; moptimize_result_display(S)
  return(moptimize_result_coefs(S))
}

// builds a constraint matrix for optimization commands
// largly based on Dave Drukker's post: https://blog.stata.com/2016/02/09/programming-an-estimation-command-in-stata-handling-factor-variables-in-optimize/
real matrix gmatch::Ct(string rowvector varlist) {
  string scalar tempmat
  real scalar ko, p, j
  real matrix Ct, mo
  tempmat = st_tempname()
  st_matrix(tempmat,J(1,length(varlist),0))
  stata("matrix colnames " +  tempmat + " = " + invtokens(varlist))
  stata("_ms_omit_info   " +  tempmat)
  mo = st_matrix("r(omit)")
  ko = sum(mo)
  p  = cols(mo)
  if (ko) {
    Ct   = J(0, p, .)
    for (j=1; j<=p; j++) {
      if (mo[j]==1) Ct = Ct \ e(j, p)
    }
    Ct = Ct, J(ko, 1, 0)
  }
  else Ct = J(0,p+1,.)
  return(Ct)
}

// sends coefficients and N back to Stata in a matrix named `gmatch_beta_out'
void gmatch::postbeta(real rowvector beta) {
  string scalar tempmatname
  st_eclear()
  tempmatname=st_tempname()
  st_matrix(tempmatname,beta)
  st_local("gmatch_beta_out",tempmatname)
  if      ((this.K==cols(beta)  ) & cols(beta)) st_matrixcolstripe(tempmatname, (J(cols(beta),1,""),this.varlist'))
  else if ((this.K==cols(beta)-1) & cols(beta)) st_matrixcolstripe(tempmatname, (J(cols(beta),1,""),(this.varlist' \ "_cons")))
  else                                                      _error("beta does not have the expected dimensions.")
  st_local("gmatch_N_out",strofreal(this.N))
}

void gmatch_logit_eval(transmorphic S, real rowvector beta, real colvector lnf)
{
  real colvector Y, pm, xb, lj
  Y  = moptimize_util_depvar(S, 1)
  xb = moptimize_util_xb(S, beta, 1)
  pm = 2*(Y :!= 0) :- 1
  lj = invlogit(pm:*xb)
  if (anyof(lj, 0)) {
    lnf = .
    return
  }
  lnf  = ln(lj)
}

// function that computes a variety of matching weights schemes, including CBPS weights (and returns them in this.W_mtch)
//    est corresponds to the options in gmatch::logitweights()
//        "ate"  computes weights for average treatment effect (the default)
//        "atet" computes weights for average treatment effect on the treated
//        "ateu" computes weights for average treatment effect on the untreated
//    fctn corresponds to the balance measure
//        "mean_sd_sq" minimizes the mean standardized difference squared
//    denominator is passed to stddiff() and related functions
//    oid=1 turns on the "over-identified" version of the CBPS model; oid=0 leaves it off
//    cvopt adds the CV of the matching weights to the optimization objective function
//         Let loss_0 be the objective function and wgt_cv() be the coefficient of variation of the matching weights
//         This function returns a vector, loss, to the optimization engine instead of loss_0.
//               (The elements of the vector will be summed because we're using a gf evaluator type.)
//         If cvopt=(a,b,c), then the loss function is modified as:
//              loss = ( loss_0 \ a * abs((wgt_cv() - b)^c) )
//         With >=6 arguments, cvopt=(a,b,c,d,e,f), the loss function also targets skewness of the weights (wgt_skewness())
//         Specifically, the loss function is modified as:
//              loss = ( loss_0 \ a * abs((wgt_cv() - b)^c) \ e * abs((wgt_skewness() - e)^f) )
//         With >=9 arguments, cvopt=(a,b,c,d,e,f,g,h,i), the loss function also targets kurtosis of the weights (wgt_kurtosis())
//         Specifically, the loss function is modified as:
//              loss = ( loss_0 \ a * abs((wgt_cv() - b)^c) \ e * abs((wgt_skewness() - e)^f) \ g * abs((wgt_kurtosis() - h)^i))
//    Other remarks: I'm using optimize() instead of _optimize(); because of this, errors that occuure
//                   during optimization will cause gmatch() to abort with error.
//                   In the case that convergance is not achieved, a warning will be printed to screen
//                   but the function will NOT abort with error. (This follows common practice for other Stata commands, like logit.)
real rowvector gmatch::gmatch(| string scalar est,
                                string scalar fctn,
                                real scalar denominator,
                                real rowvector cvopt)
{
  real rowvector beta
  real colvector pscore, cbpswgt
  real matrix ww, Ct
  real scalar oid, unnorm
  if (args()<1) est="ate"
  if (args()<2) fctn="sd_sq"
  if (args()<3) denominator=1
  if (args()<4) cvopt=J(1,0,.)
  this.reweight()
  oid = 0

  // If the user is asking for the IPW result, just call my ipw() function
  if (fctn=="ipw") {
    if (!length(cvopt)) return(this.ipw(est))
    else _error("IPW does not work with modified loss function")
  }

  // I have two implementations of the CBPS function.  Here I pick the one I need.
  // cbps_port_stata - has the gradient functions built in (so it converges faster)
  //                 - but it cannot deal with weighted data.
  //                 - was based on the Stata implementation of CBPS by Filip Premik
  // cbps_port_r     - works with weighted data
  //                 - doesn't have the gradient functions, and therefore
  //                      (1) works with cvopt and
  //                      (2) converges more slowly
  //                 - was based on the R implementation of CBPS on CRAN by Imai et al.
  // In addition, the program looks at Stata local mlopts with instructions for controlling maximization
  else if (fctn=="ipwcbps") {
    fctn="cbps"
    oid=1
  }
  if (fctn=="cbps" & all(this.W_orig:==1)) {
    fctn="cbps_port_stata"
  }
  else if (fctn=="cbps") {
    fctn="cbps_port_r"
  }

  // set up the optimization problem
  // see Mata documenation for optimize()
  transmorphic S
  S=optimize_init()
  optimize_init_evaluator(S, &gmatch_cbps_eval())
  optimize_init_which(S, "min")
  optimize_init_argument(S, 1, this)
  optimize_init_argument(S, 2, est)
  optimize_init_argument(S, 3, fctn)
  optimize_init_argument(S, 4, denominator)
  optimize_init_argument(S, 5, oid)
  optimize_init_argument(S, 6, cvopt)
  optimize_init_technique(S, "bfgs")     // The default is "nr". Since bfgs is faster, I'm going to use bfgs unless the user picked another technique in st_local("mlopts").
  optimize_init_tracelevel(S, "value" )  // "none", "value", "params"

  // the remaining optimization options depend on the method
  if (fctn=="cbps_port_r") {
    optimize_init_conv_ptol(S,  1e-13)
    optimize_init_conv_vtol(S,  1e-14)
    optimize_init_conv_nrtol(S, 1e-12)
    optimize_init_evaluatortype(S,"d0")
  }
  else if (fctn=="cbps_port_stata") {
    optimize_init_conv_ptol(S,  1e-13)
    optimize_init_conv_vtol(S,  1e-14)
    optimize_init_conv_nrtol(S, 1e-12)
    if (oid)  optimize_init_evaluatortype(S,"gf1")  // for overidentified version
    else      optimize_init_evaluatortype(S,"d1")   // d1 if I'm running plain vanilla. otherwise just use "do" (numerical gradient)
  }
  else if (fctn=="mean_sd_sq" | fctn=="sd_sq" | fctn=="stdprogdiff") {
    optimize_init_evaluatortype(S,"d0")
    optimize_init_conv_ignorenrtol(S, "off")
    optimize_init_conv_ptol(S,  1e-10)
    optimize_init_conv_vtol(S,  1e-11)
    optimize_init_conv_nrtol(S, 1e-9)
  }
  else _error(fctn + " is invalid with gmatch::gmatch()")

  // look for a local macro (created by the .ado program) with additonal maximization options, such as iterate()
  if (st_local("mlopts")!="") optimize_init_mlopts(S, st_local("mlopts"))

  // cvopt adds 1 or more elements to the loss function
  // Switch to a gf evaluator type, so the elements get summed
  // I don't have gradient functions programmed, so we use gf0
  if (length(cvopt)) optimize_init_evaluatortype(S,"gf0")

  // for certain methods,
  // -- normalize Xs to mean 0, sd 1, apply SVD
  // -- add a column with constant term
  if (fctn=="cbps_port_r") {
    real matrix sel, meansP_orig, sdP_orig, svd_s, svd_v, svd_s_inv, tmp
    unnorm = 1
    meansP_orig = mean(this.X, this.W)
    sdP_orig = sqrt(this.diagvariance(this.X, this.W))
    sel = selectindex(sdP_orig')' // if X has column of zeros (e.g,. the omited category for a factor variable), then the SD is 0 and Xstd would be ".".  Therefore, take that column out of the X matrix.
    this.Xstd = ((this.X[.,sel] :- meansP_orig[.,sel]) :/ sdP_orig[.,sel], J(this.N_raw,1,1))
    pragma unset svd_v
    pragma unset svd_s
    _svd(this.Xstd, svd_s, svd_v)
  }
  else if (fctn=="cbps_port_stata") {
    unnorm=0
    if (!length(this.XC)) this.XC = (this.X, J(this.N_raw,1,1)) // not the most efficient way to add a constant term to X -- data is copied from a view into a matrix -- but at least I only do it once
  }
  else unnorm=0

  "Step 1 (initial values from logit model):"
  real rowvector beta_logit
  if (fctn=="cbps_port_r") {
    // constant term was already added to Xstd above, so we don't add a constant again in logitbeta()
    Ct = this.Ct((this.varlist[sel],"_cons"))
    beta_logit = this.logitbeta(this.T, this.Xstd, this.W, 0, Ct)
  }
  else {
    // constant term will be added to X in last column by the logitbeta function
    Ct = this.Ct((this.varlist,"_cons"))
    beta_logit = this.logitbeta(this.T, this.X, this.W, 1, Ct)
  }
  optimize_init_constraints(S, Ct)
  optimize_init_params(S, beta_logit)

  // This is an extra matrix the can be passed to the optimization engine. I use it for different purposes.
  // It is only calculated once -- not once each time the objective function is called.
  if (fctn=="cbps_port_stata") {
    ww = this.cbps_port_stata_wgt_matrix(beta_logit, oid, est)
    ww = invsym(ww)
  }
  else if (fctn=="cbps_port_r" & !oid) {
    ww = invsym(quadcross(this.Xstd,this.W,this.Xstd))
  }
  else ww = .
  optimize_init_argument(S, 7, ww)
  ""

  "Step 2 (CBPS) :"
  (void) optimize(S)
  beta    = optimize_result_params(S)
  st_numscalar("r(converged)", optimize_result_converged(S))
  // /* */ "beta:" ; beta

  // undoing the normalization and SVD
  if (unnorm)  {
    this.Xstd = .
    svd_s_inv = svd_s:^-1
    svd_s_inv = svd_s_inv :* (svd_s :> 1e-5)
    beta = (svd_v' * diag(svd_s_inv) * beta')'
    if (length(beta)<this.K) { // deal with the columns I took out above
      tmp = J(1,this.K+1,0)
      tmp[(sel, this.K+1)] = beta
      beta = tmp
    }
    beta[sel] = (beta[sel] :/ sdP_orig[sel])
    beta[this.K+1] = beta[this.K+1] :- meansP_orig[sel] * beta[sel]'
    // /* */ "CBPS beta after undoing the normalization"; ((this.varlist,"_cons")', strofreal(beta)')
  }

  pscore = this.logitpredict(this.X, beta)
  pscore = this.trim(pscore)
  cbpswgt = this.logitweights(pscore, est)

  this.postbeta(beta)
  this.reweight(cbpswgt, pscore)
  return(beta)
}

// helper function -- note this is not a member of the class.
// It exists because optimize does not allow you to have evaluator functions that are
// members of a class (https://www.stata.com/statalist/archive/2010-10/msg01188.html)
void gmatch_cbps_eval(real todo, real beta,
                      class gmatch scalar M,
                      string est, string fctn, real denominator, real oid, real cvopt, real ww,
                      real lnf, real g, real H)
{
  M.cbpseval(todo,beta,est,fctn,denominator,oid,cvopt,ww,lnf,g,H)
}

// specify the function to be called by optimize() to evaluate f(p).
void gmatch::cbpseval( real   scalar    todo,
                       real   rowvector beta,
                       string scalar    est,
                       string scalar    fctn,
                       real   scalar    denominator,
                       real   scalar    oid,
                       real   rowvector cvopt,
                       real   matrix    ww,
                       real   matrix    lnf,
                       real   matrix    g,
                       real   matrix    H)
{
  real colvector  pscore, cbpswgt
  if      (fctn=="cbps_port_stata")  this.cbps_port_stata(todo,beta,est,oid,ww,lnf,g,H)
  else if (fctn=="cbps_port_r")      this.cbps_port_r(todo,beta,est,oid,ww,lnf,g,H)
  else if (fctn=="mean_sd_sq" | fctn=="sd_sq"  | fctn=="stdprogdiff") {
    pscore = this.logitpredict(this.X, beta)
    pscore = this.trim(pscore)
    cbpswgt = this.logitweights(pscore, est)
    this.reweight(cbpswgt)
    if      (fctn=="mean_sd_sq")      lnf = this.mean_sd_sq(denominator)
    else if (fctn=="sd_sq")           lnf = quadsum(this.sd_sq(denominator))
    else if (fctn=="stdprogdiff")     lnf = quadsum(this.stdprogdiff(denominator):^2)
    else                              _error(fctn + " is invalid with gmatch::cbpseval()")
  }

  // cvopt, a row vector, modifies the loss function as documented above
  if (!length(cvopt)) return
  else if (mod(length(cvopt),3)!=0 | length(cvopt)<3 | length(cvopt)>9) _error("cvopt() should have 0, 3, 6, or 9 elements")
  else if (todo>0) _error("cvopt is not compatable with todo>0 in gmatch::cbpseval()")
  if (fctn=="cbps_port_stata" | fctn=="cbps_port_r") {
    if (fctn=="cbps_port_r") pscore = this.logitpredict(this.Xstd, beta)
    else                     pscore = this.logitpredict(this.X, beta)
    pscore = this.trim(pscore)
    cbpswgt = this.logitweights(pscore, est)
    this.reweight(cbpswgt)
  }
  if (cvopt[1,1]) lnf = (lnf \ (cvopt[1,1]:*abs((this.wgt_cv(est):-cvopt[1,2]):^cvopt[1,3])))

  if (length(cvopt)<6) return
  if (cvopt[1,4]) lnf = (lnf \ (cvopt[1,4]:*abs((this.wgt_skewness(est):-cvopt[1,5]):^cvopt[1,6])))

  if (length(cvopt)<9) return
  if (cvopt[1,7]) lnf = (lnf \ (cvopt[1,7]:*abs((this.wgt_kurtosis(est):-cvopt[1,8]):^cvopt[1,9])))
}

// Calls CBPS model (not over-identified)
// This just calls gmatch() -- described above.
real rowvector gmatch::cbps(| string scalar est, real scalar denominator)
{
  if (args()<1) est="ate"
  if (args()<2) denominator=1
  return(gmatch(est, "cbps", denominator))
}

// Calls over-identified CBPS model
// This just calls gmatch() -- described above.
real rowvector gmatch::cbpsoid(| string scalar est, real scalar denominator)
{
  if (args()<1) est="ate"
  if (args()<2) denominator=1
  return(gmatch(est, "ipwcbps", denominator))
}

// Port of the objective function from the Stata verion of CBPS
void gmatch::cbps_port_stata( real   scalar    todo,
                              real   rowvector beta,
                              string scalar    est,
                              real   scalar    oid,
                              real   matrix    ww,
                              real   matrix    lnf,
                              real   matrix    g,
                              real   matrix    H)
{
   real colvector  pscore
   pscore = this.logitpredict(this.X, beta)
   pscore = this.trim(pscore)
   real matrix dpscore, gg, G
   dpscore = pscore:*(1:-pscore)
   gg = this.cbps_port_stata_moments(pscore, dpscore, oid, est)
   lnf = gg' * ww * gg
   if (todo==0) return
   G = this.cbps_port_stata_gradient(pscore, oid, est)
   g = (G' * ww * gg :* (2:*this.N))'
}

// Port of the moment function from the Stata verion of CBPS
real colvector gmatch::cbps_port_stata_moments(real colvector pscore, real matrix dpscore, real scalar overid, string scalar est)
{
  real colvector gg

  if (strlower(est)=="ate") {
      gg=quadcross(this.XC, (this.T-pscore):/pscore:/(1:-pscore)):/this.N_raw
  }
  else if (strlower(est)=="atet") {
      gg=quadcross(this.XC, (this.T-pscore):/(1:-pscore)):/this.N1_raw
  }
  else _error(est + " is invalid with gmatch::cbps_port_stata_moments()")

  if(overid) {
    gg = (quadcross(this.XC, dpscore:*(this.T-pscore):/pscore:/(1:-pscore)):/this.N_raw \ gg)
  }

  return(gg)
}


// Port of the gradient function from the Stata verion of CBPS
real matrix gmatch::cbps_port_stata_gradient(real colvector pscore, real scalar overid, string scalar est)
{
  real matrix G, dw
  if (strlower(est)=="ate") {
    G = -(this.XC:*((this.T:-pscore):^2):/pscore:/(1:-pscore))'this.XC
  }
  else if (strlower(est)=="atet") {
    dw=(pscore:*(this.T:-1)):/(1:-pscore):*(this.N_raw/this.N1_raw)
    G = quadcross(this.XC:*dw, this.XC)
  }
  if (overid) {
    G = ((-(this.XC:*(pscore:*(1:-pscore)))' this.XC) \ G)
  }
  G = G :/ this.N
  return(G)
}


// Port of the weighting matrix function from the Stata verion of CBPS
real matrix gmatch::cbps_port_stata_wgt_matrix(real rowvector beta, real scalar overid, string scalar est)
{
  real matrix ww
  real colvector pscore, dpscore
  pscore  = this.logitpredict(this.X, beta)
  pscore  = this.trim(pscore)
  if (!overid) {
    if (strlower(est)=="ate") {
      ww = quadcross(this.XC:/(pscore:*(1:-pscore)), this.XC)
    }
    else if (strlower(est)=="atet") {
      ww = quadcross(this.XC:*(pscore:/(1:-pscore)):*(this.N_raw/this.N1_raw):^2, this.XC)
    }
  }
  else {
    dpscore = pscore:*(1:-pscore)
    if (strlower(est)=="ate") {
      ww = ((quadcross(this.XC:*(dpscore:^2:/pscore:/(1:-pscore)),this.XC), // this seems inefficint. isn't  pscore:/(1:-pscore) = dpscore:^-1 ?
             quadcross(this.XC:*(dpscore:/pscore:/(1:-pscore)),this.XC)) \
            (quadcross(this.XC:*(dpscore:/pscore:/(1:-pscore)),this.XC),
             quadcross(this.XC:*(1:/pscore:/(1:-pscore)),this.XC)))
    }
    else if (strlower(est)=="atet") {
      ww = ((quadcross(this.XC:*(pscore:/(1:-pscore):*dpscore:^2:/pscore:^2),this.XC),
             quadcross(this.XC:*(pscore:/(1:-pscore):*dpscore:/pscore):*(this.N_raw/this.N1_raw),this.XC)) \
            (quadcross(this.XC:*(pscore:/(1:-pscore):*dpscore:/pscore):*(this.N_raw/this.N1_raw),this.XC),
             quadcross(this.XC:*(pscore:/(1:-pscore)):*((this.N_raw/this.N1_raw)^2),this.XC)))
    }
  }
  ww=ww:/this.N_raw
  return(ww)
}


// Port of the gmm.func()  function from CBPS.Binary.R and ATT.wt.func (version 0.17)
void gmatch::cbps_port_r(real   scalar    todo,
                         real   rowvector beta,
                         string scalar    est,
                         real   scalar    overid,
                         real   matrix    ww,
                         real   matrix    lnf,
                         real   matrix    g,
                         real   matrix    H)
{
  real colvector pscore, w_cbps
  pscore = this.logitpredict(this.Xstd, beta)
  pscore = this.trim(pscore)
  if (strlower(est)=="atet") {
     w_cbps = (this.N/this.N1) :* (this.T:-pscore) :/ (1:-pscore)
  }
  else if (strlower(est)=="ate") {
     w_cbps = (pscore:-1:+this.T):^-1
  }
  if (!overid) {
     w_cbps = 1/this.N :* w_cbps
     lnf = abs(quadcross(w_cbps, this.W_orig, this.Xstd) * ww * quadcross(this.Xstd, this.W_orig, w_cbps))
  }
  else {
    real colvector gbar, xx1, xx2, xx3
    real matrix V
    gbar = (quadcross(this.Xstd, this.W_orig, this.T:-pscore) \ quadcross(this.Xstd, this.W_orig, w_cbps)) :/ this.N
    if (strlower(est)=="atet") {
      xx1 = this.Xstd:*sqrt((1:-pscore):*pscore)
      xx2 = this.Xstd:*sqrt(pscore:/(1:-pscore))
      xx3 = this.Xstd:*sqrt(pscore)
      V =  (quadcross(xx1, this.W_orig, xx1), quadcross(xx3, this.W_orig, xx3) \
            quadcross(xx3, this.W_orig, xx3), quadcross(xx2, this.W_orig, xx2) :* (this.N:/this.N1_raw)) :/ this.N1_raw
    }
    else if (strlower(est)=="ate") {
      xx1 = this.Xstd:*sqrt((1:-pscore):*pscore)
      xx2 = this.Xstd:*((pscore:*(1:-pscore)):^-.5)
      xx3 = this.Xstd
      V = (quadcross(xx1, this.W_orig, xx1), quadcross(xx3, this.W_orig, xx3) \
           quadcross(xx3, this.W_orig, xx3), quadcross(xx2, this.W_orig, xx2)) :/ this.N
    }
    else _error(est + " is not allowed.")
    lnf = gbar' * invsym(V) * gbar
  }
}



// extract optimize_init_*() options parsed by Stata program -mlopts-
// I just copied moptimize_init_mlopts source code, then updated according to
// the optimize() help manual
void optimize_init_mlopts(transmorphic scalar M, string scalar mlopts)
{
        string scalar arg, arg1, tok
        transmorphic t, t1

        t = tokeninit(" ","","()")
        tokenset(t,mlopts)
        arg = tokenget(t)
        t1 = tokeninit("()")

        while (strlen(arg)) {
          arg1 = ""
          if (strmatch(arg,"trace"))               optimize_init_tracelevel(M, "value")
          else if (strmatch(arg,"gradient"))       optimize_init_tracelevel(M, "gradient")
          else if (strmatch(arg,"hessian"))        optimize_init_tracelevel(M, "hessian")
          else if (strmatch(arg,"showstep"))       optimize_init_tracelevel(M, "step")
          else if (strmatch(arg,"nonrtolerance"))  optimize_init_conv_ignorenrtol(M, "on")
          else if (strmatch(arg,"showtolerance"))  optimize_init_tracelevel(M, "tolerance")
          else if (strmatch(arg,"difficult"))      optimize_init_singularHmethod(M, "hybrid")
          else {
            arg1 = tokenget(t)
            tokenset(t1,arg1)
            tok = tokenget(t1)
            if (strmatch(arg,"technique"))         optimize_init_technique(M, tok)
            else if (strmatch(arg,"iterate"))      optimize_init_conv_maxiter(M, strtoreal(tok))
            else if (strmatch(arg,"tolerance"))    optimize_init_conv_ptol(M, strtoreal(tok))
            else if (strmatch(arg,"ltolerance"))   optimize_init_conv_vtol(M, strtoreal(tok))
            else if (strmatch(arg,"nrtolerance"))  optimize_init_conv_nrtol(M, strtoreal(tok))
            else arg = arg1
          }
          if (!strmatch(arg,arg1)) {
            arg = tokenget(t)
          }
        }
}

end
