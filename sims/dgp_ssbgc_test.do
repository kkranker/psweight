cd "C:\Users\kkranker\Documents\Stata\Ado\Devel\gmatch\sims"
clear all
cls
cap log close dgp_ssbgc_test
set linesize 180
log using "dgp_ssbgc_test.log", replace name(dgp_ssbgc_test)

// DGP is saved in a separate .ado file
which dgp_ssbgc

version 15.1
set tracedepth 1
set cformat %9.2fc
set pformat %5.3f
set sformat %7.3f
local tabstatopts stat(n mean sd skew cv min p5 p25 p50 p75 p95 max) format(%9.2fc)


//set trace on

// test DGP setup - 1 rep of each scenario
set type float
foreach SC in A B C D E F G {
  set seed 1
  di _n(3) "-----------------" _n "Scenario `SC'" _n "-----------------" _n

  // generate data
  dgp_ssbgc `SC', n(1000000)
  desc
  list in 1/10, sep(0)

  // descriptive stats
  summ, format
  tabstat y, by(a) `tabstatopts'
  psgraph, bin(30) treated(a) pscore(ps) nodraw title("Scenario `SC'") legend(off) name(SC`SC'_true)
  local true `true' SC`SC'_true

  // IPW
  logit a w1-w10, nolog
  predict ps_hat, pr
  label var ps_hat "Propensity score from logit model"
  psgraph, bin(30) treated(a) pscore(ps_hat) nodraw title("Scenario `SC'") legend(off) name(SC`SC'_logit)
  local logit `logit' SC`SC'_logit
  tabstat ps_hat, by(a)  `tabstatopts'

  gen W = cond(a, 1, ps_hat/(1-ps_hat))
  summ W if !a, mean
  replace W = W / r(mean) if !a
  label var W "Weight from IPW"
  tabstat W, by(a) `tabstatopts'
  hist W if !a, bin(30) title("Scenario `SC'") legend(off) nodraw name(SC`SC'_weights)
  local weights `weights' SC`SC'_weights

  tabstat w1,        by(a)  `tabstatopts'
  tabstat w1 [aw=W], by(a)  `tabstatopts'
  drop ps_hat W

}
graph combine `true' , name(pscore_true)
graph combine `logit', name(pscore_hat_logit)
graph combine `weights', name(weights_ipw) note("Comparison group only")

log close dgp_ssbgc_test
