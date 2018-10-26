cd "C:\Users\kkranker\Documents\Stata\Ado\Devel\gmatch"
clear all
cls

cap log close sims1
set linesize 180
log using "sims/sims1.log", replace name(sims1)

// DGP is saved in a separate .ado file
adopath ++ "./sims"
which dgp_ssbgc
which addstats

program define onerep
  version 15.1
  set cformat %9.2fc
  set pformat %5.3f
  set sformat %7.3f
  set type float
  set maxiter 50

  // this lets me keep the "parallal sim"  command simple
  // but probably slows things down
  cd "C:\Users\kkranker\Documents\Stata\Ado\Devel\gmatch"
  adopath ++ "./sims"

  syntax [anything], n(numlist min=1 >0 integer sort)
  tempname bbld

  local L : word count `n'
  dgp_ssbgc `anything', n(`: word `L' of `n'')

  local l = `L'
  while `l' > 0 {
    local thisN : word `l' of `n'
    di _n(3) "N = `thisN'" _n(2)
    keep in 1/`thisN'

    gmatch a w1-w10, atet ipw pooledvariance
    regress y i.a [aw=_weight], noheader
    addstats `bbld' ipw_`thisN' 1.a

    local --l
  }

  // gmatch a w1-w10, atet cbps pooledvariance
  ereturn clear
  ereturn post `bbld'
end


set seed 1

onerep A, n(100)
ereturn display

// simulate _b, reps(10): onerep A, n(500)
// codebook, compact

parallel setclusters `c(processors_mach)'
parallel sim, expr(_b) reps(1000) processors(1): onerep A, n(8000 4000 2000 1000 800 600 400 200 100)
codebook, compact

log close sims1
